#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
using namespace std;

#include "RecoLocalCalo/HcalRecAlgos/interface/HBHEPulseShapeFlag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"

HBHEPulseShapeFlagSetter::HBHEPulseShapeFlagSetter()
{
}

HBHEPulseShapeFlagSetter::HBHEPulseShapeFlagSetter(vector<double> LinearThreshold, vector<double> LinearCut,
   vector<double> RMS8MaxThreshold, vector<double> RMS8MaxCut,
   vector<double> LeftSlopeThreshold, vector<double> LeftSlopeCut,
   vector<double> RightSlopeThreshold, vector<double> RightSlopeCut,
   vector<double> RightSlopeSmallThreshold, vector<double> RightSlopeSmallCut)
{
   for(int i = 0; i < (int)LinearThreshold.size() && i < (int)LinearCut.size(); i++)
      mLambdaLinearCut.push_back(pair<double, double>(LinearThreshold[i], LinearCut[i]));
   sort(mLambdaLinearCut.begin(), mLambdaLinearCut.end());

   for(int i = 0; i < (int)RMS8MaxThreshold.size() && i < (int)RMS8MaxCut.size(); i++)
      mLambdaRMS8MaxCut.push_back(pair<double, double>(RMS8MaxThreshold[i], RMS8MaxCut[i]));
   sort(mLambdaRMS8MaxCut.begin(), mLambdaRMS8MaxCut.end());

   for(int i = 0; i < (int)LeftSlopeThreshold.size() && i < (int)LeftSlopeCut.size(); i++)
      mLeftSlopeCut.push_back(pair<double, double>(LeftSlopeThreshold[i], LeftSlopeCut[i]));
   sort(mLeftSlopeCut.begin(), mLeftSlopeCut.end());

   for(int i = 0; i < (int)RightSlopeThreshold.size() && i < (int)RightSlopeCut.size(); i++)
      mRightSlopeCut.push_back(pair<double, double>(RightSlopeThreshold[i], RightSlopeCut[i]));
   sort(mRightSlopeCut.begin(), mRightSlopeCut.end());

   for(int i = 0; i < (int)RightSlopeSmallThreshold.size() && i < (int)RightSlopeSmallCut.size(); i++)
      mRightSlopeSmallCut.push_back(pair<double, double>(RightSlopeSmallThreshold[i], RightSlopeSmallCut[i]));
   sort(mRightSlopeSmallCut.begin(), mRightSlopeSmallCut.end());

   Initialize();
}

HBHEPulseShapeFlagSetter::~HBHEPulseShapeFlagSetter()
{
}

void HBHEPulseShapeFlagSetter::Clear()
{
}

void HBHEPulseShapeFlagSetter::SetPulseShapeFlags(HBHERecHit& hbhe, const HBHEDataFrame &digi,
   const HcalCoder &coder, const HcalCalibrations &calib)
{
   hbhe.setFlagField(0, HcalCaloFlagLabels::HBHEFlatNoise);
   hbhe.setFlagField(0, HcalCaloFlagLabels::HBHESpikeNoise);
   hbhe.setFlagField(0, HcalCaloFlagLabels::HBHETriangleNoise);

   CaloSamples Tool;
   coder.adc2fC(digi, Tool);

   double Charge[10] = {0};
   for(int i = 0; i < (int)digi.size(); i++)
      Charge[i] = Tool[i] - calib.pedestal(digi.sample(i).capid());

   double TotalCharge = 0;
   for(int i = 0; i < 10; i++)
      TotalCharge = TotalCharge + Charge[i];

   if(TotalCharge < 20)
      return;
   
   double NominalChi2 = PerformNominalFit(Charge);
   double LinearChi2 = PerformLinearFit(Charge);
   double RMS8Max = CalculateRMS8Max(Charge);
   TriangleFitResult TriangleResult = PerformTriangleFit(Charge);
   
   double TS4Left = Charge[4] / TriangleResult.LeftSlope;
   double TS4Right = Charge[4] / -TriangleResult.RightSlope;

   if(TS4Left > 1000 || TS4Left < -1000)
      TS4Left = 1000;
   if(TS4Right > 1000 || TS4Right < -1000)
      TS4Right = 1000;

   if(CheckPassFilter(TotalCharge, log(LinearChi2) - log(NominalChi2), mLambdaLinearCut, -1) == false)
      hbhe.setFlagField(1, HcalCaloFlagLabels::HBHEFlatNoise);
   
   if(CheckPassFilter(TotalCharge, log(RMS8Max) * 2 - log(NominalChi2), mLambdaRMS8MaxCut, -1) == false)
      hbhe.setFlagField(1, HcalCaloFlagLabels::HBHESpikeNoise);

   if(CheckPassFilter(Charge[4], TS4Left, mLeftSlopeCut, 1) == false)
      hbhe.setFlagField(1, HcalCaloFlagLabels::HBHETriangleNoise);
   if(CheckPassFilter(Charge[4], TS4Right, mRightSlopeCut, 1) == false)
      hbhe.setFlagField(1, HcalCaloFlagLabels::HBHETriangleNoise);
   if(CheckPassFilter(Charge[4], TS4Right, mRightSlopeSmallCut, -1) == false)
      hbhe.setFlagField(1, HcalCaloFlagLabels::HBHETriangleNoise);
}

void HBHEPulseShapeFlagSetter::Initialize()
{
   vector<double> PulseShape;

   HcalPulseShapes Shapes;
   HcalPulseShapes::Shape HPDShape = Shapes.hbShape();

   for(int i = 0; i < 200; i++)
      PulseShape.push_back(HPDShape.at(i));
   PulseShape.insert(PulseShape.begin(), 150, 0);

   CumulativeIdealPulse.clear();
   CumulativeIdealPulse.push_back(0);
   for(unsigned int i = 1; i < PulseShape.size(); i++)
      CumulativeIdealPulse.push_back(CumulativeIdealPulse[i-1] + PulseShape[i]);
}

TriangleFitResult HBHEPulseShapeFlagSetter::PerformTriangleFit(double Charge[10])
{
   TriangleFitResult result;
   result.Chi2 = 0;
   result.LeftSlope = 0;
   result.RightSlope = 0;

   // right side, starting from TS4
   double MinimumRightChi2 = 1000000;
   for(int iTS = 6; iTS <= 10; iTS++)   // iTS is the place where first TS center in flat line
   {
      // fit a straight line to the triangle part
      double Nominator = 0;
      double Denominator = 0;

      for(int i = 5; i < iTS; i++)
      {
         Nominator = Nominator + (i - 4) * (Charge[i] - Charge[4]);
         Denominator = Denominator + (i - 4) * (i - 4);
      }

      double BestSlope = Nominator / Denominator;
      if(BestSlope > 0)
         BestSlope = 0;

      // check if the slope is reasonable
      if(iTS != 10)
      {
         if(BestSlope > -Charge[4] / (iTS - 4))
            BestSlope = -Charge[4] / (iTS - 4);
         if(BestSlope < -Charge[4] / (iTS - 1 - 4))
            BestSlope = -Charge[4] / (iTS - 1 - 4);
      }
      else
      {
         if(BestSlope < -Charge[4] / (iTS - 1 - 4)) 
            BestSlope = -Charge[4] / (iTS - 1 - 4);

      }

      // calculate partial chi2
      double Chi2 = 0;
      for(int i = 5; i < iTS; i++)
         Chi2 = Chi2 + (Charge[4] - Charge[i] + (i - 4) * BestSlope)
            * (Charge[4] - Charge[i] + (i - 4) * BestSlope);
      for(int i = iTS; i < 10; i++)
         Chi2 = Chi2 + Charge[i] * Charge[i];

      if(MinimumRightChi2 > Chi2)
      {
         MinimumRightChi2 = Chi2;
         result.RightSlope = BestSlope;
      }
   }

   // left side
   double MinimumLeftChi2 = 1000000;
   for(int iTS = 0; iTS < 4; iTS++)   // the first time after linear fit ends
   {
      // fit a straight line to the triangle part
      double Nominator = 0;
      double Denominator = 0;
      for(int i = iTS; i < 4; i++)
      {
         Nominator = Nominator + (i - 4) * (Charge[i] - Charge[4]);
         Denominator = Denominator + (i - 4) * (i - 4);
      }

      double BestSlope = Nominator / Denominator;
      if(BestSlope < 0)
         BestSlope = 0;

      // check slope
      if(iTS != 0)
      {
         if(BestSlope > Charge[4] / (4 - iTS))
            BestSlope = Charge[4] / (4 - iTS);
         if(BestSlope < Charge[4] / (4 + 1 - iTS))
            BestSlope = Charge[4] / (4 + 1 - iTS);
      }
      else
      {
         if(BestSlope > Charge[4] / (4 - iTS))
            BestSlope = Charge[4] / (4 - iTS);
      }

      // calculate minimum chi2
      double Chi2 = 0;
      for(int i = 0; i < iTS; i++)
         Chi2 = Chi2 + Charge[i] * Charge[i];
      for(int i = iTS; i < 4; i++)
         Chi2 = Chi2 + (Charge[4] - Charge[i] + (i - 4) * BestSlope)
            * (Charge[4] - Charge[i] + (i - 4) * BestSlope);

      if(MinimumLeftChi2 > Chi2)
      {
         MinimumLeftChi2 = Chi2;
         result.LeftSlope = BestSlope;
      }
   }

   result.Chi2 = MinimumLeftChi2 + MinimumRightChi2;

   return result;
}

double HBHEPulseShapeFlagSetter::PerformNominalFit(double Charge[10])
{
   double MinimumChi2 = 100000;

   double F[10] = {0};
   double SumF2 = 0;
   double SumTF = 0;
   double SumT2 = 0;

   for(int i = 0; i + 250 < (int)CumulativeIdealPulse.size(); i++)
   {
      if(CumulativeIdealPulse[i+250] - CumulativeIdealPulse[i] < 1e-5)
         continue;

      for(int j = 0; j < 10; j++)
         F[j] = CumulativeIdealPulse[i+j*25+25] - CumulativeIdealPulse[i+j*25];

      SumF2 = 0;
      SumTF = 0;
      SumT2 = 0;
      for(int j = 0; j < 10; j++)
      {
         SumF2 = SumF2 + F[j] * F[j] / fabs(Charge[j]);
         SumTF = SumTF + F[j];
         SumT2 = SumT2 + fabs(Charge[j]);
      }

      double Chi2 = SumT2 - SumTF * SumTF / SumF2;

      if(Chi2 < MinimumChi2)
         MinimumChi2 = Chi2;
   }

   return MinimumChi2;
}

double HBHEPulseShapeFlagSetter::CalculateRMS8Max(double Charge[10])
{
   double TempCharge[10];
   for(int i = 0; i < 10; i++)
      TempCharge[i] = Charge[i];

   sort(TempCharge, TempCharge + 10);

   double Total = 0;
   double Total2 = 0;
   for(int i = 0; i < 8; i++)
   {
      Total = Total + TempCharge[i];
      Total2 = Total2 + TempCharge[i] * TempCharge[i];
   }

   double RMS = sqrt(Total2 - Total * Total / 8);

   return RMS / TempCharge[9];
}

double HBHEPulseShapeFlagSetter::PerformLinearFit(double Charge[10])
{
   double SumTS2OverTi = 0;
   double SumTSOverTi = 0;
   double SumOverTi = 0;
   double SumTS = 45;
   double Sum1 = 10;

   for(int i = 0; i < 10; i++)
   {
      double Error2 = Charge[i];
      if(Charge[i] < 1)
         Error2 = 1;

      SumTS2OverTi = SumTS2OverTi + i * i / Error2;
      SumTSOverTi = SumTSOverTi + i / Error2;
      SumOverTi = SumOverTi + 1 / Error2;
   }

   double CM1 = SumTS2OverTi;
   double CM2 = SumTSOverTi;
   double CD1 = SumTSOverTi;
   double CD2 = SumOverTi;
   double C1 = SumTS;
   double C2 = Sum1;

   double Slope = (C1 * CD2 - C2 * CD1) / (CM1 * CD2 - CM2 * CD1);
   double Intercept = (C1 * CM2 - C2 * CM1) / (CD1 * CM2 - CD2 * CM1);

   double Chi2 = 0;
   for(int i = 0; i < 10; i++)
   {
      double Deviation = Slope * i + Intercept - Charge[i];
      double Error2 = Charge[i];
      if(Charge[i] < 1)
         Error2 = 1;
      Chi2 = Chi2 + Deviation * Deviation / Error2;
   }

   return Chi2;
}

bool HBHEPulseShapeFlagSetter::CheckPassFilter(double Charge, double Discriminant,
   vector<pair<double, double> > &Cuts, int Side)
{
   if(Charge <= Cuts[0].first)   // too small to cut on
      return true;

   int IndexLargerThanCharge = -1;
   for(int i = 1; i < (int)Cuts.size(); i++)
   {
      if(Cuts[i].first > Charge)
      {
         IndexLargerThanCharge = i;
         break;
      }
   }

   double limit = 10000;

   if(IndexLargerThanCharge == -1)
      limit = Cuts[Cuts.size()-1].second;
   else
   {
      double C1 = Cuts[IndexLargerThanCharge].first;
      double C2 = Cuts[IndexLargerThanCharge-1].first;
      double L1 = Cuts[IndexLargerThanCharge].second;
      double L2 = Cuts[IndexLargerThanCharge-1].second;

      limit = (Charge - C1) / (C2 - C1) * (L2 - L1) + L1;
   }

   if(Side > 0 && Discriminant > limit)
      return false;
   if(Side < 0 && Discriminant < limit)
      return false;

   return true;
}


