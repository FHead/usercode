#ifndef CALCULATEME_H_11251_ASGKAJSGKDJASKDGJASKDJGASGD
#define CALCULATEME_H_11251_ASGKAJSGKDJASKDGJASKDJGASGD

// CalculateME, copy-paste expression from Roberto Mega-Vorales and compile into code
//
// GetZZBackground: ZZ_full_CM_JHU_725.txt
// GetZABackground: ZA_full_CM_JHU_726.txt
// GetZZZABackground: ZApZZ_tuInt_CM_JHU.txt
// GetZANoIntZZBackground: ZAnointZZ_tuchannels_fullJHU_CM_OnePiece802.txt

#include <complex>

#include "AngleConversion.h"

struct GeneralScalarParameters;
double GetM2Scalar(double HMass, double ZMass, double Z2Mass, double Phi0, double CosTheta0, double Phi, double CosTheta1, double CosTheta2);
double GetM2PseudoScalar(double HMass, double ZMass, double Z2Mass, double Phi0, double CosTheta0, double Phi, double CosTheta1, double CosTheta2);
double GetPhaseSpace(double HMass, double ZMass, double Z2Mass);
double GetZZBackground(EventParameters &event, bool UpType = true);
double GetZABackground(EventParameters &event, bool UpType = true);
double GetZZZABackground(EventParameters &event, bool Uptype = true);
double GetZANoIntZZBackground(EventParameters &event, bool UpType = true);
double GetZANoIntZZBackgroundOriginal(EventParameters &event, bool UpType = true);
double GetZABackground_ThetaPhi1(EventParameters &event, bool UpType = true);
double GetZABackground_Thetaphi(EventParameters &event, bool UpType = true);
double GetZABackground_ThetaTheta1(EventParameters &event, bool UpType = true);
double GetZABackground_ThetaTheta2(EventParameters &event, bool UpType = true);
double GetZZBackground_ThetaPhi1(EventParameters &event, bool UpType = true);
double GetZZBackground_Thetaphi(EventParameters &event, bool UpType = true);
double GetZZBackground_ThetaTheta1(EventParameters &event, bool UpType = true);
double GetZZBackground_ThetaTheta2(EventParameters &event, bool UpType = true);
std::complex<double> GetZZZAAABackgroundComplex(LeptonVectors &leptons, bool UpType = true);
double GetZZZAAABackground(EventParameters &event, bool UpType = true);
double GetZZZAAABackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetGeneralScalar(LeptonVectors &leptons, GeneralScalarParameters &parameters, bool Reverse = false, bool UpType = true);
double GetHyp1GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp3GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp4GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp5GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp6GeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetGoldenChannelMoneyComponents(LeptonVectors &leptons, string Component, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_s_ZZ.txt
double GetDiffCxn4VecSZZ(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_s_ZZpZA.txt
double GetDiffCxn4VecSZZpZA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_tandu_ZA.txt
double GetDiffCxn4VecTAndUZA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_tandu_ZZ.txt
double GetDiffCxn4VecTAndUZZ(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_sandtint_ZA.txt
double GetDiffCxn4VecSAndTIntZA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_tandu_ZAZZ.txt
double GetDiffCxn4VecTAndUIntZAZZ(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_sandtint_ZZ.txt
double GetDiffCxn4VecSAndTIntZZ(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_sandtint_ZZpZA.txt
double GetDiffCxn4VecSAndTIntZZpZA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_full_ZAAA.txt
double GetDiffCxn4VecFullZAAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_full_ZAZZ.txt
double GetDiffCxn4VecFullZAZZ(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_full_ZZAA.txt
double GetDiffCxn4VecFullZZAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_s_ZApZZpAA.txt
double GetDiffCxn4VecSZApZZpAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_sandtint_ZApZZpAA.txt
double GetDiffCxn4VecSAndTIntZApZZpAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_tandu_ZApZZpAA.txt
double GetDiffCxn4VecTAndUZApZZpAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

double GetHyp2BGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2CGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2DGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2EGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2FGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2GGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetHyp2HGeneralScalar(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_full_AA.txt
double GetDiffCxn4VecFullAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_s_AA.txt
double GetDiffCxn4VecSAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_sandtint_AA.txt
double GetDiffCxn4VecSAndTIntAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

// DiffCxn4Vec_tandu_AA.txt
double GetDiffCxn4VecTAndUAA(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);


double GetGeneralScalar_aAI_aAI(LeptonVectors &leptons, double aAI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aAR(LeptonVectors &leptons, double aAI = 1, double aAR = 1, bool UpType = false);
double GetGeneralScalar_aAI_aAdI(LeptonVectors &leptons, double aAI = 1, double aAdI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aAdR(LeptonVectors &leptons, double aAI = 1, double aAdR = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZAI(LeptonVectors &leptons, double aAI = 1, double aZAI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZAR(LeptonVectors &leptons, double aAI = 1, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZAdI(LeptonVectors &leptons, double aAI = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZAdR(LeptonVectors &leptons, double aAI = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZI(LeptonVectors &leptons, double aAI = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZR(LeptonVectors &leptons, double aAI = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZdI(LeptonVectors &leptons, double aAI = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aAI_aZdR(LeptonVectors &leptons, double aAI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aAI_ahI(LeptonVectors &leptons, double aAI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aAI_ahR(LeptonVectors &leptons, double aAI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aAR(LeptonVectors &leptons, double aAR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aAdI(LeptonVectors &leptons, double aAR = 1, double aAdI = 1, bool UpType = false);
double GetGeneralScalar_aAR_aAdR(LeptonVectors &leptons, double aAR = 1, double aAdR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZAI(LeptonVectors &leptons, double aAR = 1, double aZAI = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZAR(LeptonVectors &leptons, double aAR = 1, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZAdI(LeptonVectors &leptons, double aAR = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZAdR(LeptonVectors &leptons, double aAR = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZI(LeptonVectors &leptons, double aAR = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZR(LeptonVectors &leptons, double aAR = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZdI(LeptonVectors &leptons, double aAR = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aAR_aZdR(LeptonVectors &leptons, double aAR = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aAR_ahI(LeptonVectors &leptons, double aAR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aAR_ahR(LeptonVectors &leptons, double aAR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aAdI(LeptonVectors &leptons, double aAdI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aAdR(LeptonVectors &leptons, double aAdI = 1, double aAdR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZAI(LeptonVectors &leptons, double aAdI = 1, double aZAI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZAR(LeptonVectors &leptons, double aAdI = 1, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZAdI(LeptonVectors &leptons, double aAdI = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZAdR(LeptonVectors &leptons, double aAdI = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZI(LeptonVectors &leptons, double aAdI = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZR(LeptonVectors &leptons, double aAdI = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZdI(LeptonVectors &leptons, double aAdI = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_aZdR(LeptonVectors &leptons, double aAdI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aAdI_ahI(LeptonVectors &leptons, double aAdI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aAdI_ahR(LeptonVectors &leptons, double aAdI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aAdR(LeptonVectors &leptons, double aAdR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZAI(LeptonVectors &leptons, double aAdR = 1, double aZAI = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZAR(LeptonVectors &leptons, double aAdR = 1, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZAdI(LeptonVectors &leptons, double aAdR = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZAdR(LeptonVectors &leptons, double aAdR = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZI(LeptonVectors &leptons, double aAdR = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZR(LeptonVectors &leptons, double aAdR = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZdI(LeptonVectors &leptons, double aAdR = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aAdR_aZdR(LeptonVectors &leptons, double aAdR = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aAdR_ahI(LeptonVectors &leptons, double aAdR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aAdR_ahR(LeptonVectors &leptons, double aAdR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZAI(LeptonVectors &leptons, double aZAI = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZAR(LeptonVectors &leptons, double aZAI = 1, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZAdI(LeptonVectors &leptons, double aZAI = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZAdR(LeptonVectors &leptons, double aZAI = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZI(LeptonVectors &leptons, double aZAI = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZR(LeptonVectors &leptons, double aZAI = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZdI(LeptonVectors &leptons, double aZAI = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZAI_aZdR(LeptonVectors &leptons, double aZAI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZAI_ahI(LeptonVectors &leptons, double aZAI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZAI_ahR(LeptonVectors &leptons, double aZAI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZAR(LeptonVectors &leptons, double aZAR = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZAdI(LeptonVectors &leptons, double aZAR = 1, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZAdR(LeptonVectors &leptons, double aZAR = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZI(LeptonVectors &leptons, double aZAR = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZR(LeptonVectors &leptons, double aZAR = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZdI(LeptonVectors &leptons, double aZAR = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZAR_aZdR(LeptonVectors &leptons, double aZAR = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZAR_ahI(LeptonVectors &leptons, double aZAR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZAR_ahR(LeptonVectors &leptons, double aZAR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZAdI(LeptonVectors &leptons, double aZAdI = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZAdR(LeptonVectors &leptons, double aZAdI = 1, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZI(LeptonVectors &leptons, double aZAdI = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZR(LeptonVectors &leptons, double aZAdI = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZdI(LeptonVectors &leptons, double aZAdI = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_aZdR(LeptonVectors &leptons, double aZAdI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_ahI(LeptonVectors &leptons, double aZAdI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZAdI_ahR(LeptonVectors &leptons, double aZAdI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_aZAdR(LeptonVectors &leptons, double aZAdR = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_aZI(LeptonVectors &leptons, double aZAdR = 1, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_aZR(LeptonVectors &leptons, double aZAdR = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_aZdI(LeptonVectors &leptons, double aZAdR = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_aZdR(LeptonVectors &leptons, double aZAdR = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_ahI(LeptonVectors &leptons, double aZAdR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZAdR_ahR(LeptonVectors &leptons, double aZAdR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZI_aZI(LeptonVectors &leptons, double aZI = 1, bool UpType = false);
double GetGeneralScalar_aZI_aZR(LeptonVectors &leptons, double aZI = 1, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZI_aZdI(LeptonVectors &leptons, double aZI = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZI_aZdR(LeptonVectors &leptons, double aZI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZI_ahI(LeptonVectors &leptons, double aZI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZI_ahR(LeptonVectors &leptons, double aZI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZR_aZR(LeptonVectors &leptons, double aZR = 1, bool UpType = false);
double GetGeneralScalar_aZR_aZdI(LeptonVectors &leptons, double aZR = 1, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZR_aZdR(LeptonVectors &leptons, double aZR = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZR_ahI(LeptonVectors &leptons, double aZR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZR_ahR(LeptonVectors &leptons, double aZR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZdI_aZdI(LeptonVectors &leptons, double aZdI = 1, bool UpType = false);
double GetGeneralScalar_aZdI_aZdR(LeptonVectors &leptons, double aZdI = 1, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZdI_ahI(LeptonVectors &leptons, double aZdI = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZdI_ahR(LeptonVectors &leptons, double aZdI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_aZdR_aZdR(LeptonVectors &leptons, double aZdR = 1, bool UpType = false);
double GetGeneralScalar_aZdR_ahI(LeptonVectors &leptons, double aZdR = 1, double ahI = 1, bool UpType = false);
double GetGeneralScalar_aZdR_ahR(LeptonVectors &leptons, double aZdR = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_ahI_ahI(LeptonVectors &leptons, double ahI = 1, bool UpType = false);
double GetGeneralScalar_ahI_ahR(LeptonVectors &leptons, double ahI = 1, double ahR = 1, bool UpType = false);
double GetGeneralScalar_ahR_ahR(LeptonVectors &leptons, double ahR = 1, bool UpType = false);

double GetCompABackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetCompBBackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetCompCBackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetCompDBackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetCompEBackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);
double GetCompFBackground(LeptonVectors &leptons, bool Reverse = false, bool UpType = true);

struct GeneralScalarParameters
{
   double aAR, aAI;
   double aAdR, aAdI;
   double aZAR, aZAI;
   double aZAdR, aZAdI;
   double aZR, aZI;
   double aZdR, aZdI;
   double ahR, ahI;
};

#endif

