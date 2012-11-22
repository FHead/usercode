#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

#include "TauHelperFunctions2.h"
#include "DrawRandom2.h"
#include "ProgressBar.h"

#include "ElectronEfficiencyMap.h"
#include "MuonEfficiencyMap.h"
#include "CalculateME.h"
#include "Constants.h"

#define TOTALMASSBIN 50
#define SuffixCount 7

int main(int argc, char *argv[]);
void SmearLepton(FourVector &Vector);
inline unsigned int LeptonResponseMapFindBin( const double value, const double bins[], unsigned int nbins);

int main(int argc, char *argv[])
{
   ////////////////////////////////////////
   // HEADER SECTION: DEFINES PARAMETERS //
   ////////////////////////////////////////

   double HiggsMass = 125;   // Will be overwritten by program arguments
   double Phi0Center = 1.2;   // Will be overwritten by program arguments
   double Theta0Center = 0;   // Will be overwritten by program arguments

   // Lepton efficiency functions:
   //    Efficiency1 is for Z1, Efficiency2 is for Z2
   //    Switch to GetElectronEfficiencyPtEtaPhi for electron efficiencies
   //       and GetMuonEfficiencyPtEtaPhi for muon efficiencies;
   double (*Efficiency1)(double, double, double) = &GetMuonEfficiencyPtEtaPhi;
   double (*Efficiency2)(double, double, double) = &GetMuonEfficiencyPtEtaPhi;
 
   int NBinsPhi0 = 12;
   int NBinsTheta0 = 18;
   int NBinsAngle = 12;
   int NBinsPhi = 10;
   int NBinsMass1 = 20;
   int NBinsMass1Written = 22;
   int NBinsMass2 = 20;
   int NBinsMass2Written = 22;
   int NTriesPhiOffset = 10;
   
   double PhiMin = 0, Theta1Min = -1, Theta2Min = -1, Phi0Min = 0, Theta0Min = -1;
   double PhiMax = 2 * PI, Theta1Max = 1, Theta2Max = 1, Phi0Max = 2 * PI, Theta0Max = 1;
   double Mass1Min = 0, Mass2Min = 0;
   double Mass1Max = HiggsMass, Mass2Max = HiggsMass / 2;
   
   srand(time(NULL));
   RandomBase Random;

   ////////////////////////////////////////////////
   // GETTING INPUT PARAMETERS FROM COMMAND LINE //
   ////////////////////////////////////////////////

   if(argc != 6)
   {
      cerr << "Usage: " << argv[0] << " HiggsMass Phi0 cos(Theta0) Tag SuppressHeader=(N|Y)" << endl;
      return -1;
   }

   HiggsMass = atof(argv[1]);
   Phi0Center = atof(argv[2]);
   Theta0Center = atof(argv[3]);
   
   if(HiggsMass < 10)
   {
      cerr << "Weird Higgs mass: " << HiggsMass << endl;
      return -1;
   }
   if(Theta0Center < -1 || Theta0Center> 1)
   {
      cerr << "Weird theta0 value." << endl;
      return -1;
   }
 
   cerr << "Starting job with Higgs mass = " << HiggsMass << endl;
   cerr << "Starting job with phi0 = " << Phi0Center << endl;
   cerr << "Starting job with cos(theta0) = " << Theta0Center << endl;

   /////////////////////////////
   // READ HIGGS/M4L SPECTRUM //
   /////////////////////////////

   vector<double> HiggsPTList;
   vector<double> HiggsEtaList;
   vector<double> HiggsWeightList;
   ifstream in((string("HiggsFiles/Background/HiggsFile_") + argv[1]).c_str());
   while(in)
   {
      double temp1 = -1, temp2 = 9999, temp3 = 0;
      in >> temp1 >> temp2 >> temp3;
      if(temp2 > 20 || temp1 < 0 || temp3 == 0)
         continue;

      HiggsPTList.push_back(temp1);
      HiggsEtaList.push_back(temp2);
      HiggsWeightList.push_back(temp3);
   }
   in.close();

   //////////////////////////
   // PREPARE OUTPUT FILES //
   //////////////////////////

   string suffixes[] = {"EffMEE", "Eff", "EfFMEA", "EffMEB", "EffMEC", "EffMED", "MEE"};
   vector<string> Suffixes(suffixes, suffixes + SuffixCount);

   string OutputLocation = "/wntmp/yichen";   // Change this to something!

   vector<ofstream *> out;
   for(int i = 0; i < (int)Suffixes.size(); i++)
   {
      string FileName = OutputLocation +"/IndividualMap_" + argv[4] + "_FixAngle_" + Suffixes[i] + ".map7";
      out.push_back(new ofstream(FileName.c_str()));
   }

   /////////////////////////////////
   // PREPARE AUXILIARY VARIABLES //
   /////////////////////////////////

   double PhiOffsets[1000] = {0};
   for(int iPhiOffset = 0; iPhiOffset < NTriesPhiOffset; iPhiOffset++)
      PhiOffsets[iPhiOffset] = (iPhiOffset + 0.25) * 2 * PI / NTriesPhiOffset;

   double HMass = HiggsMass;

   double Z1MassMin = Mass1Min;
   double Z1MassStep = (Mass1Max - Z1MassMin) / NBinsMass1;

   double Z2MassMin = Mass2Min;
   double Z2MassStep = (Mass2Max - Z2MassMin) / NBinsMass2;

   double Phi0Step = (Phi0Max - Phi0Min) / NBinsPhi0;
   double Theta0Step = (Theta0Max - Theta0Min) / NBinsTheta0;

   double Z1Mass[TOTALMASSBIN] = {0};
   double Z2Mass[TOTALMASSBIN] = {0};

   for(int iMass1 = 0; iMass1 < NBinsMass1; iMass1++)
   {
      for(int iMass2 = 0; iMass2 < NBinsMass2; iMass2++)
      {
         Z1Mass[iMass1] = (iMass1 + 0.5) * Z1MassStep + Z1MassMin;
         Z2Mass[iMass2] = (iMass2 + 0.5) * Z2MassStep + Z2MassMin;
      }
   }

   /////////////////////////////
   // WRITE HEADERS IF NEEDED //
   /////////////////////////////

   if(argv[5][0] != 'Y' && argv[5][0] != 'y')
   {
      for(int i = 0; i < (int)Suffixes.size(); i++)
      {
         *(out[i]) << HMass << endl;
         *(out[i]) << NBinsPhi0 << " " << NBinsTheta0 << " " << NBinsPhi
            << " " << NBinsAngle << " " << NBinsAngle << " " << NBinsMass1Written
            << " " << NBinsMass2Written << endl;
         *(out[i]) << "0 " << 2 * PI << endl;
         *(out[i]) << "-1 1" << endl;
         *(out[i]) << "0 " << 2 * PI << endl;
         *(out[i]) << "-1 1" << endl;
         *(out[i]) << "-1 1" << endl;
         *(out[i]) << Z1MassMin << " "
            << (HMass - Z1MassMin) / NBinsMass1 * NBinsMass1Written + Z1MassMin << endl;
         *(out[i]) << Z2MassMin << " "
            << (HMass / 2 - Z2MassMin) / NBinsMass2 * NBinsMass2Written + Z2MassMin << endl;
      }
   }

   ////////////////////
   // START THE LOOP //
   ////////////////////

   ProgressBar Bar(cout, NBinsPhi * NBinsAngle);
   Bar.SetStyle(1);

   for(int iPhi = 0; iPhi < NBinsPhi; iPhi++)
   {
      for(int iTheta1 = 0; iTheta1 < NBinsAngle; iTheta1++)
      {
         Bar.Update(iPhi * NBinsAngle + iTheta1 + 1);
         Bar.Print();

         for(int iTheta2 = 0; iTheta2 < NBinsAngle; iTheta2++)
         {
            // Prepare mass grid - since we need smearing and things will migrate from bin to bin
            double PassedTries[SuffixCount][TOTALMASSBIN][TOTALMASSBIN] = {0};
            for(int iS = 0; iS < SuffixCount; iS++)
            {
               for(int iMass1 = 0; iMass1 < NBinsMass1Written; iMass1++)
               {
                  for(int iMass2 = 0; iMass2 < NBinsMass2Written; iMass2++)
                     PassedTries[iS][iMass1][iMass2] = 0;
               }
            }

            for(int iMass1 = 0; iMass1 < NBinsMass1; iMass1++)
            {
               for(int iMass2 = 0; iMass2 < NBinsMass2; iMass2++)
               {
                  for(int iPhiOffset = 0; iPhiOffset < NTriesPhiOffset; iPhiOffset++)
                  {
                     // For each try, randomize parameters within the box
                     double Phi0 = Phi0Center + Random.DrawRandom(-0.5, 0.5) * Phi0Step;
                     double Theta0 = acos(Theta0Center + Random.DrawRandom(-0.499, 0.499) * Theta0Step);
                     double Phi = (iPhi + Random.DrawRandom(0, 1)) * 2 * PI / NBinsPhi;
                     double Theta1 = acos((iTheta1 + Random.DrawRandom(0, 1)) * 2 / NBinsAngle - 1);
                     double Theta2 = acos((iTheta2 + Random.DrawRandom(0, 1)) * 2 / NBinsAngle - 1);
                     double Z1 = Z1Mass[iMass1] + Random.DrawRandom(-0.5, 0.5) * Z1MassStep;
                     double Z2 = Z2Mass[iMass2] + Random.DrawRandom(-0.5, 0.5) * Z2MassStep;

                     if(Z1 + Z2 > HMass)   // don't bother if gen-level masses sum up above H mass
                        continue;
                     if(Z1 < Z2)   // avoid partial double-counting
                        continue;

                     if(iMass1 <= 1 || iMass2 <= 1)
                        continue;

                     // Prepare input parameters - event2 is for symmetrization
                     EventParameters event;                 EventParameters event2;
                     event.Phi0 = Phi0;                     event2.Phi0 = Phi0 + PI;
                     event.Theta0 = Theta0;                 event2.Theta0 = PI - Theta0;
                     event.Phi = Phi;                       event2.Phi = Phi;
                     event.Theta1 = Theta1;                 event2.Theta1 = Theta1;
                     event.Theta2 = Theta2;                 event2.Theta2 = Theta2;
                     event.HMass = HMass;                   event2.HMass = HMass;
                     event.ZMass = Z1;                      event2.ZMass = Z1;
                     event.Z2Mass = Z2;                     event2.Z2Mass = Z2;
                     event.PhiH = PhiOffsets[iPhiOffset];   event2.PhiH = PhiOffsets[iPhiOffset];
                     
                     // Get m4l PT and eta from our spectrum
                     int RandomIndex = (int)(Random.DrawRandom(0, 1) * HiggsPTList.size());

                     double HiggsEta = HiggsEtaList[RandomIndex];
                     double HiggsPT = HiggsPTList[RandomIndex];
                     double Weight = HiggsWeightList[RandomIndex];

                     // Calculate lepton vectors - one for plugging in expression, one for efficiency etc.
                     LeptonVectors Leptons0 = ConvertAnglesToVectorsRoberto(event, 0, 0);
                     LeptonVectors Leptons1 = ConvertAnglesToVectorsRoberto(event, HiggsPT, HiggsEta);

                     FourVector &FinalLepton11 = Leptons1.Lepton11;
                     FourVector &FinalLepton12 = Leptons1.Lepton12;
                     FourVector &FinalLepton21 = Leptons1.Lepton21;
                     FourVector &FinalLepton22 = Leptons1.Lepton22;

                     // Smear energy scale!  (and recalculate Z masses to see what bin we're in)
                     SmearLepton(FinalLepton11);
                     SmearLepton(FinalLepton12);
                     SmearLepton(FinalLepton21);
                     SmearLepton(FinalLepton22);

                     // Apply some basic lepton acceptance - this needs to be updated
                     if(FinalLepton11.GetAbsEta() > 2.5 || FinalLepton12.GetAbsEta() > 2.5
                           || FinalLepton21.GetAbsEta() > 2.5 || FinalLepton22.GetAbsEta() > 2.5)
                        continue;
                     if(FinalLepton11.GetPT() < 5 || FinalLepton12.GetPT() < 5
                           || FinalLepton21.GetPT() < 5 || FinalLepton22.GetPT() < 5)
                        continue;

                     // New masses after smearing
                     double NewZ1Mass = (FinalLepton11 + FinalLepton12).GetMass();
                     double NewZ2Mass = (FinalLepton21 + FinalLepton22).GetMass();

                     if(NewZ1Mass < NewZ2Mass)
                        swap(NewZ1Mass, NewZ2Mass);
                     if(NewZ1Mass < 40 || NewZ1Mass > 120)
                        continue;
                     if(NewZ2Mass < 12 || NewZ2Mass > 120)
                        continue;

                     int iZ1MassNew = (int)((NewZ1Mass - Z1MassMin) / Z1MassStep);
                     int iZ2MassNew = (int)((NewZ2Mass - Z2MassMin) / Z2MassStep);

                     if(iZ1MassNew >= TOTALMASSBIN || iZ1MassNew < 0)
                        continue;
                     if(iZ2MassNew >= TOTALMASSBIN || iZ2MassNew < 0)
                        continue;

                     // Finally, calculate lepton efficiencies
                     double LeptonEfficiency =
                        Efficiency1(FinalLepton11.GetPT(), FinalLepton11.GetAbsEta(), FinalLepton11.GetPhi())
                        * Efficiency1(FinalLepton12.GetPT(), FinalLepton12.GetAbsEta(), FinalLepton12.GetPhi())
                        * Efficiency2(FinalLepton21.GetPT(), FinalLepton21.GetAbsEta(), FinalLepton21.GetPhi())
                        * Efficiency2(FinalLepton22.GetPT(), FinalLepton22.GetAbsEta(), FinalLepton22.GetPhi());

                     // Calculate ME from expressions from Roberto
                     double MEWeightA = GetZANoIntZZBackground(event) + GetZANoIntZZBackground(event2);
                     double MEWeightB = GetZZZABackground(event) + GetZZZABackground(event2);
                     double MEWeightC = GetZZBackground(event) + GetZZBackground(event2);
                     double MEWeightD = GetZABackground(event) + GetZABackground(event2);
                     double MEWeightE =
                        GetZZZAAABackground(Leptons0, false, false) + GetZZZAAABackground(Leptons0, true, false)
                        + GetZZZAAABackground(Leptons0, false, true) + GetZZZAAABackground(Leptons0, true, true);

                     MEWeightA = MEWeightA * Z1 * Z2 * 4;   // The notorious Jacobian factor
                     MEWeightB = MEWeightB * Z1 * Z2 * 4;
                     MEWeightC = MEWeightC * Z1 * Z2 * 4;
                     MEWeightD = MEWeightD * Z1 * Z2 * 4;
                     MEWeightE = MEWeightE * Z1 * Z2 * 4;

                     // Put things in the correct mass bins
                     PassedTries[0][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight * MEWeightE;
                     PassedTries[1][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight;
                     PassedTries[2][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight * MEWeightA;
                     PassedTries[3][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight * MEWeightB;
                     PassedTries[4][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight * MEWeightC;
                     PassedTries[5][iZ1MassNew][iZ2MassNew] += LeptonEfficiency * Weight * MEWeightD;
                     PassedTries[6][iMass1][iMass2] += MEWeightE;
                  }   // end-of-phioffset loop
               }   // end-of-Z2Mass loop
            }   // end-of-Z1Mass loop

            // Write things out to file
            for(int iS = 0; iS < SuffixCount; iS++)
            {
               for(int iMass1 = 0; iMass1 < NBinsMass1Written; iMass1++)
               {
                  for(int iMass2 = 0; iMass2 < NBinsMass2Written; iMass2++)
                     *(out[iS]) << " " << PassedTries[iS][iMass1][iMass2];
                  *(out[iS]) << endl;
               }
            }
         }   // end-of-theta2 loop
      }   // end-of-theta1 loop
   }   // end-of-phi loop

   Bar.PrintLine();

   /////////////
   // CLEANUP //
   /////////////

   for(int iS = 0; iS < SuffixCount; iS++)
   {
      *(out[iS]) << endl;
      out[iS]->close();
      delete out[iS];
   }

   return 0;
}

void SmearLepton(FourVector &Vector)
{
   static RandomBase Random;

   // double PT = Vector.GetPT();
   // Vector = Vector * GenerateLeptonPtFromHist(LeptonResponseFile, "Electrons", PT, 0.3) / PT;

   static const double ptBins[15] = {5,7,8,9,10,12,14,16,18,20,25,30,35,40,50};
   // static const double etaBins[17] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442,1.566,1.8,2,2.1,2.2,2.3,2.4,2.5,2.6};
   static const double pzpBins[17] =
      {0,0.1973753202,0.3799489623,0.537049567,0.6640367703,0.761594156,0.833654607,0.894540039,0.9163871674,0.9468060128,0.9640275801,0.9704519366,0.97574313,0.9800963963,0.9836748577,0.9866142982,0.9890274022};

   int PtBin = LeptonResponseMapFindBin(Vector.GetPT(), ptBins, 14 - 1);
   int EtaBin = LeptonResponseMapFindBin(fabs(Vector[3] / Vector.GetP()), pzpBins, 16 - 1);

   static const double mean[16][14] = {
      {2.44991e-05, -0.000149334, 5.57965e-06, -0.000455782, -0.000239303, -0.000230198, -0.000109128, -0.00045267, -0.000177629, -9.83564e-05, -0.000429434, -0.000414415, -0.000375942, -0.000551152},
      {-0.000235641, -0.000291944, -0.000338815, -0.000330302, -0.000227168, -6.56378e-05, -0.000454549, -0.000429626, -0.000398811, -0.000264356, -0.000358339, -0.000510966, -0.000503078, -0.000629049},
      {-0.000125317, -0.000392274, -0.000280742, -0.000341562, -0.000252906, -0.000152172, -0.000285231, -0.00060248, -0.0004333, -0.000385277, -0.000599954, -0.000833558, -0.000498149, -0.000764647},
      {8.42627e-06, -0.000216617, -0.00031635, 9.67597e-05, -0.000184272, -0.000434312, -0.000504495, -9.0329e-05, -0.00046031, -0.000590764, -0.000497475, -0.000630796, -0.000785809, -0.000692762},
      {0.000398467, 3.28056e-05, 0.000113506, 0.000103844, -0.000346381, -7.89942e-05, -0.000347648, -0.000526769, -0.000685861, -0.000751779, -0.000864401, -0.000912099, -0.000752909, -0.00115101},
      {0.000674753, -0.00055251, -0.000221459, -0.000241277, -0.00100144, -0.000545017, -0.000906392, -0.000426533, -0.000820143, -0.00114611, -0.00131188, -0.00128459, -0.00147968, -0.00120705},
      {-0.00114604, -0.00138359, -0.00161128, -0.00140013, -0.00151796, -0.00161084, -0.0013402, -0.00153765, -0.0018934, -0.0017029, -0.00191201, -0.00180248, -0.00173663, -0.00187499},
      {-0.00109943, -0.00183758, -0.0018783, -0.00178843, -0.00192326, -0.00123929, -0.0018844, -0.00146853, -0.00133279, -0.00215737, -0.00188658, -0.00169481, -0.00223034, -0.00213486},
      {-0.00081687, -0.00135341, -0.00179217, -0.0013042, -0.00155819, -0.00163222, -0.00149117, -0.0019851, -0.00183379, -0.0018186, -0.00200237, -0.00212033, -0.00227697, -0.00246566},
      {-0.00126187, -0.000879105, -0.00142002, -0.00161185, -0.00175276, -0.00177354, -0.00159942, -0.00147886, -0.00154443, -0.00226977, -0.00233344, -0.00238358, -0.00278203, -0.00295222},
      {-0.00105092, -0.00198809, -0.000879293, -0.00240821, -0.0020562, -0.00226411, -0.00257947, -0.00259519, -0.00217108, -0.00309811, -0.00311022, -0.00161709, -0.00350283, -0.00378001},
      {-0.00150355, -0.00177692, -0.00186188, -0.00262448, -0.000894784, -0.00109178, -0.00200604, -0.00234174, -0.00210333, -0.00286123, -0.00257985, -0.00258507, -0.00222163, -0.00377099},
      {-0.000958727, -0.00253823, -0.00245147, -0.00242494, -0.00234165, -0.00233876, -0.00318113, -0.00128225, -0.00265225, -0.00295341, -0.00377778, -0.00363138, -0.0037289, -0.00438474},
      {-0.00109475, -0.00232223, -0.00240169, -0.00260532, -0.00166665, -0.00374696, -0.00297557, -0.00196759, -0.00371736, -0.00292813, -0.0040801, -0.00453383, -0.00625798, -0.00570589},
      {-0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006},
      {-0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006, -0.006}};

   static const double sigma[16][14] = {
      {0.00739344, 0.00775915, 0.0076801, 0.00825316, 0.00805092, 0.00840312, 0.0084091, 0.00878305, 0.00884362, 0.00913437, 0.00974197, 0.0101187, 0.0105767, 0.0110027},
      {0.00846152, 0.00848204, 0.00918793, 0.00914101, 0.00931898, 0.009, 0.00973429, 0.00957909, 0.0100585, 0.0102374, 0.0105348, 0.011008, 0.0113803, 0.0118701},
      {0.00936197, 0.00981436, 0.00978253, 0.00986696, 0.0100982, 0.0106991, 0.0103439, 0.0107398, 0.0105, 0.0111685, 0.011218, 0.0119209, 0.011927, 0.012761},
      {0.0100801, 0.0102789, 0.010672, 0.0106581, 0.0106458, 0.0112343, 0.0114141, 0.0105, 0.0117675, 0.0119833, 0.0120044, 0.0126498, 0.01298, 0.0136238},
      {0.0119224, 0.012466, 0.0122609, 0.0115, 0.0129602, 0.0134881, 0.0137757, 0.0136903, 0.0143067, 0.0143837, 0.0146947, 0.0151286, 0.0145, 0.0159705},
      {0.0150203, 0.0154148, 0.0157124, 0.0156906, 0.0159099, 0.0159859, 0.0162745, 0.0155, 0.0163747, 0.0170136, 0.0175251, 0.0176152, 0.0177868, 0.0184926},
      {0.0168999, 0.016, 0.0166911, 0.0168824, 0.0174103, 0.017854, 0.0180717, 0.018, 0.0187352, 0.018, 0.0192316, 0.0193895, 0.0194329, 0.0200056},
      {0.0164716, 0.0159999, 0.017239, 0.016, 0.0171992, 0.0171425, 0.0180467, 0.0179694, 0.0176833, 0.0183, 0.0186773, 0.0194139, 0.0198407, 0.0194},
      {0.015755, 0.0149999, 0.015, 0.0163016, 0.0163733, 0.0166871, 0.0168273, 0.0168261, 0.0170199, 0.0176337, 0.0180661, 0.0182993, 0.018, 0.0198938},
      {0.0181529, 0.0188446, 0.0185, 0.0186616, 0.019104, 0.0195675, 0.019, 0.0203664, 0.0202659, 0.0207058, 0.0216973, 0.0223438, 0.0235266, 0.0246652},
      {0.0195, 0.0212849, 0.02, 0.0228697, 0.0223192, 0.0226252, 0.0239364, 0.022749, 0.0231552, 0.0241103, 0.0253495, 0.0268228, 0.0273987, 0.0284956},
      {0.021, 0.0235681, 0.02, 0.0244535, 0.0245681, 0.0255102, 0.0259765, 0.0270135, 0.0259085, 0.0270402, 0.0287394, 0.0279296, 0.0319751, 0.0335993},
      {0.023, 0.024, 0.0265694, 0.0263221, 0.026445, 0.0273819, 0.0298332, 0.0296533, 0.0289494, 0.0313245, 0.033501, 0.036245, 0.0373117, 0.0407096},
      {0.0239199, 0.024, 0.0273773, 0.0278265, 0.0271103, 0.0280794, 0.0284486, 0.031538, 0.0289243, 0.0353967, 0.0372429, 0.0400654, 0.0440607, 0.0495239},
      {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},
      {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};

   static const double alphaL[16][14] = {
      {1.84171, 1.74128, 1.81004, 2.11071, 1.77622, 1.75484, 1.67844, 1.94784, 1.68338, 1.62709, 1.76951, 1.70588, 1.68638, 1.66505},
      {1.79599, 1.70845, 1.92716, 2.00897, 1.91856, 1.71708, 1.88995, 1.8755, 1.80035, 1.75653, 1.7559, 1.70044, 1.68223, 1.67752},
      {1.90071, 1.93324, 1.91286, 1.95086, 1.88634, 2.01884, 1.8317, 1.79525, 1.78691, 1.81176, 1.71882, 1.83073, 1.68421, 1.83502},
      {1.9406, 1.9357, 2.10884, 1.99562, 1.92118, 2.00463, 2.00567, 1.67408, 1.91014, 1.88462, 1.68, 1.79969, 1.78058, 1.71701},
      {1.75262, 1.96898, 1.9272, 1.72707, 1.87031, 1.85529, 1.79378, 1.79084, 1.92058, 1.88476, 1.83601, 1.80846, 1.62064, 1.70697},
      {1.81452, 2.112, 1.87269, 1.93471, 2.03071, 2.03796, 1.96604, 1.78397, 1.79616, 1.9119, 1.85886, 1.815, 1.83382, 1.77454},
      {1.91765, 1.71954, 1.94319, 1.95831, 1.91306, 2.06551, 1.83783, 1.8972, 2.02864, 1.80949, 1.87925, 1.78671, 1.668, 1.66126},
      {1.9297, 2.01695, 2.03827, 1.99023, 2.13117, 1.89603, 2.06998, 1.85131, 1.90398, 1.9401, 1.78263, 1.92103, 1.83979, 1.77912},
      {1.73587, 1.89583, 1.78891, 1.95814, 1.98076, 1.96337, 1.94845, 1.91691, 1.8957, 1.87147, 1.84245, 1.75678, 1.707, 1.82954},
      {2.02114, 1.9687, 2.03328, 2.00861, 1.99523, 2.09521, 1.95953, 1.81991, 1.95896, 1.98468, 1.86572, 1.88444, 1.91746, 1.90521},
      {1.73024, 2.20255, 1.67439, 2.27511, 2.06227, 1.87507, 2.13756, 2.06365, 2.0204, 1.90795, 1.98499, 1.92432, 1.99965, 1.72986},
      {1.56215, 1.98497, 1.56737, 2.48183, 1.78339, 2.07322, 2.05614, 2.19605, 1.97179, 1.86163, 1.8, 1.44447, 1.95935, 1.91278},
      {1.4819, 1.61011, 2.02632, 1.80286, 1.47899, 1.59406, 2.19378, 1.78214, 1.91822, 2.12497, 2.08859, 2.17866, 1.75998, 1.98447},
      {1.57813, 1.71088, 1.46528, 1.89427, 1.59939, 1.74955, 1.59327, 1.88198, 1.3405, 1.77137, 2.25784, 1.70967, 1.80034, 1.70439},
      {0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966},
      {0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966, 0.256966}};

   static const double alphaR[16][14] = {
      {2.41463, 2.35623, 2.34432, 2.3463, 2.37022, 2.30247, 2.27387, 2.24429, 2.19979, 2.11489, 2.04516, 1.97259, 1.90269, 1.82033},
      {2.32379, 2.20382, 2.57085, 2.6298, 2.49222, 2.27451, 2.41315, 2.1743, 2.40746, 2.34666, 2.22297, 2.19488, 2.07921, 2.01168},
      {2.30675, 2.50561, 2.41508, 2.40787, 2.76264, 3.64988, 2.34733, 2.64239, 2.37074, 2.53605, 2.15581, 2.16761, 2.00414, 1.98867},
      {2.41599, 2.28598, 2.71407, 3.10419, 2.68622, 2.61507, 2.53203, 2.27965, 2.56785, 2.43016, 2.09106, 2.26811, 2.02137, 2.083},
      {2.29075, 2.3475, 2.37196, 2.10631, 2.33236, 3.62446, 3.12258, 2.51386, 2.6385, 2.37136, 2.29904, 2.22426, 1.80177, 1.87165},
      {2.54717, 1.89453, 2.75265, 2.48901, 2.50645, 2.43424, 2.52224, 2.26171, 2.37718, 2.5228, 2.27394, 2.00614, 1.98382, 1.89367},
      {2.39843, 1.95657, 2.12871, 2.1061, 2.54437, 2.4914, 2.44042, 2.30988, 2.46239, 1.98787, 2.20715, 2.08923, 1.87958, 1.82635},
      {2.42135, 1.94849, 2.37163, 1.85904, 2.39982, 2.39721, 2.78338, 2.62555, 2.60456, 2.21436, 2.20978, 2.08877, 2.06338, 1.8301},
      {2.24403, 1.77663, 1.80067, 2.45953, 2.54245, 2.64688, 2.35636, 2.18515, 2.39744, 2.38382, 2.19851, 2.02431, 1.80367, 1.9385},
      {2.52494, 3.91036, 2.47823, 2.36698, 2.61094, 2.55736, 2.18797, 3.48102, 2.49575, 2.26932, 2.33065, 2.04055, 1.99951, 1.94202},
      {1.72235, 2.46104, 3.57679, 3.74017, 2.57926, 2.31559, 2.58968, 2.08667, 1.89485, 2.06586, 1.95675, 3.10592, 1.9657, 1.73133},
      {1.47297, 1.8934, 1.19878, 2.31738, 2.45106, 2.3279, 2.69042, 2.21169, 1.98521, 1.68679, 1.98223, 1.50562, 2.12069, 1.73507},
      {1.51367, 1.7328, 1.99804, 2.01842, 1.85997, 1.56187, 2.18322, 2.33444, 1.83771, 1.93094, 2.02453, 1.97574, 1.85458, 1.78109},
      {1.2987, 1.46178, 1.74994, 1.85632, 1.70377, 1.42185, 1.36178, 1.8501, 1.25476, 1.60448, 1.35859, 1.36298, 1.51924, 1.86623},
      {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5},
      {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}};

   static const double nL[16][14] = {
      {1.79606, 1.80007, 1.51751, 1.15843, 1.55675, 1.60987, 1.72512, 1.32482, 1.73657, 1.74101, 1.53237, 1.61929, 1.73179, 1.76587},
      {1.92807, 1.78899, 1.50603, 1.36267, 1.48862, 1.73511, 1.39185, 1.42861, 1.60639, 1.60578, 1.63959, 1.72184, 1.76286, 1.80339},
      {1.82351, 1.52936, 1.48729, 1.39348, 1.50671, 1.26892, 1.49437, 1.62143, 1.53995, 1.55014, 1.70991, 1.50131, 1.79651, 1.56445},
      {1.87872, 1.55646, 1.23943, 1.2554, 1.43989, 1.3853, 1.31464, 1.74955, 1.42716, 1.44448, 1.82317, 1.5432, 1.64493, 1.86023},
      {2.36596, 1.46554, 1.50881, 1.71557, 1.58026, 1.63684, 1.74033, 1.70555, 1.50716, 1.48744, 1.59277, 1.6969, 2.02954, 1.88922},
      {2.36922, 1.31879, 1.91213, 1.57517, 1.41925, 1.35006, 1.50747, 1.73148, 1.70456, 1.45337, 1.52972, 1.65184, 1.53516, 1.73215},
      {2.178, 1.99277, 1.57561, 1.37258, 1.59655, 1.32757, 1.69638, 1.47963, 1.38333, 1.60018, 1.55402, 1.73097, 2.02011, 2.02357},
      {2.11485, 1.53639, 1.36763, 1.40272, 1.21348, 1.66632, 1.17007, 1.63419, 1.4448, 1.38535, 1.77726, 1.32572, 1.51993, 1.69568},
      {2.81447, 1.55069, 1.76522, 1.60324, 1.49463, 1.52047, 1.39962, 1.46878, 1.60609, 1.50304, 1.55676, 1.79344, 1.90042, 1.66564},
      {2.07679, 1.75451, 1.47285, 1.60009, 1.55201, 1.37504, 1.555, 1.94241, 1.51396, 1.45406, 1.8607, 1.87955, 1.78681, 1.90265},
      {3.15878, 1.28948, 2.57319, 1.05958, 1.57625, 1.95849, 1.40414, 1.68443, 1.57472, 1.86512, 1.72561, 2.00096, 1.91188, 2.84171},
      {4.10865, 1.88003, 3.26386, 0.726548, 2.47953, 1.50357, 1.73276, 1.32945, 1.63617, 2.03821, 2.51721, 5.1019, 2.09192, 2.44653},
      {6.31061, 3.89993, 1.78775, 3.10899, 5.17758, 3.90903, 1.27795, 3.38921, 2.11321, 1.47293, 2.13959, 1.58809, 3.91754, 2.94092},
      {4.57202, 2.55234, 5.51059, 3.24901, 3.43546, 2.58996, 5.23344, 2.26803, 5.80744, 3.06346, 1.41924, 3.50837, 3.55432, 6.12389},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};

   static const double nR[16][14] = {
      {1.81397, 1.81223, 1.68456, 1.70087, 1.4376, 1.41583, 1.48218, 1.40893, 1.4467, 1.59835, 1.79925, 1.89469, 2.10231, 2.23337},
      {2.17749, 2.21512, 1.4146, 1.13063, 1.33212, 1.39447, 1.30122, 1.67219, 1.17466, 1.3015, 1.48741, 1.54638, 1.83819, 1.96385},
      {2.20741, 1.80801, 1.79239, 1.60169, 1.03894, 0.00880213, 1.397, 1.00886, 1.30609, 1.05735, 1.69891, 1.69803, 1.99793, 2.13221},
      {1.92571, 2.23982, 1.29064, 0.71453, 1.20404, 1.07102, 1.17691, 1.46607, 1.01924, 1.2821, 1.99734, 1.59851, 2.19301, 2.06306},
      {2.32737, 2.15987, 1.97782, 2.26115, 1.91674, 0.000111699, 0.424632, 1.17316, 0.928319, 1.50067, 1.65501, 1.78921, 2.92436, 2.79735},
      {1.86328, 5.82227, 1.3467, 1.56586, 1.65454, 1.57037, 1.23991, 1.8067, 1.4602, 1.13673, 1.73801, 2.38937, 2.45856, 3.01869},
      {1.88047, 3.64054, 2.37066, 2.75733, 1.27712, 1.40883, 1.44452, 1.49675, 1.34014, 2.49545, 1.9849, 2.22923, 3.06525, 3.35625},
      {2.58266, 3.87303, 3.60925, 3.17855, 1.55639, 1.62848, 0.938961, 1.06766, 1.16142, 2.171, 1.99161, 2.60985, 2.84174, 3.57068},
      {3.0321, 5.13776, 3.71342, 1.80296, 1.32248, 1.17131, 1.72412, 2.20212, 1.5456, 1.62884, 2.059, 2.56852, 3.78731, 3.13215},
      {2.25545, 0.000162555, 1.96968, 2.09914, 1.19579, 1.29054, 2.21669, 0.0209986, 1.35812, 1.66484, 1.60743, 2.6981, 2.88518, 4.16322},
      {8.63217, 2.05431, 0.000111756, 1.18242e-05, 1.24248, 1.68071, 1.09283, 2.10755, 2.7589, 2.10072, 2.91698, 0.166039, 3.09951, 5.64214},
      {11.0025, 5.17021, 29.9948, 1.90244, 1.25705, 1.85415, 0.671118, 1.87423, 2.06837, 4.71302, 2.32824, 6.91078, 2.2891, 5.52502},
      {6.17842, 4.62347, 2.47783, 2.28668, 2.93152, 6.49879, 1.94519, 1.16147, 2.82246, 2.63832, 2.13021, 2.7626, 3.85205, 5.11197},
      {13.8713, 5.1993, 7.0047, 1.9152, 2.32122, 5.55505, 4.98285, 2.05274, 7.04314, 3.59938, 6.8678, 7.41387, 4.68534, 1.91494},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};

   static const double L[16][14] = {{0.224705, 0.28372, 0.314845, 0.373415, 0.325066, 0.322566, 0.346549, 0.314114, 0.339587, 0.384313, 0.339917, 0.357749, 0.338542, 0.346226},
      {0.230575, 0.308407, 0.241138, 0.248604, 0.252084, 0.314742, 0.315059, 0.306139, 0.291012, 0.322649, 0.312485, 0.330453, 0.333721, 0.327663},
      {0.191353, 0.230623, 0.25608, 0.27072, 0.26606, 0.304564, 0.308316, 0.29009, 0.323362, 0.301312, 0.319899, 0.306165, 0.324258, 0.280476},
      {0.167616, 0.221937, 0.265645, 0.336275, 0.269117, 0.240489, 0.278749, 0.343388, 0.28218, 0.291994, 0.321479, 0.312572, 0.293505, 0.288405},
      {0.212753, 0.230115, 0.24024, 0.312425, 0.253279, 0.247811, 0.262261, 0.27156, 0.244685, 0.274085, 0.271271, 0.262416, 0.327142, 0.289958},
      {0.18383, 0.210563, 0.19385, 0.217827, 0.212065, 0.237215, 0.218728, 0.270238, 0.268405, 0.269592, 0.276056, 0.268917, 0.291118, 0.276124},
      {0.153324, 0.266153, 0.213236, 0.276491, 0.224436, 0.232436, 0.244867, 0.268869, 0.22725, 0.286639, 0.255313, 0.268615, 0.295379, 0.299422},
      {0.152747, 0.18576, 0.228638, 0.241513, 0.275296, 0.218575, 0.390105, 0.250821, 0.278478, 0.282189, 0.261866, 0.334752, 0.292486, 0.281449},
      {0.198067, 0.246237, 0.260321, 0.199551, 0.214518, 0.216529, 0.269329, 0.260297, 0.231796, 0.277108, 0.277987, 0.274966, 0.288029, 0.256543},
      {0.123771, 0.170097, 0.193866, 0.176583, 0.192532, 0.194874, 0.209649, 0.216194, 0.220725, 0.225136, 0.203286, 0.192086, 0.188411, 0.180173},
      {0.189289, 0.17883, 0.240461, 0.587556, 0.15818, 0.187867, 0.165493, 0.141819, 0.176162, 0.183058, 0.167063, 0.163097, 0.141999, 0.19978},
      {0.249745, 0.150083, 0.269312, -0.049215, 0.191587, 0.167906, 0.138895, 0.164823, 0.186696, 0.186428, 0.182408, 0.303363, 0.143418, 0.14193},
      {0.267453, 0.228491, 0.143751, 0.160986, 0.280703, 0.236613, 0.188928, 0.162644, 0.157202, 0.153287, 0.101505, 0.115489, 0.162132, 0.106581},
      {0.233479, 0.222392, 0.284983, 0.126805, 0.245458, 0.201516, 0.218059, 0.161734, 0.366947, 0.174563, 0.117197, 0.189708, 0.152864, 0.164081},
      {5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777},
      {5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777, 5.64777}};

   static const double M[16][14] = {{2.40477, 2.38118, 2.39464, 2.43925, 2.38948, 2.38055, 2.36095, 2.41107, 2.35608, 2.33347, 2.35918, 2.33546, 2.32013, 2.30032},
      {2.39053, 2.36238, 2.42629, 2.4401, 2.42174, 2.37014, 2.41316, 2.39332, 2.3965, 2.38388, 2.37463, 2.35971, 2.34354, 2.33405},
      {2.40836, 2.42462, 2.41702, 2.42251, 2.42518, 2.45178, 2.39896, 2.4053, 2.39169, 2.40482, 2.36031, 2.38464, 2.33468, 2.3647},
      {2.4214, 2.41243, 2.45448, 2.44662, 2.429, 2.43905, 2.43615, 2.36031, 2.42347, 2.41316, 2.34435, 2.38727, 2.35846, 2.35218},
      {2.37923, 2.42158, 2.41683, 2.35707, 2.40496, 2.42661, 2.41308, 2.39977, 2.42753, 2.40989, 2.39651, 2.38548, 2.2852, 2.31978},
      {2.40579, 2.39027, 2.42263, 2.42412, 2.43835, 2.43584, 2.43023, 2.38362, 2.39394, 2.42199, 2.39882, 2.36329, 2.3638, 2.33838},
      {2.41686, 2.33629, 2.39976, 2.39961, 2.42305, 2.44196, 2.40541, 2.408, 2.43607, 2.35974, 2.39695, 2.36792, 2.31176, 2.30051},
      {2.42002, 2.38749, 2.43239, 2.36928, 2.44462, 2.41329, 2.45169, 2.41542, 2.42377, 2.4074, 2.37908, 2.392, 2.37519, 2.32809},
      {2.372, 2.33917, 2.32442, 2.42626, 2.43315, 2.43428, 2.41913, 2.40119, 2.41324, 2.40835, 2.38967, 2.35386, 2.30722, 2.35638},
      {2.4379, 2.44512, 2.43741, 2.42828, 2.43764, 2.4481, 2.40796, 2.41981, 2.42806, 2.41836, 2.40404, 2.38029, 2.38038, 2.37016},
      {2.29533, 2.45464, 2.38831, 2.4777, 2.44511, 2.40465, 2.45379, 2.41142, 2.37947, 2.38726, 2.3844, 2.43618, 2.38772, 2.29726},
      {2.18201, 2.37446, 2.07092, 2.46457, 2.39538, 2.43385, 2.44784, 2.4376, 2.38662, 2.31324, 2.35709, 2.15473, 2.40133, 2.33304},
      {2.17014, 2.26787, 2.39578, 2.36255, 2.25342, 2.21932, 2.43485, 2.38844, 2.35474, 2.39749, 2.40679, 2.40944, 2.32857, 2.3536},
      {2.11988, 2.21723, 2.22717, 2.35423, 2.25827, 2.21177, 2.15023, 2.35104, 2.01828, 2.27463, 2.25818, 2.1805, 2.25533, 2.31823},
      {1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748},
      {1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748, 1.50748}};

   static const double R[16][14] = {{0.0500158, 0.0589867, 0.0672438, 0.065952, 0.0835317, 0.104406, 0.101899, 0.123716, 0.130975, 0.134953, 0.135964, 0.153421, 0.164022, 0.189748},
      {0.0534819, 0.0729384, 0.048724, 0.103657, 0.0720984, 0.116982, 0.0973578, 0.10762, 0.154029, 0.117197, 0.116024, 0.115959, 0.12145, 0.133899},
      {0.0554072, 0.0386903, 0.0507022, 0.0608959, 0.2126, -3.11455e-06, 0.0953585, 1.31289, 0.108339, 0.291699, 0.110393, 0.107107, 0.13408, 0.131088},
      {0.0465083, 0.0579475, 0.0411419, -0.00651751, 0.0595522, 0.188778, 0.106501, 0.102651, 0.763246, 0.0976073, 0.107587, 0.0899299, 0.117896, 0.106436},
      {0.0555147, 0.0504391, 0.0511809, 0.0926099, 0.0590557, -4.32767e-08, -0.00180409, 0.114369, -0.151087, 0.0759699, 0.0782089, 0.0858981, 0.166381, 0.14428},
      {0.0330522, 0.105912, 0.0319318, 0.0502057, 0.0436013, 0.0584481, 0.0851395, 0.0767285, 0.0791239, 0.136733, 0.0780522, 0.114597, 0.118758, 0.131449},
      {0.0501756, 0.103921, 0.0843025, 0.0810897, 0.071156, 0.0620902, 0.0677821, 0.0905363, 0.0771779, 0.116387, 0.0799244, 0.0978844, 0.134986, 0.147144},
      {0.0359337, 0.103655, 0.0350334, 0.139414, 0.0654606, 0.0610834, -0.114864, 0.191414, 0.0929456, 0.0721269, 0.0790964, 0.0876062, 0.0889743, 0.142215},
      {0.0536152, 0.144213, 0.150225, 0.0443459, 0.063677, 0.0777727, 0.0629245, 0.0770122, 0.0667381, 0.0634022, 0.0788974, 0.104252, 0.148102, 0.115762},
      {0.0293629, -1.98826e-08, 0.038018, 0.0490011, 0.0774057, 0.0660102, 0.0760245, -1.44019e-05, 0.0674759, 0.0840411, 0.0750983, 0.0970913, 0.103689, 0.102822},
      {0.149, 0.0383159, -5.20917e-08, -2.89895e-09, 0.0713721, 0.0730356, 0.158976, 0.103386, 0.137487, 0.109359, 0.114644, -0.000515361, 0.108794, 0.156836},
      {0.252381, 0.109054, 0.420659, 0.0620527, 0.0989525, 0.0620726, -0.0203305, 0.0840004, 0.135928, 0.181414, 0.123987, 0.249987, 0.0883729, 0.1562},
      {0.250681, 0.164093, 0.114013, 0.114832, 0.1447, 0.223462, 0.0869609, 0.202007, 0.155721, 0.129277, 0.119926, 0.112665, 0.130442, 0.142889},
      {0.35706, 0.290995, 0.144181, 0.201267, 0.24154, 0.312136, 0.363487, 0.190349, 0.422728, 0.238237, 0.342335, 0.334995, 0.263901, 0.196568},
      {1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06},
      {1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06, 1.118e-06}};


   double Mean = mean[EtaBin][PtBin];
   double Sigma = sigma[EtaBin][PtBin];
   double AlphaL = alphaL[EtaBin][PtBin];
   double AlphaR = alphaR[EtaBin][PtBin];
   double NL = nL[EtaBin][PtBin];
   double NR = nR[EtaBin][PtBin];

   double RandomNumber = -1;
   do
   {
      RandomNumber = Random.DrawDoubleSidedCBShape(Mean, Sigma, AlphaL, AlphaR, NL, NR,
         L[EtaBin][PtBin], M[EtaBin][PtBin], R[EtaBin][PtBin]);
      // RandomNumber = DrawDoubleSidedCBShape(Mean, Sigma, AlphaL, AlphaR, NL, NR);
   } while(RandomNumber < -0.25 || RandomNumber > 0.25);

   Vector = Vector * (1 + RandomNumber);
}

inline unsigned int LeptonResponseMapFindBin( const double value, const double bins[], unsigned int nbins)
{
   unsigned int nbinboundaries = nbins + 1;
   unsigned int bin = 0;
   for(unsigned int i = 0; i < nbinboundaries; i++)
   {
      if(i < nbinboundaries - 1)
      {
         if(value >= bins[i] && value < bins[i+1])
         {
            bin = i + 1;
            break;
         }
      }
      else if(i == nbinboundaries - 1)
      {
         if(value >= bins[i])
         {
            bin = nbinboundaries - 1;
            break;
         }
      }
   }

   return bin;
}

