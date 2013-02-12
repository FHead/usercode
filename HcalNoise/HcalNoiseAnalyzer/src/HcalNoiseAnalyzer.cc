//---------------------------------------------------------------------------
// -*- C++ -*-
//
// Package:    HcalNoiseAnalyzer
// Class:      HcalNoiseAnalyzer
// 
/**\class HcalNoiseAnalyzer HcalNoiseAnalyzer.cc HcalNoise/HcalNoiseAnalyzer/src/HcalNoiseAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yi Chen,40 3-B12,+41227675736,
//         Created:  Wed Oct 27 11:08:10 CEST 2010
// $Id: HcalNoiseAnalyzer.cc,v 1.1 2013/02/05 12:36:37 chenyi Exp $
//
//
//---------------------------------------------------------------------------
#include <memory>
#include <string>
#include <map>
using namespace std;
//---------------------------------------------------------------------------
#include "TTree.h"
#include "TFile.h"
//---------------------------------------------------------------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/HcalNoiseRBX.h"
#include "RecoMET/METAlgorithms/interface/HcalHPDRBXMap.h"

#include "RecoMET/METAlgorithms/interface/HcalHPDRBXMap.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//---------------------------------------------------------------------------
class HcalNoiseAnalyzer;
//---------------------------------------------------------------------------
class HcalNoiseAnalyzer : public edm::EDAnalyzer
{
public:
   explicit HcalNoiseAnalyzer(const edm::ParameterSet&);
   ~HcalNoiseAnalyzer();

private:
   virtual void beginJob();
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();

private:
   bool FillHBHE;                  // Whether to store HBHE information or not
   bool FillHF;                    // Whether to store HF information or not
   double TotalChargeThreshold;    // To avoid trees from overweight, only store digis above some threshold
   string sHBHERecHitCollection;   // Name of the HBHE rechit collection
   edm::Service<TFileService> FileService;

   // Basic event coordinates
   long long RunNumber;
   long long EventNumber;
   long long LumiSection;
   long long Bunch;
   long long Orbit;
   long long Time;

   // Trigger bits
   bool TTrigger[64];
   bool L1Trigger[128];
   bool HLTrigger[256];

   // Sum \vec{ET}
   double EBET[2];
   double EEET[2];
   double HBET[2];
   double HEET[2];
   double HFET[2];
   double NominalMET[2];

   // Sum |E|
   double EBSumE;
   double EESumE;
   double HBSumE;
   double HESumE;
   double HFSumE;

   // Sum |ET|
   double EBSumET;
   double EESumET;
   double HBSumET;
   double HESumET;
   double HFSumET;

   // Summary variables for tracks and muons (and PV)
   int NumberOfGoodTracks;
   int NumberOfGoodTracks15;
   int NumberOfGoodTracks30;
   double TotalPTTracks[2];
   double SumPTTracks;
   double SumPTracks;
   int NumberOfGoodPrimaryVertices;
   int NumberOfMuonCandidates;
   int NumberOfCosmicMuonCandidates;

   // HBHE rechits and digis
   int PulseCount;
   double Charge[5184][10];
   double Pedestal[5184][10];
   double Energy[5184];
   int IEta[5184];
   int IPhi[5184];
   int Depth[5184];
   double RecHitTime[5184];
   uint32_t FlagWord[5184];
   uint32_t AuxWord[5184];

   // HF rechits and digis
   int HFPulseCount;
   int HFADC[5184][10];
   double HFCharge[5184][10];
   double HFPedestal[5184][10];
   int HFIEta[5184];
   int HFIPhi[5184];
   int HFDepth[5184];
   double HFEnergy[5184];

   // HBHE RBX energy and mega-pulse shape
   double RBXCharge[72][10];
   double RBXEnergy[72];
   double RBXCharge15[72][10];
   double RBXEnergy15[72];

   // Summary variables for baseline Hcal noise filter
   int HPDHits;
   int HPDNoOtherHits;
   int MaxZeros;
   double MinE2E10;
   double MaxE2E10;

   // Official decision from the baseline hcal noise filter
   bool OfficialDecision;

   // Summary variables for (PF) jets
   double LeadingJetEta;
   double LeadingJetPhi;
   double LeadingJetPt;
   double LeadingJetHad;
   double LeadingJetEM;
   double FollowingJetEta;
   double FollowingJetPhi;
   double FollowingJetPt;
   double FollowingJetHad;
   double FollowingJetEM;
   int JetCount20;
   int JetCount30;
   int JetCount50;
   int JetCount100;

private:
   TTree *OutputTree;

   const CaloGeometry *Geometry;

   void ClearVariables();
   void CalculateTotalEnergiesHBHE(const HBHERecHitCollection &RecHits);
   void CalculateTotalEnergiesHF(const HFRecHitCollection &RecHits);
   void CalculateTotalEnergiesEB(const EcalRecHitCollection &RecHits);
   void CalculateTotalEnergiesEE(const EcalRecHitCollection &RecHits);

   void Initialize();
};
//---------------------------------------------------------------------------
HcalNoiseAnalyzer::HcalNoiseAnalyzer(const edm::ParameterSet& iConfig)
{
   // Get stuff and initialize here
   FillHBHE = iConfig.getUntrackedParameter<bool>("FillHBHE", true);
   FillHF = iConfig.getUntrackedParameter<bool>("FillHF", false);
   TotalChargeThreshold = iConfig.getUntrackedParameter<double>("TotalChargeThreshold", 10);

   sHBHERecHitCollection = iConfig.getUntrackedParameter<string>("HBHERecHits", "hbhereco");

   Initialize();
}
//---------------------------------------------------------------------------
HcalNoiseAnalyzer::~HcalNoiseAnalyzer()
{
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   ClearVariables();

   // get stuff
   Handle<HBHERecHitCollection> hRecHits;
   iEvent.getByLabel(InputTag(sHBHERecHitCollection), hRecHits);

   Handle<HFRecHitCollection> hHFRecHits;
   iEvent.getByLabel(InputTag("hfreco"), hHFRecHits);

   Handle<HBHEDigiCollection> hHBHEDigis;
   iEvent.getByType(hHBHEDigis);

   Handle<HFDigiCollection> hHFDigis;
   iEvent.getByType(hHFDigis);

   Handle<EcalRecHitCollection> hEBRecHits;
   iEvent.getByLabel(InputTag("ecalRecHit", "EcalRecHitsEB"), hEBRecHits);

   Handle<EcalRecHitCollection> hEERecHits;
   iEvent.getByLabel(InputTag("ecalRecHit", "EcalRecHitsEE"), hEERecHits);

   Handle<CaloMETCollection> hCaloMET;
   iEvent.getByLabel(InputTag("met"), hCaloMET);

   ESHandle<HcalDbService> hConditions;
   iSetup.get<HcalDbRecord>().get(hConditions);

   ESHandle<CaloGeometry> hGeometry;
   iSetup.get<CaloGeometryRecord>().get(hGeometry);
   Geometry = hGeometry.product();

   Handle<TriggerResults> hTrigger;
   iEvent.getByLabel(InputTag("TriggerResults::HLT"), hTrigger);

   Handle<L1GlobalTriggerReadoutRecord> hL1GlobalTrigger;
   iEvent.getByLabel(InputTag("gtDigis"), hL1GlobalTrigger);

   Handle<VertexCollection> hVertices;
   iEvent.getByLabel(InputTag("offlinePrimaryVertices"), hVertices);

   Handle<TrackCollection> hTracks;
   iEvent.getByLabel(InputTag("generalTracks"), hTracks);

   Handle<MuonCollection> hStandardMuon;
   iEvent.getByLabel("muons", hStandardMuon);
   
   Handle<MuonCollection> hCosmicMuon;
   iEvent.getByLabel("muonsFromCosmics", hCosmicMuon);

   Handle<HcalNoiseSummary> hSummary;
   iEvent.getByType(hSummary);

   Handle<CaloJetCollection> hCaloJets;
   iEvent.getByLabel("ak5CaloJets", hCaloJets);

   Handle<bool> hNoiseResult;
   iEvent.getByLabel(InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"), hNoiseResult);
   OfficialDecision = *hNoiseResult;

   // basic event coordinates
   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();
   LumiSection = iEvent.luminosityBlock();
   Bunch = iEvent.bunchCrossing();
   Orbit = iEvent.orbitNumber();
   Time = iEvent.time().value();

   // event cleaning related - vertices
   NumberOfGoodPrimaryVertices = 0;
   for(int i = 0; i < (int)hVertices->size(); i++)
   {
      if((*hVertices)[i].ndof() <= 4)
         continue;
      if((*hVertices)[i].z() > 15)
         continue;
      if((*hVertices)[i].position().rho() > 2)
         continue;
      NumberOfGoodPrimaryVertices = NumberOfGoodPrimaryVertices + 1;
   }

   // event cleaning related - tracks
   for(int i = 0; i < (int)hTracks->size(); i++)
   {
      if((*hTracks)[i].numberOfValidHits() < 6)
         continue;
      if((*hTracks)[i].pt() < 0.5 || (*hTracks)[i].pt() > 500)
         continue;
      if((*hTracks)[i].eta() < -2.4 || (*hTracks)[i].eta() > 2.4)
         continue;
      if((*hTracks)[i].chi2() / (*hTracks)[i].ndof() > 20)
         continue;

      NumberOfGoodTracks = NumberOfGoodTracks + 1;
      TotalPTTracks[0] = TotalPTTracks[0] + (*hTracks)[i].px();
      TotalPTTracks[1] = TotalPTTracks[1] + (*hTracks)[i].py();
      SumPTTracks = SumPTTracks + (*hTracks)[i].pt();
      SumPTracks = SumPTracks + (*hTracks)[i].p();
      
      if((*hTracks)[i].pt() > 1.5)
         NumberOfGoodTracks15 = NumberOfGoodTracks15 + 1;
      if((*hTracks)[i].pt() > 3.0)
         NumberOfGoodTracks30 = NumberOfGoodTracks30 + 1;
   }

   // HBHE rechit maps - we want to link rechits and digis together
   map<HcalDetId, int> RecHitIndex;
   for(int i = 0; i < (int)hRecHits->size(); i++)
   {
      HcalDetId id = (*hRecHits)[i].id();
      RecHitIndex.insert(pair<HcalDetId, int>(id, i));
   }

   // HF rechit maps
   map<HcalDetId, int> HFRecHitIndex;
   for(int i = 0; i < (int)hHFRecHits->size(); i++)
   {
      HcalDetId id = (*hHFRecHits)[i].id();
      HFRecHitIndex.insert(pair<HcalDetId, int>(id, i));
   }

   // triggers
   int TTriggerSize = hL1GlobalTrigger->technicalTriggerWord().size();
   for(int i = 0; i < 64 && i < TTriggerSize; i++)
      TTrigger[i] = hL1GlobalTrigger->technicalTriggerWord()[i];
   
   int L1TriggerSize = hL1GlobalTrigger->decisionWord().size();
   for(int i = 0; i < 128 && i < L1TriggerSize; i++)
      L1Trigger[i] = hL1GlobalTrigger->decisionWord()[i];

   const TriggerNames &Names = iEvent.triggerNames(*hTrigger);
   int HLTriggerSize = Names.triggerNames().size();
   for(int i = 0; i < 256 && i < HLTriggerSize; i++)
      HLTrigger[i] = hTrigger->accept(i);
   
   // total calorimeter energies - from rechit
   CalculateTotalEnergiesHBHE(*hRecHits);
   CalculateTotalEnergiesHF(*hHFRecHits);
   CalculateTotalEnergiesEB(*hEBRecHits);
   CalculateTotalEnergiesEE(*hEERecHits);

   if(hCaloMET->size() > 0)
   {
      NominalMET[0] = (*hCaloMET)[0].px();
      NominalMET[1] = (*hCaloMET)[0].py();
   }

   // loop over digis
   for(HBHEDigiCollection::const_iterator iter = hHBHEDigis->begin(); iter != hHBHEDigis->end(); iter++)
   {
      HcalDetId id = iter->id();
     
      int RBXIndex = HcalHPDRBXMap::indexRBX(id);

      // First let's convert ADC to deposited charge
      const HcalCalibrations &Calibrations = hConditions->getHcalCalibrations(id);
      const HcalQIECoder *ChannelCoder = hConditions->getHcalCoder(id);
      const HcalQIEShape *Shape = hConditions->getHcalShape();
      HcalCoderDb Coder(*ChannelCoder, *Shape);
      CaloSamples Tool;
      Coder.adc2fC(*iter, Tool);

      // Total charge of the digi
      double TotalCharge = 0;
      for(int i = 0; i < (int)iter->size(); i++)
         TotalCharge = TotalCharge + Tool[i] - Calibrations.pedestal(iter->sample(i).capid());

      // Add this rechit/digi into RBX total charge and total energy
      for(int i = 0; i < (int)iter->size(); i++)
      {
         const HcalQIESample &QIE = iter->sample(i);
         RBXCharge[RBXIndex][i] = RBXCharge[RBXIndex][i] + Tool[i] - Calibrations.pedestal(QIE.capid());

         if((*hRecHits)[RecHitIndex[id]].energy() > 1.5)
            RBXCharge15[RBXIndex][i] = RBXCharge15[RBXIndex][i] + Tool[i] - Calibrations.pedestal(QIE.capid());
      }
      RBXEnergy[RBXIndex] = RBXEnergy[RBXIndex] + (*hRecHits)[RecHitIndex[id]].energy();
  
      if((*hRecHits)[RecHitIndex[id]].energy() > 1.5)
         RBXEnergy15[RBXIndex] = RBXEnergy15[RBXIndex] + (*hRecHits)[RecHitIndex[id]].energy();
      
      // If total charge is smaller than threshold, don't store this rechit/digi into the tree
      if(TotalCharge < TotalChargeThreshold)
         continue;

      // Safety check - there are only 5184 channels in HBHE, but just in case...
      if(PulseCount >= 5184)
      {
         PulseCount = PulseCount + 1;
         continue;
      }

      // Fill things into the tree
      for(int i = 0; i < (int)iter->size(); i++)
      {
         const HcalQIESample &QIE = iter->sample(i);

         Pedestal[PulseCount][i] = Calibrations.pedestal(QIE.capid());
         Charge[PulseCount][i] = Tool[i] - Pedestal[PulseCount][i];
      }

      Energy[PulseCount] = (*hRecHits)[RecHitIndex[id]].energy();
      RecHitTime[PulseCount] = (*hRecHits)[RecHitIndex[id]].time();

      FlagWord[PulseCount] = (*hRecHits)[RecHitIndex[id]].flags();
      AuxWord[PulseCount] = (*hRecHits)[RecHitIndex[id]].aux();

      IEta[PulseCount] = id.ieta();
      IPhi[PulseCount] = id.iphi();
      Depth[PulseCount] = id.depth();

      PulseCount = PulseCount + 1;
   }

   // Loop over HF digis
   HFPulseCount = 0;
   for(HFDigiCollection::const_iterator iter = hHFDigis->begin(); iter != hHFDigis->end(); iter++)
   {
      HcalDetId id = iter->id();

      HFIEta[HFPulseCount] = id.ieta();
      HFIPhi[HFPulseCount] = id.iphi();
      HFDepth[HFPulseCount] = id.depth();

      // ADC -> fC
      const HcalCalibrations &Calibrations = hConditions->getHcalCalibrations(id);
      const HcalQIECoder *ChannelCoder = hConditions->getHcalCoder(id);
      const HcalQIEShape *Shape = hConditions->getHcalShape();
      HcalCoderDb Coder(*ChannelCoder, *Shape);
      CaloSamples Tool;
      Coder.adc2fC(*iter, Tool);

      // Fill!
      for(int i = 0; i < iter->size(); i++)
      {
         const HcalQIESample &QIE = iter->sample(i);

         HFADC[HFPulseCount][i] = iter->sample(i).adc();

         HFPedestal[HFPulseCount][i] = Calibrations.pedestal(QIE.capid());
         HFCharge[HFPulseCount][i] = Tool[i] - HFPedestal[HFPulseCount][i];
      }

      HFEnergy[HFPulseCount] = (*hHFRecHits)[HFRecHitIndex[id]].energy();

      HFPulseCount = HFPulseCount + 1;
   }

   // muons
   NumberOfMuonCandidates = hStandardMuon->size();
   NumberOfCosmicMuonCandidates = hCosmicMuon->size();

   // hcal sumamry objects
   HPDHits = hSummary->maxHPDHits();
   HPDNoOtherHits = hSummary->maxHPDNoOtherHits();
   MaxZeros = hSummary->maxZeros();
   MinE2E10 = hSummary->minE2Over10TS();
   MaxE2E10 = hSummary->maxE2Over10TS();

   // jets
   int JetCollectionCount = hCaloJets->size();
   map<double, int, greater<double> > JetPTMap;
   for(int i = 0; i < JetCollectionCount; i++)
   {
      JetPTMap.insert(pair<double, int>((*hCaloJets)[i].pt(), i));

      if((*hCaloJets)[i].pt() > 20)
         JetCount20 = JetCount20 + 1;
      if((*hCaloJets)[i].pt() > 30)
         JetCount30 = JetCount30 + 1;
      if((*hCaloJets)[i].pt() > 50)
         JetCount50 = JetCount50 + 1;
      if((*hCaloJets)[i].pt() > 100)
         JetCount100 = JetCount100 + 1;
   }

   map<double, int, greater<double> >::iterator iter = JetPTMap.begin();
   if(JetPTMap.size() > 0)
   {
      if(iter->second < (int)hCaloJets->size())
      {
         LeadingJetEta = (*hCaloJets)[iter->second].eta();
         LeadingJetPhi = (*hCaloJets)[iter->second].phi();
         LeadingJetPt = (*hCaloJets)[iter->second].pt();
         LeadingJetHad = (*hCaloJets)[iter->second].hadEnergyInHB() + (*hCaloJets)[iter->second].hadEnergyInHE() + (*hCaloJets)[iter->second].hadEnergyInHF();
         LeadingJetEM = (*hCaloJets)[iter->second].emEnergyInEB() + (*hCaloJets)[iter->second].emEnergyInEE() + (*hCaloJets)[iter->second].emEnergyInHF();
      }
   }

   if(JetPTMap.size() > 1)
   {
      iter++;
      if(iter->second < (int)hCaloJets->size())
      {
         FollowingJetEta = (*hCaloJets)[iter->second].eta();
         FollowingJetPhi = (*hCaloJets)[iter->second].phi();
         FollowingJetPt = (*hCaloJets)[iter->second].pt();
         FollowingJetHad = (*hCaloJets)[iter->second].hadEnergyInHB() + (*hCaloJets)[iter->second].hadEnergyInHE() + (*hCaloJets)[iter->second].hadEnergyInHF();
         FollowingJetEM = (*hCaloJets)[iter->second].emEnergyInEB() + (*hCaloJets)[iter->second].emEnergyInEE() + (*hCaloJets)[iter->second].emEnergyInHF();
      }
   }

   // finally actually fill the tree
   OutputTree->Fill();
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::beginJob()
{
   // Make branches in the output trees
   OutputTree = FileService->make<TTree>("HcalNoiseTree", "Hcal noise tree version 1,1632");

   OutputTree->Branch("RunNumber", &RunNumber, "RunNumber/LL");
   OutputTree->Branch("EventNumber", &EventNumber, "EventNumber/LL");
   OutputTree->Branch("LumiSection", &LumiSection, "LumiSection/LL");
   OutputTree->Branch("Bunch", &Bunch, "Bunch/LL");
   OutputTree->Branch("Orbit", &Orbit, "Orbit/LL");
   OutputTree->Branch("Time", &Time, "Time/LL");

   OutputTree->Branch("TTrigger", &TTrigger, "TTrigger[64]/O");
   OutputTree->Branch("L1Trigger", &L1Trigger, "L1Trigger[128]/O");
   OutputTree->Branch("HLTrigger", &HLTrigger, "HLTrigger[256]/O");

   OutputTree->Branch("EBET", &EBET, "EBET[2]/D");
   OutputTree->Branch("EEET", &EEET, "EEET[2]/D");
   OutputTree->Branch("HBET", &HBET, "HBET[2]/D");
   OutputTree->Branch("HEET", &HEET, "HEET[2]/D");
   OutputTree->Branch("HFET", &HFET, "HFET[2]/D");
   OutputTree->Branch("NominalMET", &NominalMET, "NominalMET[2]/D");
   
   OutputTree->Branch("EBSumE", &EBSumE, "EBSumE/D");
   OutputTree->Branch("EESumE", &EESumE, "EESumE/D");
   OutputTree->Branch("HBSumE", &HBSumE, "HBSumE/D");
   OutputTree->Branch("HESumE", &HESumE, "HESumE/D");
   OutputTree->Branch("HFSumE", &HFSumE, "HFSumE/D");

   OutputTree->Branch("EBSumET", &EBSumET, "EBSumET/D");
   OutputTree->Branch("EESumET", &EESumET, "EESumET/D");
   OutputTree->Branch("HBSumET", &HBSumET, "HBSumET/D");
   OutputTree->Branch("HESumET", &HESumET, "HESumET/D");
   OutputTree->Branch("HFSumET", &HFSumET, "HFSumET/D");

   OutputTree->Branch("NumberOfGoodTracks", &NumberOfGoodTracks, "NumberOfGoodTracks/I");
   OutputTree->Branch("NumberOfGoodTracks15", &NumberOfGoodTracks15, "NumberOfGoodTracks15/I");
   OutputTree->Branch("NumberOfGoodTracks30", &NumberOfGoodTracks30, "NumberOfGoodTracks30/I");
   OutputTree->Branch("TotalPTTracks", &TotalPTTracks, "TotalPTTracks[2]/D");
   OutputTree->Branch("SumPTTracks", &SumPTTracks, "SumPTTracks/D");
   OutputTree->Branch("SumPTracks", &SumPTracks, "SumPTracks/D");
   OutputTree->Branch("NumberOfGoodPrimaryVertices", &NumberOfGoodPrimaryVertices, "NumberOfGoodPrimaryVertices/I");
   OutputTree->Branch("NumberOfMuonCandidates", &NumberOfMuonCandidates, "NumberOfMuonCandidates/I");
   OutputTree->Branch("NumberOfCosmicMuonCandidates", &NumberOfCosmicMuonCandidates,
      "NumberOfCosmicMuonCandidates/I");

   if(FillHBHE == true)
   {
      OutputTree->Branch("PulseCount", &PulseCount, "PulseCount/I");
      OutputTree->Branch("Charge", &Charge, "Charge[5184][10]/D");
      OutputTree->Branch("Pedestal", &Pedestal, "Pedestal[5184][10]/D");
      OutputTree->Branch("Energy", &Energy, "Energy[5184]/D");
      OutputTree->Branch("IEta", &IEta, "IEta[5184]/I");
      OutputTree->Branch("IPhi", &IPhi, "IPhi[5184]/I");
      OutputTree->Branch("Depth", &Depth, "Depth[5184]/I");
      OutputTree->Branch("RecHitTime", &RecHitTime, "RecHitTime[5184]/D");
      OutputTree->Branch("FlagWord", &FlagWord, "FlagWord[5184]/i");
      OutputTree->Branch("AuxWord", &AuxWord, "AuxWord[5184]/i");
   }

   if(FillHF == true)
   {
      OutputTree->Branch("HFPulseCount", &HFPulseCount, "HFPulseCount/I");
      OutputTree->Branch("HFADC", &HFADC, "HFADC[5184][10]/I");
      OutputTree->Branch("HFCharge", &HFCharge, "HFCharge[5184][10]/D");
      OutputTree->Branch("HFPedestal", &HFPedestal, "HFPedestal[5184][10]/D");
      OutputTree->Branch("HFIPhi", &HFIPhi, "HFIPhi[5184]/I");
      OutputTree->Branch("HFIEta", &HFIEta, "HFIEta[5184]/I");
      OutputTree->Branch("HFDepth", &HFDepth, "HFDepth[5184]/I");
      OutputTree->Branch("HFEnergy", &HFEnergy, "HFEnergy[5184]/D");
   }

   OutputTree->Branch("RBXCharge", &RBXCharge, "RBXCharge[72][10]/D");
   OutputTree->Branch("RBXEnergy", &RBXEnergy, "RBXEnergy[72]/D");
   OutputTree->Branch("RBXCharge15", &RBXCharge15, "RBXCharge15[72][10]/D");
   OutputTree->Branch("RBXEnergy15", &RBXEnergy15, "RBXEnergy15[72]/D");

   OutputTree->Branch("HPDHits", &HPDHits, "HPDHits/I");
   OutputTree->Branch("HPDNoOtherHits", &HPDNoOtherHits, "HPDNoOtherHits/I");
   OutputTree->Branch("MaxZeros", &MaxZeros, "MaxZeros/I");
   OutputTree->Branch("MinE2E10", &MinE2E10, "MinE2E10/D");
   OutputTree->Branch("MaxE2E10", &MaxE2E10, "MaxE2E10/D");

   OutputTree->Branch("LeadingJetEta", &LeadingJetEta, "LeadingJetEta/D");
   OutputTree->Branch("LeadingJetPhi", &LeadingJetPhi, "LeadingJetPhi/D");
   OutputTree->Branch("LeadingJetPt", &LeadingJetPt, "LeadingJetPt/D");
   OutputTree->Branch("LeadingJetHad", &LeadingJetHad, "LeadingJetHad/D");
   OutputTree->Branch("LeadingJetEM", &LeadingJetEM, "LeadingJetEM/D");
   OutputTree->Branch("FollowingJetEta", &FollowingJetEta, "FollowingJetEta/D");
   OutputTree->Branch("FollowingJetPhi", &FollowingJetPhi, "FollowingJetPhi/D");
   OutputTree->Branch("FollowingJetPt", &FollowingJetPt, "FollowingJetPt/D");
   OutputTree->Branch("FollowingJetHad", &FollowingJetHad, "FollowingJetHad/D");
   OutputTree->Branch("FollowingJetEM", &FollowingJetEM, "FollowingJetEM/D");
   OutputTree->Branch("JetCount20", &JetCount20, "JetCount20/I");
   OutputTree->Branch("JetCount30", &JetCount30, "JetCount30/I");
   OutputTree->Branch("JetCount50", &JetCount50, "JetCount50/I");
   OutputTree->Branch("JetCount100", &JetCount100, "JetCount100/I");

   OutputTree->Branch("OfficialDecision", &OfficialDecision, "OfficialDecision/O");
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::endJob()
{
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::ClearVariables()
{
   RunNumber = 0;
   EventNumber = 0;
   LumiSection = 0;
   Bunch = 0;
   Orbit = 0;
   Time = 0;

   for(int i = 0; i < 64; i++)
      TTrigger[i] = false;
   for(int i = 0; i < 128; i++)
      L1Trigger[i] = false;
   for(int i = 0; i < 256; i++)
      HLTrigger[i] = false;

   EBET[0] = 0;   EBET[1] = 0;
   EEET[0] = 0;   EEET[1] = 0;
   HBET[0] = 0;   HBET[1] = 0;
   HEET[0] = 0;   HEET[1] = 0;
   HFET[0] = 0;   HFET[1] = 0;
   NominalMET[0] = 0;   NominalMET[1] = 0;

   EBSumE = 0;
   EESumE = 0;
   HBSumE = 0;
   HESumE = 0;
   HFSumE = 0;

   EBSumET = 0;
   EESumET = 0;
   HBSumET = 0;
   HESumET = 0;
   HFSumET = 0;

   NumberOfGoodTracks = 0;
   NumberOfGoodTracks15 = 0;
   NumberOfGoodTracks30 = 0;
   TotalPTTracks[0] = 0;   TotalPTTracks[1] = 0;
   SumPTTracks = 0;
   SumPTracks = 0;
   NumberOfGoodPrimaryVertices = 0;
   NumberOfMuonCandidates = 0;
   NumberOfCosmicMuonCandidates = 0;

   PulseCount = 0;
   for(int i = 0; i < 5184; i++)
   {
      for(int j = 0; j < 10; j++)
      {
         Charge[i][j] = 0;
         Pedestal[i][j] = 0;
      }

      Energy[i] = 0;
      IEta[i] = 0;
      IPhi[i] = 0;
      Depth[i] = 0;
      RecHitTime[i] = 0;
      FlagWord[i] = 0;
      AuxWord[i] = 0;
   }

   HFPulseCount = 0;
   for(int i = 0; i < 5184; i++)
   {
      for(int j = 0; j < 10; j++)
      {
         HFADC[i][j] = 0;
         HFCharge[i][j] = 0;
         HFPedestal[i][j] = 0;
      }
      HFIPhi[i] = 0;
      HFIEta[i] = 0;
      HFDepth[i] = 0;
      HFEnergy[i] = 0;
   }

   for(int i = 0; i < 72; i++)
   {
      for(int j = 0; j < 10; j++)
      {
         RBXCharge[i][j] = 0;
         RBXCharge15[i][j] = 0;
      }

      RBXEnergy[i] = 0;
      RBXEnergy15[i] = 0;
   }

   HPDHits = 0;
   HPDNoOtherHits = 0;
   MaxZeros = 0;
   MinE2E10 = 0;
   MaxE2E10 = 0;

   LeadingJetEta = 0;
   LeadingJetPhi = 0;
   LeadingJetPt = 0;
   FollowingJetEta = 0;
   FollowingJetPhi = 0;
   FollowingJetPt = 0;
   JetCount20 = 0;
   JetCount30 = 0;
   JetCount50 = 0;
   JetCount100 = 0;

   OfficialDecision = false;
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::CalculateTotalEnergiesHBHE(const HBHERecHitCollection &RecHits)
{
   for(int i = 0; i < (int)RecHits.size(); i++)
   {
      bool IsHB = true;
      if(RecHits[i].id().subdet() == HcalEndcap)
         IsHB = false;

      double eta = Geometry->getPosition(RecHits[i].id()).eta();
      double phi = Geometry->getPosition(RecHits[i].id()).phi();
      double energy = RecHits[i].energy();
      double et = energy / cosh(eta);

      if(IsHB == true)
      {
         HBET[0] = HBET[0] + et * cos(phi);
         HBET[1] = HBET[1] + et * sin(phi);
         HBSumE = HBSumE + energy;
         HBSumET = HBSumET + et;
      }
      else   // is HE
      {
         HEET[0] = HEET[0] + et * cos(phi);
         HEET[1] = HEET[1] + et * sin(phi);
         HESumE = HESumE + energy;
         HESumET = HESumET + et;
      }
   }
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::CalculateTotalEnergiesHF(const HFRecHitCollection &RecHits)
{
   for(int i = 0; i < (int)RecHits.size(); i++)
   {
      double eta = Geometry->getPosition(RecHits[i].id()).eta();
      double phi = Geometry->getPosition(RecHits[i].id()).phi();
      double energy = RecHits[i].energy();
      double et = energy / cosh(eta);
         
      HFET[0] = HFET[0] + et * cos(phi);
      HFET[1] = HFET[1] + et * sin(phi);
      HFSumE = HFSumE + energy;
      HFSumET = HFSumET + et;
   }
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::CalculateTotalEnergiesEB(const EcalRecHitCollection &RecHits)
{
   for(int i = 0; i < (int)RecHits.size(); i++)
   {
      double eta = Geometry->getPosition(RecHits[i].id()).eta();
      double phi = Geometry->getPosition(RecHits[i].id()).phi();
      double energy = RecHits[i].energy();
      double et = energy / cosh(eta);
         
      EBET[0] = EBET[0] + et * cos(phi);
      EBET[1] = EBET[1] + et * sin(phi);
      EBSumE = EBSumE + energy;
      EBSumET = EBSumET + et;
   }
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::CalculateTotalEnergiesEE(const EcalRecHitCollection &RecHits)
{
   for(int i = 0; i < (int)RecHits.size(); i++)
   {
      double eta = Geometry->getPosition(RecHits[i].id()).eta();
      double phi = Geometry->getPosition(RecHits[i].id()).phi();
      double energy = RecHits[i].energy();
      double et = energy / cosh(eta);
         
      EEET[0] = EEET[0] + et * cos(phi);
      EEET[1] = EEET[1] + et * sin(phi);
      EESumE = EESumE + energy;
      EESumET = EESumET + et;
   }
}
//---------------------------------------------------------------------------
void HcalNoiseAnalyzer::Initialize()
{
   /*
   // Reads in ideal pulse shape - for fitting purposes
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
   */
}
//---------------------------------------------------------------------------
DEFINE_FWK_MODULE(HcalNoiseAnalyzer);
