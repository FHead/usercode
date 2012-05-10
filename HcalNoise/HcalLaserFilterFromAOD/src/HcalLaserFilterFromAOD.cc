// -*- C++ -*-
//
// Package:    HcalLaserFilterFromAOD
// Class:      HcalLaserFilterFromAOD
// 
/**\class HcalLaserFilterFromAOD HcalLaserFilterFromAOD.cc HcalNoise/HcalLaserFilterFromAOD/src/HcalLaserFilterFromAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yi Chen,40 3-B12,+41227675736,
//         Created:  Thu May 10 01:23:03 CEST 2012
// $Id$
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

class HcalLaserFilterFromAOD : public edm::EDFilter
{
   public:
      explicit HcalLaserFilterFromAOD(const edm::ParameterSet&);
      ~HcalLaserFilterFromAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      int RecHitCount;
      std::string NoiseSummaryLabel;
      bool TaggingMode;
};

HcalLaserFilterFromAOD::HcalLaserFilterFromAOD(const edm::ParameterSet& iConfig)
{
   RecHitCount = iConfig.getUntrackedParameter<int>("RecHitCount", 5000);
   NoiseSummaryLabel = iConfig.getUntrackedParameter<std::string>("NoiseSummaryLabel", "hcalnoise");
   TaggingMode = iConfig.getUntrackedParameter<bool>("TaggingMode", false);

   produces<bool>();
}

HcalLaserFilterFromAOD::~HcalLaserFilterFromAOD()
{
}

bool HcalLaserFilterFromAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   Handle<HcalNoiseSummary> hSummary;
   iEvent.getByLabel(NoiseSummaryLabel, hSummary);

   bool FilterDecision = true;

   if(hSummary.isValid() == false)   // not valid!  skip event
   {
      iEvent.put(std::auto_ptr<bool>(new bool(true)));
      return true;
   }

   if(hSummary->GetRecHitCount() >= RecHitCount)
   {
      std::cout << "Event rejected due to hcal laser!!!" << std::endl;
      FilterDecision = false;
   }

   iEvent.put(std::auto_ptr<bool>(new bool(FilterDecision)));

   if(TaggingMode == false && FilterDecision == false)
      return false;
   return true;
}

void HcalLaserFilterFromAOD::beginJob()
{
}

void HcalLaserFilterFromAOD::endJob()
{
}

bool HcalLaserFilterFromAOD::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

bool 
HcalLaserFilterFromAOD::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

bool 
HcalLaserFilterFromAOD::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

bool 
HcalLaserFilterFromAOD::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

void HcalLaserFilterFromAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalLaserFilterFromAOD);
