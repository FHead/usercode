#ifndef MixedChargeAnalysis_h_
#define MixedChargeAnalysis_h_

//
// A simple example class which uses HistogramManager and related tools
// to analyze data in a root TTree. Implementation of the class methods
// can be found in the "MixedChargeAnalysis.icc" file.
//
// I. Volobouev
// March 2013
//

#include "RootChainProcessor.h"
#include "HistogramManager.h"
#include "ChargeMixingManager.h"
#include "MixedChargeInfo.h"
#include "HBHEChannelMap.h"
#include "ChannelChargeMix.h"
#include "HcalChargeFilter.h"
#include "AbsChannelSelector.h"
#include "HBHEChannelGeometry.h"
#include "JetSummary.h"

#include "geners/AbsArchive.hh"

#include "npstat/stat/AbsNtuple.hh"

// The class template parameters are:
// 
// Options -- an (almost) arbitrary class used to pass a collection
//            of command line arguments converted into relevant C++
//            types. See comments in the MixedChargeAnalysisOptions.h
//            header file for the requirements which this class must
//            satisfy in order to play nicely with the main executable.
//
// RootMadeClass -- class generated by the root TTree "MakeClass"
//            facility. We are not going to modify that class directly
//            so that it can be trivially updated if the tree structure
//            changes in the future.
//
template <class Options, class RootMadeClass>
class MixedChargeAnalysis : public RootChainProcessor<RootMadeClass>
{
public:
    enum {
        nTimeSlices = ChannelChargeInfo::nTimeSlices
    };

    // Publish the Options type definition. This is needed for
    // command line procesing.
    typedef Options options_type;
    typedef MixedChargeAnalysis<Options, RootMadeClass> MyType;

    //
    // The constructor arguments are as follows:
    //
    // tree         -- The root tree to process (normally a TChain).
    //
    // outputfile   -- The root file into which the results will be written.
    //
    // histoRequest -- A collection of "tags" defining which output to produce.
    //                 Useful for conditional processing of various pieces of
    //                 information (and for histogram creation in particular).
    //
    // maxEvents    -- Maximum number of events to process (counted after cuts).
    //
    // verbose      -- The switch to print various information to std::cout.
    //
    // opts         -- Options filled from the command line arguments.
    //
    // Arguments "outputfile" and "histoRequest" can be simply passed to
    // an instance of HistogramManager, while "tree" and "maxEvents" will be
    // handled by the base class, RootChainProcessor.
    //
    MixedChargeAnalysis(TTree *tree, const std::string& outputfile,
                        const std::set<std::string>& histoRequest,
                        const unsigned long maxEvents, const bool verbose,
                        const Options& opts);

    // The destructor
    virtual ~MixedChargeAnalysis();

    // Simple inspectors for analysis options and "verbosity".
    // They return the arguments provided in the constructor.
    inline const Options& getOptions() const {return options_;}
    inline bool isVerbose() const {return verbose_;}

    // The following two methods are inherited from the original class
    // generated by root "MakeClass". They have the same meaning. Override
    // "Notify" if you want a notification when a new file is opened in
    // the TChain, and override "Cut" if some events should be skipped
    // (the "event" method will not be called if "Cut" returns a negative
    // number).
    virtual Bool_t Notify();
    virtual Int_t Cut(Long64_t entryNumber);

    // Channel number. Note that calling this method only makes sense
    // after "channelNumber" array has been filled.
    inline unsigned getHBHEChannelNumber(const unsigned pulseNumber) const
        {return channelNumber_[pulseNumber];}

protected:
    //
    // The methods "beginJob", "event", and "endJob" must be implemented
    // (they are pure virtual in the base class, RootChainProcessor).
    // These methods are called, well, when you expect them to be called.
    // The "entryNumber" argument of "event" (and also of "Cut") will be
    // set to the result returned by the LoadTree method of TChain. As far
    // as I understand, LoadTree returns the event number in that particular
    // tree and not in the whole TChain, so do not expect entryNumber to
    // increase monotonously (however, it should insrease monotonously
    // between two sequential "Notify" calls).
    //
    // It is expected that all of these methods will return a status: 0 on
    // success and possibly a meaningful error code on failure. If a non-0
    // status is returned, the event processing will terminate immediately
    // and the program will exit with that status.
    //
    virtual int beginJob();
    virtual int event(Long64_t entryNumber);
    virtual int endJob();

    // The "bookManagedHistograms" and "fillManagedHistograms" methods
    // just collect together various pieces of histogram creation/filling
    // code. "bookManagedHistograms" will be called from "beginJob" and
    // "fillManagedHistograms" will be called from "event".
    //
    virtual void bookManagedHistograms();
    virtual void fillManagedHistograms();

private:
    MixedChargeAnalysis();
    MixedChargeAnalysis(const MixedChargeAnalysis&);
    MixedChargeAnalysis& operator=(const MixedChargeAnalysis&);

    // Options passed to us from the main program
    const Options options_;
    const bool verbose_;

    // The histogram manager
    HistogramManager manager_;

    // Channel number mapping tool
    HBHEChannelMap channelMap_;

    // HCAL geometry tool
    HBHEChannelGeometry channelGeometry_;

    // Channel selector
    AbsChannelSelector<MyType>* channelSelector_;

    // Mask to be used for selection of good channels
    std::vector<unsigned char> channelSelectionMask_;

    // The charge mixing manager
    ChargeMixingManager<RootMadeClass>* mixManager_;

    // Object which will contain the information about mixed charge
    MixedChargeInfo mixInfo_;

    // Random number generator
    npstat::AbsRandomGenerator* rng_;

    // Filters for calculating the charge (depending on the operating
    // mode, this vector can be empty)
    std::vector<HcalChargeFilter> filters;

    // Linearized channel number (index valid up to this->PulseCount)
    unsigned channelNumber_[HBHEChannelMap::ChannelCount];

    // Number of readouts mixed with the channel (up to this->PulseCount)
    unsigned readoutsMixed_[HBHEChannelMap::ChannelCount];

    // Charge before mixing and charge added (up to this->PulseCount)
    double chargeBeforeMixing_[HBHEChannelMap::ChannelCount];
    double chargeAdded_[HBHEChannelMap::ChannelCount];

    // Similar thing for the sideband charge (up to this->PulseCount)
    double preCharge_[HBHEChannelMap::ChannelCount];
    double preChargeAdded_[HBHEChannelMap::ChannelCount];
    double postCharge_[HBHEChannelMap::ChannelCount];
    double postChargeAdded_[HBHEChannelMap::ChannelCount];

    // Charge reconstructed by the filter (up to this->PulseCount)
    double chargeReconstructed_[HBHEChannelMap::ChannelCount];

    // Summary for the locally reconstructed jets
    JetSummary jetSummary_;

    // Archive for writing out the channel data
    gs::AbsArchive* channelAr_;

    // Ntuple for the channel data
    npstat::AbsNtuple<ChannelChargeMix>* channelNtuple_;

    // Method to call in order to perform charge mixing.
    // Returns the number of channels in the combined event.
    int mixExtraCharge();

    // Method to fill the jet summary
    void fillJetSummary(JetSummary* summary);

    // Some methods used for ntuple filling
    double leadingJetHadFraction() const;
    double followingJetHadFraction() const;
};

#include "MixedChargeAnalysis.icc"

#endif // MixedChargeAnalysis_h_
