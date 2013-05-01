#ifndef NoiseTreeAnalysisOptions_h_
#define NoiseTreeAnalysisOptions_h_

#include <iostream>

#include "CmdLine.hh"

//
// Class NoiseTreeAnalysisOptions must have
//
// 1) Default constructor
//
// 2) Copy constructor (usually auto-generated)
//
// 3) Method "void parse(CmdLine& cmdline)"
//
// 4) Method "void usage(std::ostream& os) const" for printing usage
//    instructions
//
// Preferably, this class should also have "operator<<" for printing
// the option values actually used.
//
// This class works in tandem with the analysis class.
// NoiseTreeAnalysisOptions object is a "const" member in the analysis
// class, so it is safe to make NoiseTreeAnalysisOptions a struct.
//
// The "parse" method must use normal methods of "CmdLine"
// ("option", "has", and "require") to fill the members of
// this class. Note that, if you find yourself using method
// "option" to assign values to some members, you should
// initialize these members in the default constructor.
//
// Do not use here switches reserved for use by the main program.
// These switches are:
//   "-h", "--histogram"
//   "-n", "--maxEvents"
//   "-s", "--noStats"
//   "-t", "--treeName"
//   "-v", "--verbose"
//
struct NoiseTreeAnalysisOptions
{
    NoiseTreeAnalysisOptions()
        : hbGeometryFile("Geometry/hb.ctr"),
          heGeometryFile("Geometry/he.ctr"),
          maxLogContribution(10.0),
          correctionPhaseNS(6.0),
          nPhiBins(144),
          minTSlice(4),
          maxTSlice(6),
          hpdShapeNumber(105)
    {
    }

    void parse(CmdLine& cmdline)
    {
        cmdline.option(NULL, "--converters") >> convertersGSSAFile;
        cmdline.option(NULL, "--hbgeo") >> hbGeometryFile;
        cmdline.option(NULL, "--hegeo") >> heGeometryFile;
        cmdline.option(NULL, "--maxLogContribution") >> maxLogContribution;
        cmdline.option(NULL, "--correctionPhaseNS") >> correctionPhaseNS;
        cmdline.option(NULL, "--nPhiBins") >> nPhiBins;
        cmdline.option(NULL, "--minTSlice") >> minTSlice;
        cmdline.option(NULL, "--maxTSlice") >> maxTSlice;
        cmdline.option(NULL, "--hpdShapeNumber") >> hpdShapeNumber;

        if (minTSlice > 10 || maxTSlice > 10 || minTSlice >= maxTSlice)
            throw CmdLineError("Invalid specification for time slice integration");

        if (maxLogContribution < 0.0)
            throw CmdLineError("Invalid specification for maxLogContribution");
    }

    void usage(std::ostream& os) const
    {
        os << "[--converters converterFile]"
           << " [--hbgeo filename]"
           << " [--hegeo filename]"
           << " [--maxLogContribution value]"
           << " [--correctionPhaseNS value]"
           << " [--nPhiBins nBins]"
           << " [--minTSlice tSlice]"
           << " [--maxTSlice tSlice]"
           << " [--hpdShapeNumber value]"
            ;
    }

    std::string convertersGSSAFile;
    std::string hbGeometryFile;
    std::string heGeometryFile;

    double maxLogContribution;
    double correctionPhaseNS;

    unsigned nPhiBins;
    unsigned minTSlice;
    unsigned maxTSlice;

    int hpdShapeNumber;
};

std::ostream& operator<<(std::ostream& os, const NoiseTreeAnalysisOptions& o)
{
    os << "converters = \"" << o.convertersGSSAFile << '"'
       << ", hbgeo = \"" << o.hbGeometryFile << '"'
       << ", hegeo = \"" << o.heGeometryFile << '"'
       << ", maxLogContribution = " << o.maxLogContribution
       << ", correctionPhaseNS = " << o.correctionPhaseNS
       << ", nPhiBins = " << o.nPhiBins
       << ", minTSlice = " << o.minTSlice
       << ", maxTSlice = " << o.maxTSlice
       << ", hpdShapeNumber = " << o.hpdShapeNumber
        ;
    return os;
}

#endif // NoiseTreeAnalysisOptions_h_
