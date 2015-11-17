#ifndef AMSimulation_StubCleanerPU_h_
#define AMSimulation_StubCleanerPU_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HelperMath.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Picky.h"
using namespace slhcl1tt;


class StubCleanerPU {
  public:
    // Constructor
    StubCleanerPU(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {

        // Initialize
        picky_ = new Picky();
    }

    // Destructor
    ~StubCleanerPU() {
        if (picky_)  delete picky_;
    }

    // Main driver
    int run();


  private:
    // Member functions
    // Select one unique stub per layer
    int cleanStubs(TString src, TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Picky
    Picky * picky_;
};

#endif
