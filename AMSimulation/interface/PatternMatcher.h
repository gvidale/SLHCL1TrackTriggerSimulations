#ifndef AMSimulation_PatternMatcher_h_
#define AMSimulation_PatternMatcher_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/AssociativeMemory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HitBuffer.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ModuleOverlapMap.h"
using namespace slhcl1tt;


class PatternMatcher {
  public:
    typedef PatternBankReaderT<kPatternMatcher> PatternBankReader;
    typedef TTTrackReaderT<kPatternMatcher> TTStubPlusTPReader;
    typedef TTTrackWriterT<kPatternMatcher> TTRoadWriter;

    // Constructor
    PatternMatcher(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {

        // Initialize
        ttmap_ = new TriggerTowerMap();
        ttmap_->read(po_.datadir);

        arbiter_ = new SuperstripArbiter();
        arbiter_->setDefinition(po_.superstrip, po_.tower, ttmap_,po_.flower_charge, po.flower_pt);

        momap_ = 0;
        if (po_.removeOverlap) {
            momap_   = new ModuleOverlapMap();
            momap_->readModuleOverlapMap(po_.datadir);
        }
    }

    // Destructor
    ~PatternMatcher() {
        if (ttmap_)     delete ttmap_;
        if (arbiter_)   delete arbiter_;
        if (momap_)     delete momap_;
    }

    // Main driver
    int run();


  private:
    // Member functions

    // Load pattern bank
    int loadPatterns(TString bank);

    // Do pattern recognition, write roads (patterns that fired)
    int makeRoads(TString src, TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Operators
    TriggerTowerMap   * ttmap_;
    SuperstripArbiter * arbiter_;
    ModuleOverlapMap  * momap_;

    // Associative memory
    AssociativeMemory associativeMemory_;

    // Hit buffer
    HitBuffer hitBuffer_;
};

#endif
