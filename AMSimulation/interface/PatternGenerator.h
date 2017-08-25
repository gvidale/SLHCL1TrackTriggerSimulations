#ifndef AMSimulation_PatternGenerator_h_
#define AMSimulation_PatternGenerator_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternAttribute.h"
using namespace slhcl1tt;


class PatternGenerator {
  public:
    typedef PatternBankWriterT<kPatternGenerator> PatternBankWriter;
    typedef TTStubReaderT<kPatternGenerator> TTStubReader;

    // Constructor
    PatternGenerator(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {

        // Initialize
        ttmap_ = new TriggerTowerMap();
        ttmap_->read(po_.datadir);

        arbiter_ = new SuperstripArbiter();
        arbiter_->setDefinition(po_.superstrip, po_.tower, ttmap_, po.flower_charge, po.flower_pt);
    }

    // Destructor
    ~PatternGenerator() {
        if (ttmap_)     delete ttmap_;
        if (arbiter_)   delete arbiter_;
    }

    // Main driver
    int run();


  private:
    // Member functions

    // Generate pattern bank
    int makePatterns(TString src, int flower_charge);

    // Write pattern bank
    int writePatterns(TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Operators
    TriggerTowerMap   * ttmap_;
    SuperstripArbiter * arbiter_;

    // Pattern bank data
    std::map<pattern_type, unsigned>                patternBank_map_;
    std::vector<std::pair<pattern_type, unsigned> > patternBank_pairs_;

    std::map<pattern_type, PatternAttribute>        patternAttributes_map_;

    // Bookkeepers
    float coverage_;
    unsigned coverage_count_;
};

#endif
