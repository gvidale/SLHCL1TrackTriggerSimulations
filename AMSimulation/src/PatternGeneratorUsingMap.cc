#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternGenerator.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/CSVFileReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"

static const unsigned MIN_NGOODSTUBS = 3;
static const unsigned MAX_NGOODSTUBS = 8;
static const unsigned MAX_NTOWERS = 6 * 8;
static const unsigned MAX_NMODULES = 15428;
static const unsigned MAX_FREQUENCY = 0xffff;  // 65535

namespace {
bool sortByFrequency(const std::pair<pattern_type, unsigned>& lhs, const std::pair<pattern_type, unsigned>& rhs) {
    return lhs.second > rhs.second;
}
}


// _____________________________________________________________________________
// Make the patterns
int PatternGenerator::makePatterns_map(TString src) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and generating patterns." << std::endl;

    // _________________________________________________________________________
    // For reading events
    TTStubReader reader(verbose_);
    if (reader.init(src, false)) {
        std::cout << Error() << "Failed to initialize TTStubReader." << std::endl;
        return 1;

    } else {
        TChain* tchain = reader.getChain();
        tchain->SetBranchStatus("*"                 , 0);
        tchain->SetBranchStatus("genParts_pt"       , 1);
      //tchain->SetBranchStatus("genParts_eta"      , 1);
      //tchain->SetBranchStatus("genParts_phi"      , 1);
      //tchain->SetBranchStatus("genParts_vx"       , 1);
      //tchain->SetBranchStatus("genParts_vy"       , 1);
      //tchain->SetBranchStatus("genParts_vz"       , 1);
      //tchain->SetBranchStatus("genParts_charge"   , 1);
      //tchain->SetBranchStatus("TTStubs_x"         , 1);
      //tchain->SetBranchStatus("TTStubs_y"         , 1);
      //tchain->SetBranchStatus("TTStubs_z"         , 1);
        tchain->SetBranchStatus("TTStubs_r"         , 1);
        tchain->SetBranchStatus("TTStubs_eta"       , 1);
        tchain->SetBranchStatus("TTStubs_phi"       , 1);
        tchain->SetBranchStatus("TTStubs_coordx"    , 1);
        tchain->SetBranchStatus("TTStubs_coordy"    , 1);
      //tchain->SetBranchStatus("TTStubs_roughPt"   , 1);
      //tchain->SetBranchStatus("TTStubs_trigBend"  , 1);
      //tchain->SetBranchStatus("TTStubs_clusWidth" , 1);
        tchain->SetBranchStatus("TTStubs_modId"     , 1);
      //tchain->SetBranchStatus("TTStubs_tpId"      , 1);
    }

    id_type caloSuperstrip, muonSuperstrip, fakeSuperstrip;
    if (po.mode == 2) {  // luciano superstrip
        caloSuperstrip = arbiter_ -> superstrip_luciano(25, 0, 0, po.unitPhi, po.unitEta);
        muonSuperstrip = arbiter_ -> superstrip_luciano(26, 0, 0, po.unitPhi, po.unitEta);
        fakeSuperstrip = arbiter_ -> superstrip_luciano(27, 0, 0, po.unitPhi, po.unitEta);
    } else {
        caloSuperstrip = arbiter_ -> superstrip(25, 0, 0, 0, 0);
        muonSuperstrip = arbiter_ -> superstrip(26, 0, 0, 0, 0);
        fakeSuperstrip = arbiter_ -> superstrip(27, 0, 0, 0, 0);
    }

    // Allocate memory
    allPatterns_map_.clear();  // std::map does not have reserve()


    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<id_type> stubLayers;
    pattern_type patt;
    pattern_type pattEmpty;  // for sanity check

    std::map<unsigned, unsigned> towerCountMap;

    // Bookkeepers
    float bankSize_f = 0., bankOldSize_f = 0., nKept_f = 0., nKeptOld_f = 0.;
    int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nstubs = reader.vb_modId->size();
        if (verbose_>1 && ievt%100000==0) {
            // Coverage info
            bankSize_f = allPatterns_map_.size();
            nKept_f = nKept;
            coverage_ = 1. - (bankSize_f - bankOldSize_f) / (nKept_f - nKeptOld_f);
            std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7i, # patterns: %7.0f, coverage: %7.5f", ievt, nKept, bankSize_f, coverage_) << std::endl;

            bankOldSize_f = bankSize_f;
            nKeptOld_f = nKept_f;
        }

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << std::endl;

        if (nstubs < MIN_NGOODSTUBS) {  // skip if not enough stubs
            ++nRead;
            continue;
        }

        if (nstubs > MAX_NGOODSTUBS) {
            std::cout << Error() << "Too many stubs: " << nstubs << std::endl;
            return 1;
        }

        float simPt = reader.vp_pt->front();
        if (simPt < po.minPt || po.maxPt < simPt) {
            ++nRead;
            continue;
        }

        // _____________________________________________________________________
        // Start generating patterns
        bool keep = true;

        stubLayers.clear();
        patt.fill(0);
        patt.at(6) = caloSuperstrip;  // dummy for now
        patt.at(7) = muonSuperstrip;  // dummy for now

        // Quick loop over reconstructed stubs
        id_type ssId, moduleId, lay, lad, mod, col, row;  // declare the usual suspects
        float stub_r, stub_eta, stub_phi;
        for (unsigned l=0; (l<nstubs) && keep; ++l) {
            moduleId = reader.vb_modId->at(l);

            // Skip if moduleId not in any trigger tower
            if (po.requireTriggerTower && triggerTowerReverseMap_.find(moduleId) == triggerTowerReverseMap_.end()) {
                if (verbose_>2)  std::cout << Debug() << "... ... skip moduleId: " << moduleId << " not in any trigger tower." << std::endl;
                continue;
            }

            lay = decodeLayer(moduleId);
            const unsigned& count = std::count(stubLayers.begin(), stubLayers.end(), lay);
            if (count != 0) {
                std::cout << Error() << "There should be only one stub in any layer" << std::endl;
                return 1;
            }
            stubLayers.push_back(lay);
        }

        // Decide how to lay out the pattern
        //const std::vector<unsigned>& indices = stitcher_ -> stitch(stubLayers);
        const std::vector<unsigned>& indices = stitcher_ -> stitch_layermap(stubLayers);

        if (indices.empty())
            keep = false;

        if (verbose_>2) {
            std::cout << Debug() << "... evt: " << ievt << " layers: ";
            std::copy(stubLayers.begin(), stubLayers.end(), std::ostream_iterator<unsigned>(std::cout, " "));
            std::cout << "  indices: ";
            std::copy(indices.begin(), indices.end(), std::ostream_iterator<unsigned>(std::cout, " "));
            std::cout << std::endl;
        }
        assert(!keep || indices.size() == nLayers_);

        // _____________________________________________________________________
        // Encode superstrip id
        towerCountMap.clear();  // count trigger towers for all superstrips

        // Loop over reconstructed stubs
        for (unsigned k=0, l=0; (k<nLayers_) && keep; ++k) {
            // Fake superstrip
            if (indices.at(k) == 999999) {
                // Arbitrarily put fake superstrip in layer 27
                patt.at(k) = fakeSuperstrip;
                continue;
            }

            // Real superstrip (tracker only)
            assert(k < nstubs && indices.at(k) < nstubs);
            l = indices.at(k);

            // Break moduleId into lay, lad, mod
            moduleId = reader.vb_modId->at(l);
            lay = decodeLayer(moduleId);
            lad = decodeLadder(moduleId);
            mod = decodeModule(moduleId);

            // col <-- coordy, row <-- coordx
            // use half-strip unit
            col = halfStripRound(reader.vb_coordy->at(l));
            row = halfStripRound(reader.vb_coordx->at(l));

            // global r, eta, phi
            stub_r = reader.vb_r->at(l);
            stub_eta = reader.vb_eta->at(l);
            stub_phi = reader.vb_phi->at(l);

            // Find superstrip address
            if (po.mode == 2) {  // luciano superstrip
                ssId = arbiter_ -> superstrip_luciano(lay, stub_phi, stub_eta, po.unitPhi, po.unitEta);
            } else {
                ssId = arbiter_ -> superstrip(lay, lad, mod, col, row);
            }
            patt.at(k) = ssId;

            // Find associated trigger towers
            if (po.requireTriggerTower) {
                const std::vector<unsigned>& towerIds = triggerTowerReverseMap_.at(moduleId);
                for (unsigned i=0; i<towerIds.size(); ++i) {
                    ++towerCountMap[towerIds.at(i)];
                }
                ++towerCountMap[999999];  // total count
            }

            if (verbose_>2) {
                std::cout << Debug() << "... ... stub: " << l << " moduleId: " << moduleId << " col: " << col << " row: " << row << " r: " << stub_r << " eta: " << stub_eta << " phi: " << stub_phi << std::endl;
                std::cout << Debug() << "... ... stub: " << l << " ssId: " << ssId << " towers: ";
                if (po.requireTriggerTower) {
                    const std::vector<unsigned>& towerIds = triggerTowerReverseMap_.at(moduleId);
                    std::copy(towerIds.begin(), towerIds.end(), std::ostream_iterator<unsigned>(std::cout, " "));
                }
                std::cout << std::endl;
            }
        }

        if (keep && po.requireTriggerTower) {
            // Drop patterns that are not within any trigger tower
            std::map<unsigned, unsigned>::const_iterator ittower;
            const unsigned totalTowerCount = towerCountMap.at(999999);
            towerCountMap.erase(999999);

            unsigned highestTowerCount = 0;
            for (ittower = towerCountMap.begin(); ittower != towerCountMap.end(); ++ittower) {
                if (highestTowerCount < ittower->second)
                    highestTowerCount = ittower->second;
            }

            keep = (highestTowerCount == totalTowerCount);
            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " highestTowerCount: " << highestTowerCount << " totalTowerCount: " << totalTowerCount << std::endl;
        }

        assert (!keep || patt != pattEmpty);

        // _____________________________________________________________________
        // Insert pattern
        if (keep) {
            ++allPatterns_map_[patt];
            ++nKept;

            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " patt: " << patt << std::endl;
        }
        ++nRead;
    }
    if (verbose_)  std::cout << Info() << Form("Read: %7i, kept: %7i, # patterns: %7lu, coverage: %7.5f", nRead, nKept, allPatterns_map_.size(), coverage_) << std::endl;

    // Keep this number to calculate sorted coverage
    coverage_count_ = nKept;

    // _________________________________________________________________________
    // Sort by frequency

    // Convert map to vector of pairs
    const unsigned origSize = allPatterns_map_.size();
    //allPatterns_map_pairs_.reserve(allPatterns_map_.size());  // can cause bad_alloc
    //allPatterns_map_pairs_.insert(allPatterns_map_pairs_.end(), allPatterns_map_.begin(), allPatterns_map_.end());

    for (std::map<pattern_type, unsigned>::const_iterator it = allPatterns_map_.begin();
         it != allPatterns_map_.end(); ) {  // should not cause bad_alloc
        allPatterns_map_pairs_.push_back(*it);
        it = allPatterns_map_.erase(it);
    }
    assert(allPatterns_map_pairs_.size() == origSize);

    // Clear the map and release memory
    std::map<pattern_type, unsigned> mapEmpty;
    allPatterns_map_.clear();
    allPatterns_map_.swap(mapEmpty);

    // Sort by frequency
    std::stable_sort(allPatterns_map_pairs_.begin(), allPatterns_map_pairs_.end(), sortByFrequency);
    assert(allPatterns_map_pairs_.front().second <= MAX_FREQUENCY);

    if (verbose_>2) {
        for (unsigned i=0; i<allPatterns_map_pairs_.size(); ++i) {
            const auto& pair = allPatterns_map_pairs_.at(i);
            std::cout << Debug() << "... patt: " << i << "  " << pair.first << " freq: " << pair.second << std::endl;
        }
    }

    unsigned highest_freq = allPatterns_map_pairs_.size() ? allPatterns_map_pairs_.front().second : 0;
    if (verbose_)  std::cout << Info() << "Generated " << allPatterns_map_pairs_.size() << " patterns, highest freq: " << highest_freq << std::endl;

    return 0;
}


// _____________________________________________________________________________
// Output patterns into a TTree
int PatternGenerator::writePatterns_map(TString out) {

    // _________________________________________________________________________
    // For writing
    PatternBankWriter writer(verbose_);
    if (writer.init(out)) {
        std::cout << Error() << "Failed to initialize TTRoad writer." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // Save trigger tower maps
    writer.pb_ttmap = &triggerTowerMap_;
    writer.pb_ttrmap = &triggerTowerReverseMap_;
    writer.fillTriggerTowerMaps();

    // _________________________________________________________________________
    // Save pattern bank statistics
    *(writer.pb_coverage) = coverage_;
    *(writer.pb_count) = coverage_count_;
    writer.fillPatternBankStats();

    // _________________________________________________________________________
    // Save pattern bank
    const long long nentries = allPatterns_map_pairs_.size();

    unsigned integralFreq = 0;
    count_type freq = MAX_FREQUENCY, oldFreq = MAX_FREQUENCY;

    for (long long ievt=0; ievt<nentries; ++ievt) {
        if (verbose_>1 && ievt%10000==0) {
            float coverage = float(integralFreq) / coverage_count_ * coverage_;
            std::cout << Debug() << Form("... Writing event: %7lld, sorted coverage: %7.5f", ievt, coverage) << std::endl;
        }

        freq = allPatterns_map_pairs_.at(ievt).second;

        // Make sure patterns are indeed sorted by frequency
        assert(oldFreq >= freq);
        oldFreq = freq;
        integralFreq += freq;

        if (freq < minFrequency_)  // cut off
            break;

        *(writer.pb_frequency) = freq;
        writer.pb_superstripIds->clear();
        const pattern_type& patt = allPatterns_map_pairs_.at(ievt).first;
        for (pattern_type::const_iterator it = patt.begin(); it != patt.end(); ++it) {
            writer.pb_superstripIds->push_back(*it);
        }
        writer.fillPatternBank();
    }

    long long nentries2 = writer.writeTree();
    assert(nentries == nentries2);

    return 0;
}