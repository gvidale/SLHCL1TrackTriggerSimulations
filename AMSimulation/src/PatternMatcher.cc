#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternMatcher.h"


// _____________________________________________________________________________
int PatternMatcher::loadPatterns(TString bank) {
    if (verbose_)  std::cout << Info() << "Loading patterns from " << bank << std::endl;

    // _________________________________________________________________________
    // For reading pattern bank
    PatternBankReader pbreader(verbose_);
    pbreader.init(bank);

    long long npatterns = pbreader.getEntries();
    if (npatterns > po_.maxPatterns)
        npatterns = po_.maxPatterns;
    assert(npatterns > 0);

    // Setup hit buffer
    const unsigned nss = arbiter_ -> nsuperstripsPerLayer();

    if (hitBuffer_.init(po_.nLayers, nss)) {
        std::cout << Error() << "Failed to initialize HitBuffer." << std::endl;
        return 1;
    }

    // Setup associative memory
    if (associativeMemory_.init(npatterns)) {
        std::cout << Error() << "Failed to initialize AssociativeMemory." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << "Assume " << nss << " possible superstrips per layer." << std::endl;

    // _________________________________________________________________________
    // Load the patterns

    float pattInvPt = 0.;

    for (long long ipatt=0; ipatt<npatterns; ++ipatt) {
        pbreader.getPattern(ipatt);
        if (pbreader.pb_frequency < po_.minFrequency)
            break;

        assert(pbreader.pb_superstripIds->size() == po_.nLayers);

        if (verbose_>3) {
            for (unsigned i=0; i<pbreader.pb_superstripIds->size(); ++i) {
                std::cout << Debug() << "... patt: " << ipatt << "  ";
                std::copy(pbreader.pb_superstripIds->begin(), pbreader.pb_superstripIds->end(),
                          std::ostream_iterator<unsigned>(std::cout, " "));
                std::cout << " freq: " << (unsigned) pbreader.pb_frequency << std::endl;
            }
        }

        // Fill the associative memory
        pbreader.getPatternInvPt(ipatt, pattInvPt);

        associativeMemory_.insert(pbreader.pb_superstripIds->begin(), pbreader.pb_superstripIds->end(), pattInvPt, pbreader.pb_frequency);
    }

    associativeMemory_.freeze();
    assert(associativeMemory_.size() == npatterns);

    if (verbose_)  std::cout << Info() << "Successfully loaded " << npatterns << " patterns." << std::endl;

    return 0;
}

// _____________________________________________________________________________
int PatternMatcher::makeRoads(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and matching patterns." << std::endl;

    // _________________________________________________________________________
    // Get trigger tower reverse map
    const std::map<unsigned, bool>& ttrmap = ttmap_ -> getTriggerTowerReverseMap(po_.tower);


    // _________________________________________________________________________
    // For reading
    TTStubPlusTPReader reader(verbose_);
    reader.init(src);

    // For writing
    TTRoadWriter writer(verbose_);
    writer.init(reader.getChain(), out);


    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<TTRoad> roads;
    roads.reserve(300);

    std::vector<std::string> stubs_bitString;
    std::vector<unsigned>    stubs_superstripId;

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nstubs = reader.vb_modId->size();
        if (verbose_>1 && ievt%100==0)  std::cout << Debug() << Form("... Processing event: %7lld, triggering: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << std::endl;

        if (!nstubs) {  // skip if no stub
            ++nRead;
            writer.fillRoads(std::vector<TTRoad>(), std::vector<std::string>(), std::vector<unsigned>());
            continue;
        }

        if (nstubs > 500000) {
            std::cout << Error() << "Way too many stubs: " << nstubs << std::endl;
            return 1;
        }


        // _____________________________________________________________________
        // Skip stubs

        std::vector<bool> stubsNotInTower;  // true: not in this trigger tower
        std::vector<bool> stubsInOverlapping(nstubs,false);  // true: stub is in overlapping region and has TO BE removed
        for (unsigned istub=0; istub<nstubs; ++istub) {
            unsigned moduleId = reader.vb_modId   ->at(istub);

            // Skip if not in this trigger tower
            bool isNotInTower = (ttrmap.find(moduleId) == ttrmap.end());
            stubsNotInTower.push_back(isNotInTower);

            // RR // Skip if in overlapping regions
            if (po_.removeOverlap) {
                float    stub_coordx = reader.vb_coordx->at(istub);
                float    stub_coordy = reader.vb_coordy->at(istub);
                std::map<unsigned,ModuleOverlap>::iterator it_mo = momap_->moduleOverlap_map_.find(moduleId);
                if (it_mo != momap_->moduleOverlap_map_.end()) {
                    float minx = it_mo->second.x1;
                    if (stub_coordx < minx) {
                        if (verbose_>2)  std::cout << Info() << "Removing stub in module " << ievt << "\t" << moduleId << "\t x1: " <<  stub_coordx << std::endl;
                        stubsInOverlapping.at(istub)=true;
                        continue;
                    }
                    float maxx = it_mo->second.x2;
                    if (stub_coordx > maxx) {
                        if (verbose_>2)  std::cout << Info() << "Removing stub in module " << ievt << "\t"  << moduleId << "\t x2: " <<  stub_coordx << std::endl;
                        stubsInOverlapping.at(istub)=true;
                        continue;
                    }
                    float miny = it_mo->second.y1;
                    if (stub_coordy < miny) {
                        if (verbose_>2)  std::cout << Info() << "Removing stub in module " << ievt << "\t"  << moduleId << "\t y1: " <<  stub_coordy << std::endl;
                        stubsInOverlapping.at(istub)=true;
                        continue;
                    }
                    float maxy = it_mo->second.y2;
                    if (stub_coordy > maxy) {
                        if (verbose_>2)  std::cout << Info() << "Removing stub in module " << ievt << "\t"  << moduleId << "\t y2: " <<  stub_coordy << std::endl;
                        stubsInOverlapping.at(istub)=true;
                        continue;
                    }
                }
            }  // end if removeOverlap
        }

        // Null stub information for those that are not in this trigger tower
        reader.nullStubs(stubsNotInTower);

        // _____________________________________________________________________
        // Skip tracking particles

        std::vector<bool> trkPartsNotPrimary;  // true: not primary

        const unsigned nparts = reader.vp2_primary->size();
        for (unsigned ipart=0; ipart<nparts; ++ipart) {

            // Skip if not primary
            bool  primary         = reader.vp2_primary->at(ipart);
            int   simCharge       = reader.vp2_charge->at(ipart);
            trkPartsNotPrimary.push_back(!(simCharge!=0 && primary));
        }

        // Null trkPart information for those that are not primary
        reader.nullParticles(trkPartsNotPrimary);


        // _____________________________________________________________________
        // Start pattern recognition
        hitBuffer_.reset();

        stubs_bitString.clear();
        stubs_superstripId.clear();

        // Loop over reconstructed stubs
        for (unsigned istub=0; istub<nstubs; ++istub) {
            bool isNotInTower = stubsNotInTower.at(istub);
            if (isNotInTower) {
                stubs_bitString.push_back("");
                stubs_superstripId.push_back(0);
                continue;
            }

            if (po_.removeOverlap && stubsInOverlapping.at(istub)) {
                stubs_bitString.push_back("");
                stubs_superstripId.push_back(0);
                continue;
            }

            unsigned moduleId = reader.vb_modId   ->at(istub);
            float    strip    = reader.vb_coordx  ->at(istub);  // in full-strip unit
            float    segment  = reader.vb_coordy  ->at(istub);  // in full-strip unit

            float    stub_r   = reader.vb_r       ->at(istub);
            float    stub_phi = reader.vb_phi     ->at(istub);
            float    stub_z   = reader.vb_z       ->at(istub);
            float    stub_ds  = reader.vb_trigBend->at(istub);  // in full-strip unit

            // Find superstrip ID
            unsigned ssId = 0;
            if (!arbiter_ -> useGlobalCoord()) {  // local coordinates
                ssId = arbiter_ -> superstripLocal(moduleId, strip, segment);

            } else {                              // global coordinates
                ssId = arbiter_ -> superstripGlobal(moduleId, stub_r, stub_phi, stub_z, stub_ds);
            }

            unsigned lay16    = compressLayer(decodeLayer(moduleId));

            // Push into hit buffer
            hitBuffer_.insert(lay16, ssId, istub);

            if (verbose_>2) {
                std::cout << Debug() << "... ... stub: " << istub << " moduleId: " << moduleId << " strip: " << strip << " segment: " << segment << " r: " << stub_r << " phi: " << stub_phi << " z: " << stub_z << " ds: " << stub_ds << " clusRef0: " << reader.vb_clusRef0->at(istub) << " clusRef1: " << reader.vb_clusRef1->at(istub) << std::endl;
                std::cout << Debug() << "... ... stub: " << istub << " ssId: " << ssId << std::endl;
            }

            std::string bitString = "";
            stubs_bitString.push_back(bitString);
            stubs_superstripId.push_back(ssId);
        }

        hitBuffer_.freeze(po_.maxStubs);

        // _____________________________________________________________________
        // Perform associative memory lookup
        const std::vector<unsigned>& firedPatterns = associativeMemory_.lookup(hitBuffer_, po_.nLayers, po_.maxMisses);


        // _____________________________________________________________________
        // Create roads
        roads.clear();

        // Collect stubs
        for (std::vector<unsigned>::const_iterator it = firedPatterns.begin(); it != firedPatterns.end(); ++it) {
            // Create and set TTRoad
            TTRoad aroad;
            aroad.patternRef   = (*it);
            aroad.tower        = po_.tower;
            aroad.nstubs       = 0;
            aroad.patternInvPt = 0.;
            aroad.patternFreq  = 0;

            // Retrieve the superstripIds and other attributes
            pattern_type patt;
            associativeMemory_.retrieve(aroad.patternRef, patt, aroad.patternInvPt, aroad.patternFreq);

            aroad.superstripIds.clear();
            aroad.stubRefs.clear();

            aroad.superstripIds.resize(po_.nLayers);
            aroad.stubRefs.resize(po_.nLayers);

            for (unsigned layer=0; layer<po_.nLayers; ++layer) {
                const unsigned ssId     = patt.at(layer);

                if (hitBuffer_.isHit(layer, ssId)) {
                    const std::vector<unsigned>& stubRefs = hitBuffer_.getHits(layer, ssId);
                    aroad.superstripIds.at(layer) = ssId;
                    aroad.stubRefs     .at(layer) = stubRefs;
                    aroad.nstubs                 += stubRefs.size();

                } else {
                    aroad.superstripIds.at(layer) = ssId;
                }
            }

            roads.push_back(aroad);  // save aroad

            if (verbose_>2)  std::cout << Debug() << "... ... road: " << roads.size() - 1 << " " << aroad << std::endl;

            if (roads.size() >= (unsigned) po_.maxRoads)
                break;
        }

        if (! roads.empty())
            ++nKept;

//        // Sort roads by pT
//        std::stable_sort(roads.begin(), roads.end(), [](const TTRoad& lhs, const TTRoad& rhs) {
//          return (std::abs(lhs.patternInvPt) < std::abs(rhs.patternInvPt));  // higher pT first
//        }); GV - 07/25 - No roads sorting by pT, keep original pattern sorting by popularity.

        assert(reader.vb_modId->size() == stubs_bitString.size());
        assert(reader.vb_modId->size() == stubs_superstripId.size());

        writer.fillRoads(roads, stubs_bitString, stubs_superstripId);
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;

    writer.write();

    return 0;
}


// _____________________________________________________________________________
// Main driver
int PatternMatcher::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = loadPatterns(po_.bankfile);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = makeRoads(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
