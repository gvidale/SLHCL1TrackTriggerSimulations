#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternGenerator.h"

static const unsigned MAX_FREQUENCY = 0xffffffff;  // unsigned

namespace {
// Comparator
bool sortByFrequency(const std::pair<pattern_type, unsigned>& lhs, const std::pair<pattern_type, unsigned>& rhs) {
    return lhs.second > rhs.second;
}
}

const int iCompressedLayer[23] = {
		-1,-1,-1,-1,-1,
		0 , 1, 2, 3, 4, 5,
		6 , 7, 8, 9,10,
		-1,-1,
		11,12,13,14,15
};

//RR Parameters for the stub selection to pass exactly 6 stubs to the pattern generator
static const unsigned nLayers = 16;
static const unsigned nValidConfs = 15; //RR chosen such that ~98% of the physical stub layer configurations are covered.

int selectLayerConfiguration(bool (&bConfig)[nLayers], bool (&goodLayers)[nLayers]) {

	// One extra region added wrt M. De Mattia proposal
	const bool bValidConfigs[nLayers*nValidConfs] = {
			1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0,
			1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
			1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0,
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1,
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
	};


	int iGoodConf = -1;
	for (unsigned int iConfs=0; iConfs<nValidConfs; ++iConfs) {
		unsigned int nGoodLayers = 0;
		for (unsigned int iL=0; iL<nLayers; ++iL) {
			bool isGood = bConfig[iL]&&bValidConfigs[iL+nLayers*iConfs];
			nGoodLayers += isGood;
			goodLayers[iL] = isGood;
		}
		assert(nGoodLayers<=6);
		if (nGoodLayers==6) {
			iGoodConf = iConfs;
			break;
		}
	}
//	std::cout <<   Info() << "goodLayers[] " ;
//	for (unsigned int iL=0; iL<nLayers; ++iL) std::cout << iL << "\t" << bConfig[iL] << "\t" << goodLayers[iL] << std::endl;
	return iGoodConf;
}
// _____________________________________________________________________________
// Make the patterns
int PatternGenerator::makePatterns(TString src, int flower_charge) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and generating patterns." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTStubReader reader(verbose_);
    reader.init(src);

    // _________________________________________________________________________
    // Get trigger tower reverse map
    const std::map<unsigned, bool>& ttrmap = ttmap_ -> getTriggerTowerReverseMap(po_.tower);


    // _________________________________________________________________________
    // Loop over all events

    patternBank_map_.clear();
    pattern_type patt;
    patt.fill(0);

    PatternAttribute zero_attr;
    zero_attr.reset();

    // Bookkeepers
    float coverage = 0.;
    long int bankSize = 0, bankSizeOld = -100000, nKeptOld = -100000;
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);


        // Running estimate of coverage (10000 evt steps, to reduce load. these number  are requiered if you want a display of estimate coverage)
        if (ievt%10000 ==0){
        bankSize = patternBank_map_.size();
        coverage = 1.0 - float(bankSize - bankSizeOld) / float(nKept - nKeptOld);
        bankSizeOld = bankSize;
        nKeptOld = nKept;
        }

        if (verbose_>1 && ievt%100000==0) {
            std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7ld, # patterns: %7ld, coverage: %7.5f", ievt, nKept, bankSize, coverage) << std::endl;
        }

        const unsigned nstubs = reader.vb_modId->size();
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << std::endl;

////      GV -  Ignore event if #stubs > nLayer (AM can't handle more than that anyway')
//        if (nstubs > po_.nLayers) {
//        	if (verbose_>1) std::cout << "ERROR: #stubs > nLayers -> event ignored " << "... evt: " << ievt << " #stubs: " << nstubs << std::endl;
//        	continue;
//        }
        // Get sim info
        float simPt           = reader.vp_pt->front();
        float simEta          = reader.vp_eta->front();
        float simPhi          = reader.vp_phi->front();
        //float simVx           = reader.vp_vx->front();
        //float simVy           = reader.vp_vy->front();
        float simVz           = reader.vp_vz->front();
        int   simCharge       = reader.vp_charge->front();

        float simCotTheta     = std::sinh(simEta);
        float simChargeOverPt = float(simCharge)/simPt;

        // FOR HYBRID, FLOWER ONLY: select charge. If the option is different from 1 o -1, pass all. Default is 0 (all)
        if(flower_charge*simCharge == -1){
        	++nRead;
        	continue;
        }

        // Apply track pt requirement
        if (simPt < po_.minPt || po_.maxPt < simPt) {
            ++nRead;
            continue;
        }

        // Apply trigger tower acceptance and check for a valid layer configuration
//        std::vector<unsigned> vLayers;
//        unsigned ngoodstubs = 0;
//        for (unsigned istub=0; istub<nstubs; ++istub) {
//        	unsigned moduleId = reader.vb_modId   ->at(istub);
//        	if (ttrmap.find(moduleId) != ttrmap.end()) {
//        		++ngoodstubs;
//        	}
//        }
        bool bConfig[nLayers] = {false};
        for (unsigned istub=0; istub<nstubs; ++istub) {
            unsigned moduleId = reader.vb_modId   ->at(istub);
            if (ttrmap.find(moduleId) != ttrmap.end()) {
                unsigned iLayer = (unsigned)moduleId/10000; //HACK. FIXME there must be a standard function that returns the layerID
                assert(iCompressedLayer[iLayer]>=0);
              	bConfig[iCompressedLayer[iLayer]] = true;
//              	vLayers.push_back(iLayer);
            }
        }
        bool goodLayers[nLayers] = {false};
        int iGoodConf = selectLayerConfiguration(bConfig,goodLayers);
//        for (unsigned int iL=0; iL<nLayers; ++iL) bConfig[iL] = false;
        if (verbose_>2)  {
        	std::cout << Debug() << "... iGoodConf: " << iGoodConf << std::endl;
        	for (unsigned int iL=0; iL<nLayers; ++iL) std::cout << goodLayers[iL] << "\t";
        	std::cout << std::endl;
        }

//        if (ngoodstubs != po_.nLayers) {
        if (iGoodConf<0) {
            ++nRead;
            continue;
        }
//        assert(nstubs == po_.nLayers);
        if (verbose_>2)  {
        	std::cout << "ievt: " << ievt << "\t\t" << iGoodConf << std::endl;
        	for (unsigned int iL=0; iL<nLayers; ++iL) {
        		std::cout << goodLayers[iL] << "\t";
        	}
        	std::cout << std::endl;
        }
        // _____________________________________________________________________
        // Start generating patterns

        patt.fill(0);

        // Loop over reconstructed stubs
        unsigned igoodstub = 0;
        for (unsigned istub=0; istub<nstubs; ++istub) {
            unsigned moduleId = reader.vb_modId   ->at(istub);
            unsigned iLayer = (unsigned)moduleId/10000;//HACK. FIXME there must be a standard function that returns the layerID
            if (!goodLayers[iCompressedLayer[iLayer]]) continue;
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
//            patt.at(istub) = ssId;
            patt.at(igoodstub++) = ssId;

            if (verbose_>2) {
                std::cout << Debug() << "... ... stub: " << istub << " moduleId: " << moduleId << " strip: " << strip << " segment: " << segment << " r: " << stub_r << " phi: " << stub_phi << " z: " << stub_z << " ds: " << stub_ds << std::endl;
                std::cout << Debug() << "... ... stub: " << istub << "...igoodstub: " << igoodstub << " ssId: " << ssId << std::endl;
            }
        }
//        std::cout << std::endl;
//        std::cout << "igoodstub: " << igoodstub << std::endl;
        assert(igoodstub==6);

        // Insert pattern into the bank
        ++patternBank_map_[patt];

        // Update the attributes
        if (po_.speedup<1) {
            std::pair<std::map<pattern_type, PatternAttribute>::iterator, bool> ins = patternAttributes_map_.insert(std::make_pair(patt, zero_attr));
            PatternAttribute& attr = ins.first->second;  // pass by reference
            attr.fill(simChargeOverPt, simCotTheta, simPhi, simVz);
        }

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " patt: " << patt << std::endl;

        ++nKept;
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, kept: %7ld, # patterns: %7lu, coverage: %7.5f", nRead, nKept, patternBank_map_.size(), coverage) << std::endl;

    // Save these numbers
    coverage_       = coverage;
    coverage_count_ = nKept;


    // _________________________________________________________________________
    // Sort by frequency

    // Convert map to vector of pairs
    const unsigned origSize = patternBank_map_.size();
    //patternBank_pairs_.reserve(patternBank_map_.size());  // can cause bad_alloc
    //patternBank_pairs_.insert(patternBank_pairs_.end(), patternBank_map_.begin(), patternBank_map_.end());

    for (std::map<pattern_type, unsigned>::const_iterator it = patternBank_map_.begin();
         it != patternBank_map_.end(); ) {  // should not cause bad_alloc
        patternBank_pairs_.push_back(*it);
        it = patternBank_map_.erase(it);
    }
    assert(patternBank_pairs_.size() == origSize);

    // Clear map and release memory
    std::map<pattern_type, unsigned> mapEmpty;
    patternBank_map_.clear();
    patternBank_map_.swap(mapEmpty);

    // Sort by frequency
    std::stable_sort(patternBank_pairs_.begin(), patternBank_pairs_.end(), sortByFrequency);

    if (verbose_>2) {
        for (unsigned i=0; i<patternBank_pairs_.size(); ++i) {
            const std::pair<pattern_type, unsigned>& apair = patternBank_pairs_.at(i);
            std::cout << Debug() << "... patt: " << i << "  " << apair.first << " freq: " << apair.second << std::endl;
        }
    }

    unsigned highest_freq = patternBank_pairs_.size() ? patternBank_pairs_.front().second : 0;
    if (verbose_)  std::cout << Info() << "Generated " << patternBank_pairs_.size() << " patterns, highest freq: " << highest_freq << std::endl;
    assert(patternBank_pairs_.front().second <= MAX_FREQUENCY);

    return 0;
}


// _____________________________________________________________________________
// Output patterns into a TTree
int PatternGenerator::writePatterns(TString out) {

    // _________________________________________________________________________
    // For writing
    PatternBankWriter writer(verbose_);
    writer.init(out);

    // _________________________________________________________________________
    // Save pattern bank statistics
    *(writer.pb_coverage)   = coverage_;
    *(writer.pb_count)      = coverage_count_;
    *(writer.pb_tower)      = po_.tower;
    *(writer.pb_superstrip) = po_.superstrip;
    writer.fillPatternBankInfo();

    // _________________________________________________________________________
    // Save pattern bank
    const long long npatterns = patternBank_pairs_.size();

    // Bookkeepers
    unsigned nKept = 0;
    unsigned freq = MAX_FREQUENCY, oldFreq = MAX_FREQUENCY;
    int n90=0, n95=0, n99=0;

    for (long long ipatt=0; ipatt<npatterns; ++ipatt) {
        freq = patternBank_pairs_.at(ipatt).second;

        // Check whether patterns are indeed sorted by frequency
        assert(oldFreq >= freq);
        oldFreq = freq;
        nKept += freq;

        float coverage = float(nKept) / coverage_count_ * coverage_;
        if (!(coverage >= 0.90))
            n90 = ipatt;
        if (!(coverage >= 0.95))
            n95 = ipatt;
        if (!(coverage >= 0.99))
            n99 = ipatt;

        if (verbose_>1 && ipatt%1000==0) {
            std::cout << Debug() << Form("... Writing event: %7lld, sorted coverage: %7.5f", ipatt, coverage) << std::endl;
        }

        if (freq < (unsigned) po_.minFrequency)  // cut off
            break;

        writer.pb_superstripIds->clear();
        const pattern_type& patt = patternBank_pairs_.at(ipatt).first;
        for (unsigned ilayer=0; ilayer<po_.nLayers; ++ilayer) {
            writer.pb_superstripIds->push_back(patt.at(ilayer));
        }
        *(writer.pb_frequency) = freq;

        if (po_.speedup<1) {
            const PatternAttribute& attr = patternAttributes_map_.at(patt);
            assert(freq == attr.n);

            *(writer.pb_invPt_mean)     = attr.invPt_mean;
            *(writer.pb_invPt_sigma)    = std::sqrt(attr.invPt_variance);
            *(writer.pb_cotTheta_mean)  = attr.cotTheta_mean;
            *(writer.pb_cotTheta_sigma) = std::sqrt(attr.cotTheta_variance);
            *(writer.pb_phi_mean)       = attr.phi_mean;
            *(writer.pb_phi_sigma)      = std::sqrt(attr.phi_variance);
            *(writer.pb_z0_mean)        = attr.z0_mean;
            *(writer.pb_z0_sigma)       = std::sqrt(attr.z0_variance);
        }

        writer.fillPatternBank();
        writer.fillPatternAttributes();
    }

    writer.write();
    assert(coverage_count_ == nKept);

    if (verbose_)  {
      std::cout << Info() << "After sorting by frequency: " << std::endl;
      std::cout << Info() << " N(90% cov) = " << n90 << "\tPopularity = " << patternBank_pairs_.at(n90).second << std::endl;
      std::cout << Info() << " N(95% cov) = " << n95 << "\tPopularity = " << patternBank_pairs_.at(n95).second << std::endl;
      std::cout << Info() << " N(99% cov) = " << n99 << "\tPopularity = " << patternBank_pairs_.at(n99).second << std::endl;
    }

    return 0;
}


// _____________________________________________________________________________
// Main driver
int PatternGenerator::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = makePatterns(po_.input,po_.flower_charge);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = writePatterns(po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
} // GV test eclipse git sync
// GV eclipse git sync test
