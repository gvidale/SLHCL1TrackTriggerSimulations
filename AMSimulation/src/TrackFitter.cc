#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"

namespace {
// Combination arrangement
// groups[i][j] is the j-th element in the i-th group
// combinations[i][j] is the j-th element in the i-th combination
template<typename T>
std::vector<std::vector<T> > arrangeCombinations(const std::vector<std::vector<T> >& groups) {
    std::vector<T> combination;
    std::vector<std::vector<T> > combinations;

    const int ngroups = groups.size();
    std::vector<unsigned> indices(ngroups, 0);  // init to zeroes

    int i=0, j=0;
    while (true) {
        combination.clear();
        for (i=0; i<ngroups; ++i) {
            if (groups.at(i).size())
                combination.push_back(groups.at(i).at(indices.at(i)));
            else  // empty group
                combination.push_back(999999);
        }
        combinations.push_back(combination);

        for (i=ngroups-1; i>=0; --i)
            if (groups.at(i).size())
                if (indices.at(i) != groups.at(i).size() - 1)
                    break;  // take the last index that has not reached the end
        if (i == -1)  break;

        indices[i] += 1;  // increment that index
        for (j=i+1; j<ngroups; ++j)
            indices[j] = 0;  // set indices behind that index to zeroes
    }

    return combinations;
}
}

int TrackFitter::loadConstants(TString txt) {
    if (po_.algo == "PCA4" || po_.algo == "PCA5") {
        if (fitterPCA_ -> loadConstants(txt)) {
            return 1;
        }
        if (verbose_) {
            std::cout << Info() << "The matrices are: " << std::endl;
            fitterPCA_ -> print();
        }
    }

    return 0;
}


// _____________________________________________________________________________
// Do track fitting
int TrackFitter::makeTracks(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and fitting tracks." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTRoadReader reader(verbose_);

    if (reader.init(src, prefixRoad_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTRoadReader." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // For writing
    TTTrackWriter writer(verbose_);
    if (writer.init(reader.getChain(), out, prefixTrack_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTTrackWriter." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<TTTrack2> tracks;
    tracks.reserve(300);

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nroads = reader.vr_patternRef->size();
        if (verbose_>1 && ievt%1000==0)  std::cout << Debug() << Form("... Processing event: %7lld, fitting: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

        if (!nroads) {  // skip if no road
            writer.fill(std::vector<TTTrack2>());
            ++nRead;
            continue;
        }

        tracks.clear();
        int fitstatus = 0;

        if (po_.algo == "RET") {
            // _________________________________________________________________
            // Track fitters taking the entire road

            // Loop over the roads
            for (unsigned iroad=0; iroad<nroads; ++iroad) {

                // Create and set TTRoadComb
                TTRoadComb acomb;
                acomb.roadRef    = iroad;
                acomb.patternRef = reader.vr_patternRef->at(iroad);
                acomb.tower      = reader.vr_tower->at(iroad);

                acomb.stubRefs.clear();
                for (unsigned isuperstrip=0; isuperstrip<reader.vr_stubRefs->at(iroad).size(); ++isuperstrip)
                    for (unsigned istub=0; istub<reader.vr_stubRefs->at(iroad).at(isuperstrip).size(); ++istub)
                        acomb.stubRefs.push_back(reader.vr_stubRefs->at(iroad).at(isuperstrip).at(istub));
                assert(acomb.stubRefs.size() == reader.vr_nstubs->at(iroad));

                acomb.nstubs     = 0;
                acomb.stubs_r   .clear();
                acomb.stubs_phi .clear();
                acomb.stubs_z   .clear();
                acomb.stubs_bool.clear();
                for (unsigned istub=0; istub<acomb.stubRefs.size(); ++istub) {
                    const unsigned stubRef = acomb.stubRefs.at(istub);
                    if (stubRef != 999999) {
                        ++acomb.nstubs;
                        acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));
                        acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));
                        acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));
                        acomb.stubs_bool.push_back(true);
                    } else {
                        acomb.stubs_r   .push_back(0.);
                        acomb.stubs_phi .push_back(0.);
                        acomb.stubs_z   .push_back(0.);
                        acomb.stubs_bool.push_back(false);
                    }
                }

                if (verbose_>2) {
                    std::cout << Debug() << "... ... road: " << iroad << " patternRef: " << acomb.patternRef << " tower: " << acomb.tower << " # stubs: " << acomb.nstubs << std::endl;
                    std::cout << Debug() << "... ... road: " << iroad << " stubRefs: ";
                    std::copy(acomb.stubRefs.begin(), acomb.stubRefs.end(), std::ostream_iterator<unsigned>(std::cout, " "));
                    std::cout << std::endl;
                }

                // _____________________________________________________________
                // Fit

                fitstatus = fitterRetina_->fit(acomb, ievt, tracks);
                for (unsigned itrack=0; itrack<tracks.size(); ++itrack) {
                    TTTrack2& atrack = tracks.at(itrack);
                    atrack.setTower(acomb.tower);
                    atrack.setRoadRef(acomb.roadRef);
                }

                if (verbose_>2)  std::cout << Debug() << "... ... road: " << iroad << " # tracks: " << tracks.size() << " status: " << fitstatus << std::endl;

                /// Debug
                if (verbose_>3 && tracks.size()>0) {
                    std::cout << " ----------------------------------------------------------------" << std::endl;
                    std::cout << " Generated particles:" << std::endl;

                    for (unsigned ipart=0; ipart<reader.vp_pt->size(); ++ipart) {
                        float pt      = reader.vp_pt     ->at(ipart);
                        float phi     = reader.vp_phi    ->at(ipart);
                        float eta     = reader.vp_eta    ->at(ipart);
                        float vx      = reader.vp_vx     ->at(ipart);
                        float vy      = reader.vp_vy     ->at(ipart);
                        float vz      = reader.vp_vz     ->at(ipart);
                        int   charge  = reader.vp_charge ->at(ipart);

                        // curvature and helix radius
                        float R       = pt/(0.003*3.8);

                        // helix center
                        float x0      = vx - charge*R*std::sin(phi);
                        float y0      = vy + charge*R*std::cos(phi);

                        // transverse and longitudinal impact parameters
                        float d0      = charge*(std::sqrt(x0*x0+y0*y0)-R);
                        float theta   = 2.*std::atan(std::exp(-eta));
                        float z0      = vz - 2.*R/std::tan(theta)*
                                        std::asin(0.5/R*std::sqrt((vx*vx+vy*vy-d0*d0)/(1.+d0/R)));

                        std::cout << "  " << ipart << "  -  " << pt << " " << phi << " " << eta << " " << z0 << std::endl;
                    }

                    std::cout << " Fitted tracks:" << std::endl;
                    for (unsigned itrk=0; itrk<tracks.size(); ++itrk) {
                        const TTTrack2& track = tracks.at(itrk);
                        const std::vector<unsigned>& stubRefs = track.stubRefs();

                        std::cout << "  " << itrk << "  -  " << track.pt() << " " << track.phi0() << " " << track.eta() << " " << track.z0() << " " << stubRefs.size() << std::endl;
                    }
                }
            }  // loop over the roads


        } else if (po_.algo == "ATF4" || po_.algo == "ATF5" || po_.algo == "PCA4" || po_.algo == "PCA5") {
            // _________________________________________________________________
            // Track fitters taking fit combinations

            // Loop over the roads
            for (unsigned iroad=0; iroad<nroads; ++iroad) {

                // Get combinations of stubRefs
                const unsigned nstubs = reader.vr_nstubs->at(iroad);
                std::vector<std::vector<unsigned> > combinations = arrangeCombinations(reader.vr_stubRefs->at(iroad));

                for (unsigned icomb=0; icomb<combinations.size(); ++icomb)
                    assert(combinations.at(icomb).size() == reader.vr_stubRefs->at(iroad).size());

                if (combinations.size() > (unsigned) po_.maxCombs)
                    combinations.resize(po_.maxCombs);

                if (verbose_>2) {
                    std::cout << Debug() << "... ... road: " << iroad << " # stubs: " << nstubs << " # combinations: " << combinations.size() << std::endl;
                }

                // Loop over the combinations
                for (unsigned icomb=0; icomb<combinations.size(); ++icomb) {

                    // Create and set TTRoadComb
                    TTRoadComb acomb;
                    acomb.roadRef    = iroad;
                    acomb.patternRef = reader.vr_patternRef->at(iroad);
                    acomb.tower      = reader.vr_tower->at(iroad);
                    acomb.stubRefs   = combinations.at(icomb);

                    acomb.nstubs     = 0;
                    acomb.stubs_r   .clear();
                    acomb.stubs_phi .clear();
                    acomb.stubs_z   .clear();
                    acomb.stubs_bool.clear();
                    for (unsigned istub=0; istub<acomb.stubRefs.size(); ++istub) {
                        const unsigned stubRef = acomb.stubRefs.at(istub);
                        if (stubRef != 999999) {
                            ++acomb.nstubs;
                            acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));
                            acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));
                            acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));
                            acomb.stubs_bool.push_back(true);
                        } else {
                            acomb.stubs_r   .push_back(0.);
                            acomb.stubs_phi .push_back(0.);
                            acomb.stubs_z   .push_back(0.);
                            acomb.stubs_bool.push_back(false);
                        }
                    }

                    if (verbose_>2) {
                        std::cout << Debug() << "... ... ... comb: " << icomb << " patternRef: " << acomb.patternRef << " tower: " << acomb.tower << " # stubs: " << acomb.nstubs << std::endl;
                        std::cout << Debug() << "... ... ... comb: " << icomb << " stubRefs: ";
                        std::copy(acomb.stubRefs.begin(), acomb.stubRefs.end(), std::ostream_iterator<unsigned>(std::cout, " "));
                        std::cout << std::endl;
                    }

                    // _________________________________________________________
                    // Fit
                    TTTrack2 atrack;

                    if (po_.algo=="ATF4" || po_.algo=="ATF5")
                        fitstatus = fitterATF_->fit(acomb, atrack);
                    else if (po_.algo=="PCA4" || po_.algo=="PCA5")
                        fitstatus = fitterPCA_->fit(acomb, atrack);

                    atrack.setTower    (acomb.tower);
                    atrack.setHitBits  (acomb.hitbits());
                    atrack.setPtSegment(acomb.ptsegment());
                    atrack.setRoadRef  (acomb.roadRef);
                    atrack.setStubRefs (acomb.stubRefs);
                    tracks.push_back(atrack);

                    if (verbose_>2)  std::cout << Debug() << "... ... ... track: " << icomb << " status: " << fitstatus << std::endl;
                }
            }  // loop over the roads

        }

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # tracks: " << tracks.size() << std::endl;
        if (verbose_>3) {
            for (unsigned itrack=0; itrack!=tracks.size(); ++itrack) {
                std::cout << "... ... track: " << itrack << " " << tracks.at(itrack) << std::endl;
            }
        }

        if (tracks.size() > (unsigned) po_.maxTracks)
            tracks.resize(po_.maxTracks);

        // _____________________________________________________________________
        // Remove fails and duplicates

        // FIXME: implement this

        if (! tracks.empty())
            ++nKept;

        writer.fill(tracks);
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;


    // _________________________________________________________________________
    // Write histograms

    for (std::map<TString, TH1F *>::const_iterator it=fitterPCA_->histograms.begin();
         it!=fitterPCA_->histograms.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    long long nentries = writer.writeTree();
    assert(nentries == nRead);

    return 0;
}


// _____________________________________________________________________________
// Main driver
int TrackFitter::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = loadConstants(po_.matrixfile);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = makeTracks(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
