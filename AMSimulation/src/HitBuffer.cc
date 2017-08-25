#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HitBuffer.h"
using namespace slhcl1tt;

#include <cassert>
#include <iostream>
#include <stdexcept>

namespace {
// Join 'layer' and 'superstrip' into one number. FIXME this need to be modified to accept SSid that already contains layer info. number with 12bit
inline unsigned simpleHash(unsigned nss, unsigned layer, unsigned ss) {  // as simple as it gets
    return nss * layer + ss;
}
}  // namespace


// _____________________________________________________________________________
int HitBuffer::init(unsigned nlayers, unsigned nss) {
    nlayers_ = nlayers;
    nss_ = nss; //initialized with arbitrer constructor, e.g. 4096. GV-> nss not in use anymore!!

    superstripHits_.clear();

    superstripBools_.clear();
//    superstripBools_.resize(nlayers_ * nss_); //with 6 e 4096 is 24596
    superstripBools_.resize(16<<12); //SSid encoding: layer - iz - iphi -> (in bit) **** - *** - *********. So in total we need 1111-111-111111111 places
    assert(nlayers_ != 0);
    assert(nss_ != 0);
    assert(superstripBools_.size() != 0);

    return 0;
}

// _____________________________________________________________________________
void HitBuffer::reset() {
    superstripHits_.clear();

    std::fill(superstripBools_.begin(), superstripBools_.end(), 0);
}

// _____________________________________________________________________________
void HitBuffer::insert(unsigned layer, superstrip_type ss, unsigned stubRef) { //GV USING DIRECTLY SSid with layer encoded.
//    const unsigned hash = simpleHash(nss_, layer, ss);
    superstripHits_[ss].push_back(stubRef);

    superstripBools_.at(ss) = true;
}

// _____________________________________________________________________________
void HitBuffer::freeze(unsigned maxStubs) {
    for (std::map<superstrip_type, std::vector<unsigned> >::iterator it=superstripHits_.begin();
         it!=superstripHits_.end(); ++it) {

        if (it->second.size() > maxStubs) {
            //it->second.resize(maxStubs);  // keep first N stubs

            it->second.erase(it->second.begin(), it->second.end() - maxStubs);  // keep last N stubs
            assert(it->second.size() == maxStubs);
        }
    }

    frozen_ = true;
}

// _____________________________________________________________________________
bool HitBuffer::isHit(unsigned layer, superstrip_type ss) const {
//    const unsigned hash = simpleHash(nss_, layer, ss);
//    return superstripBools_.at(hash);
    return superstripBools_.at(ss); //GV Just use ssId to index hash table
}

// _____________________________________________________________________________
std::vector<unsigned> HitBuffer::getHits(unsigned layer, superstrip_type ss) const {
//    const unsigned hash = simpleHash(nss_, layer, ss);
//    return superstripHits_.at(hash);
	return superstripHits_.at(ss);

}

// _____________________________________________________________________________
void HitBuffer::print() {
    std::cout << "nbins: " << superstripBools_.size() << std::endl;
}
