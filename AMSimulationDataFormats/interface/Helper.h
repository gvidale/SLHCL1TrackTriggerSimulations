#ifndef AMSimulationDataFormats_Helper_h_
#define AMSimulationDataFormats_Helper_h_

#include <cmath>
#include <stdint.h>  // consistent with DataFormats/DetId/interface/DetId.h
#include <array>
//#include <tr1/array>
//#include <stdexcept>

namespace slhcl1tt {

// Typically, short is 16 bits, int is 32 bits, long is 32 bits and long long is 64 bits
// But we need to guarantee on the size of the integer
typedef uint16_t count_type; // for frequency: 0 - 2^16-1 (=65535 max)
typedef uint16_t bit_type;   // for DC bits: 0 - 2^16-1 (=65535 max for nDCBits_=4)
typedef uint32_t id_type;    // generic

typedef std::array<id_type ,8> pattern_type;     // maximum number of superstrips in a pattern set to 8 due to hardware design
typedef std::array<bit_type,8> pattern_bit_type; // one DC bit for one superstrip


// The following need to be hidden from dictionary generation
#ifndef __GCCXML__


/// Constants
// Current working assumption is:
// - PS module: 32 cols,  960 rows
// - 2S module:  2 cols, 1016 rows
// For full definition, we'll need 2^10 for submodule, 2^5 for subladder, 2^18 for moduleId,
// so 2^33, just 1 bit more than a 32-bit integer
// For almost full definition, we thus truncate 1 bit in submodule
// ModuleId is as defined by Sebastien Viret
//static const id_type iSubModuleStartBit_ = 0;
//static const id_type iSubLadderStartBit_ = 9;
//static const id_type iModuleIdStartBit_  = 14;

static const id_type iSubModuleMask_     = 0x1FF;   // 0-511 (full: 0-1023)
static const id_type iSubLadderMask_     = 0x1F;    // 0-31
static const id_type iModuleIdMask_      = 0x3FFFF; // 0-262143

// Assign a fake superstrip id
// 2^32 - 1 = 4294967295
static const id_type fakeSuperstripId_  = 0xffffffff;


/// Functions
inline id_type halfStripRound(float x) {
    static const float p = 10.;
    return floor((x*2)*p + 0.5)/p;
}

// Find the most significant bit for 32-bit word
inline id_type mostSigBit(id_type v) {
    id_type r = 0;
    while (v >>= 1)
        r++;
    return r;
}

// Retrieve layer, ladder, module from a moduleId
inline id_type decodeLayer(id_type moduleId) {
    //return (moduleId / 10000) % 100;
    return (moduleId / 10000);
}

inline id_type decodeLadder(id_type moduleId) {
    return (moduleId / 100) % 100;
}

inline id_type decodeModule(id_type moduleId) {
    return (moduleId) % 100;
}

inline id_type compressLayer(const id_type& lay) {
    if (lay <  5) return 255;
    if (lay < 16) return lay-5;  // 5-10 = barrel, 11-15 = endcap +
    if (lay < 18) return 255;
    if (lay < 23) return lay-7;  // 18-22 = endcap -
    if (lay < 25) return 255;
    if (lay < 28) return lay-9;  // 25 = calo, 26 = muon, 27 = fake
    return 255;
}

inline bool isPSModule(id_type moduleId) {
    id_type lay = decodeLayer(moduleId);
    if (5 <= lay && lay <= 7)
        return true;
    id_type lad = decodeLadder(moduleId);
    if (11 <= lay && lay <= 22 && lad <= 8)
        return true;
    return false;
}

inline bool isBarrelModule(id_type moduleId) {
    id_type lay = decodeLayer(moduleId);
    if (5 <= lay && lay <= 10)
        return true;
    return false;
}

inline id_type encodeModuleId(id_type lay, id_type lad, id_type mod) {
    return (10000*lay + 100*lad + mod);
}

#endif  // if not defined __GCCXML__

}  // namespace slhcl1tt

#endif