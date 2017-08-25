#ifndef AMSimulation_TrackParametersToTT_h_
#define AMSimulation_TrackParametersToTT_h_

#include <cmath>
#include <cassert>

// Codes originally written by Olmo Cerri (SNS Pisa)
// See https://github.com/ocerri/SLHCL1TrackTriggerSimulations/blob/dev_tt25/AMSimulation/src/TrackParametersToTT.h
// Modified for better coherence with the rest of AMSimulation

namespace slhcl1tt {

  class TrackParametersToTT {
  public:
    int get_tt(double phi, double invPt, double eta, double z0, bool apply_pt_cut=true) {
      const float range_eta[7] = {-2.4,-1.46,-0.8,0,0.8,1.46,2.4}; //GV custom size tower. (to maximize gain using flower)
      constexpr double max_eta = 2.4; //GV 08/21 FIX to include foreward disks (Range L1)
      constexpr double max_z0 = 15.; // [cm]
      constexpr double max_invPt = 1./3;  // [1/GeV]
      double etaStar = get_etaStar_from_eta(eta, z0, invPt);
      double phiStar = get_phiStar_from_phi(phi, invPt);
      if (std::abs(etaStar) > max_eta || std::abs(z0) > max_z0 || (apply_pt_cut && std::abs(invPt) > max_invPt)) {
        return -1;
      }

//      int tt_eta = (etaStar + max_eta) / (max_eta*2./6); //GV 0822: using slightly bigger foreward tower in eta, reducing size of hybrid
      int tt_eta = -1;
      for (int index = 0; index < 7; ++index){
    	  if(etaStar <= range_eta[index] ) {
    		  tt_eta = index - 1;
    		  break;
    	  }
      }
      assert(tt_eta != -1);
      int tt_phi = (phiStar + M_PI) / (M_PI*2./8);
      tt_phi = (tt_phi + 4) % 8;
      int tt = tt_eta * 8 + tt_phi;
      assert(0<=tt && tt<48);
      return tt;
    }

    double get_phiStar_from_phi(double phi, double invPt, double rStar=90.) {
      constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
      double dphi = - asin(mPtFactor * rStar * invPt);
      double phiStar = phi + dphi;
      while (phiStar < -M_PI)
        phiStar += M_PI*2.;
      while (phiStar >= M_PI)
        phiStar -= M_PI*2.;
      return phiStar;
    }

    double get_phi_from_phiStar(double phiStar, double invPt, double rStar=90.) {
      constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
      double dphi = - asin(mPtFactor * rStar * invPt);
      double phi = phiStar - dphi;
      while (phi < -M_PI)
        phi += M_PI*2.;
      while (phi >= M_PI)
        phi -= M_PI*2.;
      return phi;
    }

    double get_etaStar_from_eta(double eta, double z0, double invPt, double rStar=60.) {
      constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
      if (std::abs(rStar) < 1e-10)  return eta;
      if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
      if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
      double cot = sinh(eta);
      double cotStar = (cot * (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt)) + z0) / rStar;
      return asinh(cotStar);
    }

    double get_eta_from_etaStar(double etaStar, double z0, double invPt, double rStar=60.) {
      constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
      if (std::abs(rStar) < 1e-10)  return etaStar;
      if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
      if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
      double cotStar = sinh(etaStar);
      double cot = (rStar * cotStar - z0) / (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt));
      return asinh(cot);
    }
  };

}  // namespace slhcl1tt

#endif
