#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
using namespace slhcl1tt;

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <boost/tokenizer.hpp>
#include <string>


static const unsigned MAX_NSTRIPS   = 1024;
static const unsigned MAX_NSEGMENTS = 32;
static const unsigned MODULE_NBITS  = 8;

namespace {
bool is_power_of_two(int x) {
    return x > 0 && !(x & (x-1));
}

unsigned most_sig_bit(unsigned v) {
    unsigned r = 0;
    while (v >>= 1)
        r++;
    return r;
}

unsigned round_to_uint(float x) {
    return std::floor(x + 0.5);
}
}


// _____________________________________________________________________________
SuperstripArbiter::SuperstripArbiter()
: sstype_(SuperstripType::UNKNOWN),
  nsuperstripsPerLayer_(999999),
  useGlobalCoord_(false),
  fixedwidth_nstrips_(0),
  fixedwidth_nz_(0),
  fixedwidth_bit_lshift1_(0),
  fixedwidth_bit_lshift2_(0),
  fixedwidth_bit_rshift1_(0),
  fixedwidth_bit_rshift2_(0),
  projective_nx_(0),
  projective_nz_(0),
  projective_max_nx_(0),
  fountain_sf_(0.),
  fountain_nz_(0),
  fountain_max_nx_(0),
  flower_arbitrer_charge (1),
  flower_reference_pt_ (5), //reference curve pT = 5 GeV/c
  flower_firstSS_phizero_ (0),
  flower_lastSS_phizero_ (0), //only for hybrid (size 10, starts from layer 11)
  flower_nr_(0.),
  flower_sf_(0.),
  flower_max_nx_(0.),
  flower_opt_(false),
  fountainopt_pt_(0) {

    // phiWidths for 6 barrel layer [0..5], 5 +z endcap disks [6..10], 5 -z endcap disks [11..15]
    // CUIDADO: dummy values are used for the endcap layers
    phiWidths_ = {
        0.00381*2, 0.00439, 0.00459, 0.00485, 0.00523, 0.00575,
        0.0048   , 0.0050 , 0.0058 , 0.0064 , 0.0070 , // RR guestimates based on the 1st barrel 2S and scattering studies on disks a-la fountain
        0.0048   , 0.0050 , 0.0058 , 0.0064 , 0.0070
    };
    assert(phiWidths_.size() == 16);

    phiWidths_floweropt_ = { //GV optimized for 2s ratio
    		0.00381*2, 0.00439, 0.00459, 0.00485, 0.00523, 0.00575,
			0.0051,	0.0062,	0.0077,	0.0088,	0.098, // GV estimates on scattering studies on disks (2s), and ratio ps 2s for TOWER 41
			0.0051,	0.0062,	0.0077,	0.0088,	0.098
    };
    //  lay 6    7	8	9	10
//  ps	0.0051	0.0056	0.006	0.0066	0.0076
//  ss	0.0052	0.0083	0.0095	0.01	0.0108
//   2s 	0	0.25	0.5	 0.67	0.7
//   opt(ps-2s-%2s)  0.0051,	0.0062,	0.0077,	0.0088,	0.098


    assert(phiWidths_.size() == 16);

    phiWidths_opt_lowpt_ = {
        0.01091, 0.00518, 0.00409, 0.00309, 0.00275, 0.00571,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999
    };
    assert(phiWidths_opt_lowpt_.size() == 16);

    phiWidths_opt_highpt_ = {
        0.01079, 0.00612, 0.00439, 0.00330, 0.00260, 0.00450,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999
    };
    assert(phiWidths_opt_highpt_.size() == 16);

    // Average radii [cm] in the 6 barrel layers
    rMeans_ = {
        22.5913, 35.4772, 50.5402, 68.3101, 88.5002, 107.71
    };
    assert(rMeans_.size() == 6);

    // Strip pitchs [um] in the 6 barrel layers
    float pitches[6] = {100., 100., 100., 90., 90., 90.};
    // Module thicknesses [cm] in the 6 barrel layers
    float thicknesses[6] = {0.26, 0.16, 0.16, 0.18, 0.18, 0.18};

    // Constants used for global phi corrections
    // r * dphi / dr ~= pitch * ds / thickness
    // dphi ~= (pitch [um] / thickness [cm] / r [cm]) * (1e2/1e6) [cm/um] * ds * dr [cm]
    // --> drCorr constant = pitch * thickness /r * 1e-4
    drCorrs_.clear();
    drCorrs_.resize(6);
    for (unsigned i=0; i<6; ++i) {
        drCorrs_[i] = pitches[i] / thicknesses[i] / rMeans_[i] * 1e-4;
    }
}

// _____________________________________________________________________________
void SuperstripArbiter::setDefinition(TString definition, unsigned tt, const TriggerTowerMap* ttmap, int flower_charge,float flower_pt) {

    // Parse the definition string -->  ff_sfb:XX_nz:YY_nr:ZZ
	typedef boost::tokenizer<boost::char_separator<char> > Tokenizer; //Class
	std::string definition_string = definition.Data();
	boost::char_separator<char> sep("_:");// default constructed
	Tokenizer tok(definition_string, sep);
	std::vector<TString> token;
//	std::cout << "tokens:" << std::endl;
	for(Tokenizer::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter){
		token.push_back(*tok_iter);
//		std::cout<< *tok_iter << std::endl;//maybe need to be stringed. (what the iterator cointains)
	}
	

	int token_size = token.size();
//	std::cout << "#token pinolinolup " << token_size << std::endl;
//	OLD
//    if (definition.Contains("_")) {
//        unsigned pos = definition.Index("_");
//        //unsigned len = definition.Length();
//        TString token1 = definition(0,2); //type SS definition
//        TString token2 = definition(2,pos-2); //number scale factor (barrel)
//        TString token3 = definition(pos+1,2);  //"nz"
//	// TString token4 = definition(pos+3,len-(pos+3));
//        TString token4 = definition(pos+3,1); // number nz

	TString token5, token6;

	int posL5 = definition.Index("L5x");
	if (posL5 != -1) token5 = definition(posL5+3,1);
	int posL10 = definition.Index("L10x");
	if (posL10 != -1) token6 = definition(posL10+4,1);
	
	//Need 9 tokens for fountain and flower, barrel + disk independent
        if (token[0] == "fw") { //factor1 ss
            sstype_ = SuperstripType::FIXEDWIDTH;
        } else if (token[0] == "proj") {//factor1 nx
            sstype_ = SuperstripType::PROJECTIVE;
        } else if (token[0] == "f") {//factor1 sf
            sstype_ = SuperstripType::FOUNTAIN;
        } else if (token[0] == "fopgv") {//factor1 sf
        	sstype_ = SuperstripType::FOUNTAIN;
        	flower_opt_ = true;
        } else if (token[0] == "fop") {//sf
            sstype_ = SuperstripType::FOUNTAINOPT;
        } else if (token[0] == "ff") {//sf
        	sstype_ = SuperstripType::FOUNTAIN_FLOWER;
        } else if (token[0] == "ffopgv") {//sf
          	sstype_ = SuperstripType::FOUNTAIN_FLOWER;
          	flower_opt_ = true;
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }

        //Check validity of arguments
        if (token_size!= 5 &&( token[0]=="fw"||token[0]=="proj")){
        	throw std::invalid_argument("Incorrect superstrip definition. Wrong arguments for proj or fw.");
        }
        if (token_size!= 9 &&( token[0]!="fw"&&token[0]!="proj")){
        	throw std::invalid_argument("Incorrect superstrip definition. Wrong arguments for fountain or flower.");
        }

        if(token[1]=="sf" && (token[0]=="f"||token[0]=="ff"||token[0]=="fop" ||token[0]=="fopgv"||token[0]=="ffopgv") ){
        	// Do nothing
        } else {
        	throw std::invalid_argument("Incorrect superstrip definition iwith scale factor ff (did you mean ss or nx?).");
        }

        if (token[0]=="fopgv" || token[0]=="ffopgv") std::cout << "using optimized 2s ratio disks phi widths" << std::endl;

        token[2].ReplaceAll("p", ".");
        if (token[2].IsFloat()) {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }
        float token2f = token[2].Atof();

        if (token[3] == "nz") {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }

        if (token[4].IsFloat()) {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }
        float token4f = token[4].Atof();

        float tokenFFnr6f=0;
        float tokenFFsf7f=0;

        //Specific for fountain and flower, 9 tokens.
        if(token_size == 9){
        	if (token[5] == "nr") {
        		// Do nothing
        	} else {
        		throw std::invalid_argument("Incorrect superstrip definition.");
        	}

        	if (token[6].IsFloat()) { //scale factor fountain
        		// Do nothing
        	} else {
        		throw std::invalid_argument("Incorrect superstrip definition.");
        	}
        	tokenFFnr6f = token[6].Atof();

        	if (token[7] == "sf") {
        		// Do nothing
        	} else {
        		throw std::invalid_argument("Incorrect superstrip definition.");
        	}

        	if (token[8].IsFloat()) { //scale factor fountain
        		// Do nothing
        	} else {
        		throw std::invalid_argument("Incorrect superstrip definition.");
        	}
        	tokenFFsf7f = token[8].Atof();

        }


	if (token5.IsFloat()) {
            // Do nothing
        } else {
	  token5="1";
        }
        float token5f = token5.Atof();

	if (token6.IsFloat()) {
            // Do nothing
        } else {
        	token6="1";
        }
	float token6f = token6.Atof();

	int n_flowerstrips=0; //debug: counts SS for flower

	switch (sstype_) {
	case SuperstripType::FIXEDWIDTH:
		fixedwidth_nstrips_     = token2f;
		fixedwidth_nz_          = token4f;
		useGlobalCoord_         = false;

		if (!is_power_of_two(fixedwidth_nstrips_) || !is_power_of_two(fixedwidth_nz_) ||
				fixedwidth_nstrips_ == 0 || fixedwidth_nz_ == 0 ||
				fixedwidth_nstrips_ > MAX_NSTRIPS || fixedwidth_nz_ > MAX_NSEGMENTS)
			throw std::invalid_argument("Incorrect fixedwidth superstrip definition.");
		break;

	case SuperstripType::PROJECTIVE:
		projective_nx_  = token2f;
		projective_nz_  = token4f;
		useGlobalCoord_ = true;

		if (projective_nx_ == 0 || projective_nz_ == 0)
			throw std::invalid_argument("Incorrect projective superstrip definition.");
		break;

	case SuperstripType::FOUNTAIN:
		fountain_sf_    = token2f;
		fountain_nz_    = token4f;
		flower_nr_	    = tokenFFnr6f;
		flower_sf_      = tokenFFsf7f;
		useGlobalCoord_ = true;
		fountain_xfactor_.clear();
		fountain_xfactor_.resize(16,1);

		fountain_xfactor_ .at(0)=token5f; //L5
		fountain_xfactor_ .at(5)=token6f; //L10

		if (fountain_sf_ <= 0. || fountain_nz_ == 0 || flower_sf_ <=0. || flower_nr_ == 0)
			throw std::invalid_argument("Incorrect fountain superstrip definition.pippolo");
		break;

	case SuperstripType::FOUNTAIN_FLOWER:
		fountain_sf_    = token2f;
		fountain_nz_    = token4f;
		flower_nr_	    = tokenFFnr6f;
		flower_sf_      = tokenFFsf7f;
		std::cout << " f-flower sf_nz_nr_sf: " << fountain_sf_ << fountain_nz_ << flower_nr_ << flower_sf_ <<std::endl;
		useGlobalCoord_ = true;
		fountain_xfactor_.clear();
		fountain_xfactor_.resize(16,1);

		fountain_xfactor_ .at(0)=token5f; //L5
		fountain_xfactor_ .at(5)=token6f; //L10

		if (fountain_sf_ <= 0. || fountain_nz_ == 0 || flower_sf_ <=0. || flower_nr_ == 0)
			throw std::invalid_argument("Incorrect flower superstrip definition.");
		flower_arbitrer_charge = (flower_charge>=0)? 1 : -1; //if the option charge is positive or zero, use positive configuration (flower anticlockwise), otherwise negative
		flower_reference_pt_ = flower_pt;
		break;

	case SuperstripType::FOUNTAINOPT:
		fountainopt_pt_ = token2f;
		fountain_nz_    = token4f;
		useGlobalCoord_ = true;

		if (fountainopt_pt_ <= 0. || fountain_nz_ == 0)
			throw std::invalid_argument("Incorrect fountain superstrip definition.");
		break;

	case SuperstripType::UNKNOWN:
	default:
		throw std::invalid_argument("Incorrect superstrip definition.");
		break;
	}

//        } else {
//        	throw std::invalid_argument("Incorrect superstrip definition.");
//        }

	// Learn about the trigger tower geometry
	if (!useGlobalCoord()) {
		towerModules_.clear();
		towerModules_.resize(16);

		const std::vector<unsigned>& moduleIds = ttmap->getTriggerTowerModules(tt);
		for (unsigned i=0; i<moduleIds.size(); ++i) {
			unsigned moduleId = moduleIds.at(i);
			unsigned lay16    = compressLayer(decodeLayer(moduleId));
			assert(lay16 < 16);

			towerModules_.at(lay16).push_back(moduleId);
		}

		// Sanity check
		for (std::vector<std::vector<unsigned> >::const_iterator it = towerModules_.begin();
				it != towerModules_.end(); ++it) {
			assert(towerModules_.size() < (1u << MODULE_NBITS));  // make sure number of bits is enough
		}

	} else {
		phiMins_.clear();
		phiMaxs_.clear();
		zMins_.clear();
		zMaxs_.clear();

		phiMins_.resize(16);
		phiMaxs_.resize(16);
		zMins_.resize(16);
		zMaxs_.resize(16);

		const std::vector<LayerBounds>& boundaries = ttmap->getTriggerTowerBoundaries(tt);
		for (unsigned i=0; i<boundaries.size(); ++i) {
			const LayerBounds& lb = boundaries.at(i);
			unsigned lay16        = compressLayer(lb.layer);
			assert(lay16 < 16);

			phiMins_.at(lay16) = lb.phiMin;
			phiMaxs_.at(lay16) = lb.phiMax;
			zMins_  .at(lay16) = lb.zMin;
			zMaxs_  .at(lay16) = lb.zMax;
			std::cout << " layer " << lay16 << " phiMin " << phiMins_.at(lay16)
					<< " phiMax " <<phiMaxs_.at(lay16)<<" zMin " << zMins_  .at(lay16)<< " zMax " << zMaxs_  .at(lay16) << std::endl;
		}
//		std::cout << "setting boundaries ok" <<std::endl;
	}

	// Set number of possible superstrips per layer
	switch(sstype_) {
	case SuperstripType::FIXEDWIDTH:
		fixedwidth_bit_rshift1_ = most_sig_bit(fixedwidth_nstrips_);
		fixedwidth_bit_rshift2_ = most_sig_bit(MAX_NSEGMENTS / fixedwidth_nz_);

		fixedwidth_bit_lshift1_ = most_sig_bit(MAX_NSTRIPS / fixedwidth_nstrips_);
		fixedwidth_bit_lshift2_ = most_sig_bit(fixedwidth_nz_);

		nsuperstripsPerLayer_ = (1u << (MODULE_NBITS + fixedwidth_bit_lshift2_ + fixedwidth_bit_lshift1_));
		break;

	case SuperstripType::PROJECTIVE:
		projective_phiBins_.clear();
		projective_zBins_.clear();

		projective_phiBins_.resize(16);
		projective_zBins_.resize(16);

		for (unsigned i=0; i<phiMins_.size(); ++i) {
			projective_phiBins_.at(i) = M_PI / 4. / projective_nx_;  // constant
			//projective_phiBins_.at(i) = (phiMaxs_.at(i) - phiMins_.at(i)) / projective_nx_;  // variable
			projective_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / projective_nz_;

			unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / projective_phiBins_.at(i));
			if (projective_max_nx_ < nx)
				projective_max_nx_ = nx;
		}
		nsuperstripsPerLayer_ = projective_max_nx_ * projective_nz_;
		break;

	case SuperstripType::FOUNTAIN:
		fountain_phiBins_.clear();
		fountain_zBins_.clear();

		fountain_phiBins_.resize(16);
		fountain_zBins_.resize(16);

		for (unsigned i=0; i<phiMins_.size(); ++i) {
			if(i<6) {
				fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * fountain_sf_;
				unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
				if (fountain_max_nx_ < nx)
					fountain_max_nx_ = nx;
			}
			else	{

				if(flower_opt_==false) fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * flower_sf_;
				else			       fountain_phiBins_.at(i) = phiWidths_floweropt_.at(i) * fountain_xfactor_[i] * flower_sf_;
				unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
				if (flower_max_nx_ < nx)
					flower_max_nx_ = nx;
			}

			if(i<6) fountain_zBins_.at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;
			else    fountain_zBins_.at(i) = (zMaxs_.at(i) - zMins_.at(i)) / flower_nr_;

		}
//			if(flower_opt_ == false) fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * fountain_sf_;
//			if(flower_opt_ == true)  fountain_phiBins_.at(i) = phiWidths_floweropt_.at(i) * fountain_xfactor_[i] * fountain_sf_;
//			fountain_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;
//			unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
//			if (fountain_max_nx_ < nx)
//				fountain_max_nx_ = nx;
//		}
		//nsuperstripsPerLayer_ = fountain_max_nx_ * fountain_nz_;
		nsuperstripsPerLayer_ = 1<<12;  // 12-bit superstrip
		print();
		break;

	case SuperstripType::FOUNTAIN_FLOWER:
		fountain_phiBins_.clear();
		fountainflower_zrBins_.clear();

		fountain_phiBins_.resize(16);
		fountainflower_zrBins_.resize(16);

		//fountain_sf is for barrel layers. flower sf is for disks. Fountain_phibins contains 16 layers.
		//nz is only for barrel (first 6 places in vector), next ten are for disks nr.
		for (unsigned i=0; i<phiMins_.size(); ++i) { //phiMins is always size 16. resized when getting boundaries.
			if(i<6) {
				fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * fountain_sf_;
				unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
				if (fountain_max_nx_ < nx)
					fountain_max_nx_ = nx;
			}
			else	{

				if(flower_opt_==false) fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * flower_sf_;
				else			       fountain_phiBins_.at(i) = phiWidths_floweropt_.at(i) * fountain_xfactor_[i] * flower_sf_;
				unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
				if (flower_max_nx_ < nx)
					flower_max_nx_ = nx;
			}
//			std::cout << "phi bins " << "layer " << i << " " <<fountain_phiBins_.at(i) << std::endl;
			if(i<6) fountainflower_zrBins_.at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;
			else    fountainflower_zrBins_.at(i) = (zMaxs_.at(i) - zMins_.at(i)) / flower_nr_;

//			unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
//			if (fountain_max_nx_ < nx)
//				fountain_max_nx_ = nx;
		}

		//    	Set reference curve SSindex on hybrid (5   -> 10    ;   11  ->20)
		//    	Fill disks (other = 0)  [at(0),at(5)] ; [at(6),at(15)]
		flower_firstSS_phizero_.clear();
		flower_firstSS_phizero_.resize(16);
		flower_lastSS_phizero_.clear();
		flower_lastSS_phizero_.resize(16);

		std::cout << "arbitrer charge: "<< flower_arbitrer_charge << std::endl;
		if(flower_arbitrer_charge >= 0) {//ATTENTION! STARTS FROM 6th place!! other 0.
			for (unsigned i=6; i<phiMins_.size(); ++i){ //from layer 11 (->at(6) to layer 20 -> at (15)
				//phi0 of the reference curve that passes at phi_max @ r_max tower (FIRST SS) (anticlockwise)
				flower_firstSS_phizero_.at(i)= asin(0.3*3.8*flower_arbitrer_charge/flower_reference_pt_ / 2 * zMins_.at(i)/100) + phiMins_.at(i);  //FIRST: z:max,phimins
				flower_lastSS_phizero_.at(i)=  asin(0.3*3.8*flower_arbitrer_charge/flower_reference_pt_ / 2 * zMaxs_.at(i)/100) + phiMaxs_.at(i);    //LAST   zmin, phimax
				std::cout <<"layer " << i <<" zMax " <<zMaxs_.at(i) << " zMin " << zMins_.at(i) << " phiMax " << phiMaxs_.at(i) << " phiMin " << phiMins_.at(i) << " flower_first_phi0 " << flower_firstSS_phizero_.at(i)<< "flower_last_phi0 "<<flower_lastSS_phizero_.at(i)<<std::endl;

			}
		}
		else{
			for (unsigned i=6; i<phiMins_.size(); ++i){
				//phi0 of the reference curve that passes at phi_max @ r_min tower (FIRST SS) (clockwise)
				flower_firstSS_phizero_.at(i)= asin(0.3*3.8*flower_arbitrer_charge/flower_reference_pt_ / 2 * zMaxs_.at(i)/100) + phiMins_.at(i);  //FIRST: zmin,phimins
				flower_lastSS_phizero_.at(i)=  asin(0.3*3.8*flower_arbitrer_charge/flower_reference_pt_ / 2 * zMins_.at(i)/100) + phiMaxs_.at(i);    //LAST:   zmax, phi max
			}
		}
//		std::cout << "setting first last ss ok" << std::endl;
		nsuperstripsPerLayer_ = 1<<12;  // 12-bit superstrip
		print();

		//		# SS flower in layer
		std::cout << " [# SS-flower] (if layers with them are in tower)... "<<std::endl;
		for (unsigned q=6;q<phiMins_.size();++q){
			if (zMins_.at(q)!=0){
			n_flowerstrips = (flower_lastSS_phizero_.at(q)-flower_firstSS_phizero_.at(q))/fountain_phiBins_.at(q);
			std::cout << std::flush << "layer# "<< q+5 << ": " << n_flowerstrips << "; "
					<< "% respect to fountain " << (flower_lastSS_phizero_.at(q)-flower_firstSS_phizero_.at(q))/(phiMaxs_.at(q)-phiMins_.at(q))<<" * ";
			}
		}
		std::cout <<std::endl;
		//nsuperstripsPerLayer_ = fountain_max_nx_ * fountain_nz_;

//		std::cout << "print ok, finished constructing arbitrer " <<std::endl;
		break;

	case SuperstripType::FOUNTAINOPT:
		fountain_phiBins_.clear();
		fountain_zBins_.clear();

		fountain_phiBins_.resize(16);
		fountain_zBins_.resize(16);

		for (unsigned i=0; i<phiMins_.size(); ++i) {
			if (fountainopt_pt_ == 1)
				fountain_phiBins_.at(i) = phiWidths_opt_lowpt_.at(i);
			else
				fountain_phiBins_.at(i) = phiWidths_opt_highpt_.at(i);

			fountain_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;

			unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
			if (fountain_max_nx_ < nx)
				fountain_max_nx_ = nx;
		}
		//nsuperstripsPerLayer_ = fountain_max_nx_ * fountain_nz_;
		nsuperstripsPerLayer_ = 1<<12;  // 12-bit superstrip
		print();
		break;

	default:
		break;
	}
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripLocal(unsigned moduleId, float strip, float segment) const {
    switch (sstype_) {
    case SuperstripType::FIXEDWIDTH:
        return superstripFixedwidth(moduleId, strip, segment);
        break;

    default:
        throw std::logic_error("Incompatible superstrip type.");
        break;
    }
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripGlobal(unsigned moduleId, float r, float phi, float z, float ds) const {
    switch (sstype_) {
    case SuperstripType::PROJECTIVE:
        return superstripProjective(moduleId, r, phi, z, ds);
        break;

    case SuperstripType::FOUNTAIN:
        return superstripFountain(moduleId, r, phi, z, ds);
        break;

    case SuperstripType::FOUNTAIN_FLOWER:
        return superstripFountainFlower(moduleId, r, phi, z, ds);
        break;

    case SuperstripType::FOUNTAINOPT:
    	return superstripFountain(moduleId, r, phi, z, ds);
    	break;

    default:
        throw std::logic_error("Incompatible superstrip type.");
        break;
    }
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::compressModuleId(unsigned moduleId) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));
    //assert(lay16 < 16);

    std::vector<unsigned>::const_iterator found = std::lower_bound(towerModules_.at(lay16).begin(), towerModules_.at(lay16).end(), moduleId);
    if (found == towerModules_.at(lay16).end() || *found != moduleId)
        throw std::logic_error("Unexpected module ID.");

    return found - towerModules_.at(lay16).begin();
}


// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripFixedwidth(unsigned moduleId, float strip, float segment) const {
    // The superstrip ID is an unsigned integer with 3 fields: moduleId|segment|strip
    // The length of the 'strip' field is controlled by fixedwidth_bitshift1
    // The length of the 'segment' field is controlled by (fixedwidth_bitshift2 - fixedwidth_bitshift1)

    unsigned moduleCode = compressModuleId(moduleId);
    strip    -= 0.25;  // round down the half int precision
    segment  -= 0.25;  // round down the half int precision

    // Use strip in the least significant field
    unsigned ss = round_to_uint(strip) >> fixedwidth_bit_rshift1_;
    //std::cout << strip+0.25 << " " << round_to_uint(strip) << " " << ss << std::endl;

    // Use segment in the 2nd most sig field
    if (isPSModule(moduleId))
        ss |= ((round_to_uint(segment) >> fixedwidth_bit_rshift2_) << fixedwidth_bit_lshift1_);
    else
        ss |= (((round_to_uint(segment)*MAX_NSEGMENTS/2) >> fixedwidth_bit_rshift2_) << fixedwidth_bit_lshift1_);
    //std::cout << segment+0.25 << " " << round_to_uint(segment) << " " << ss << std::endl;

    // Use module ID in the most sig field
    ss |= (moduleCode << (fixedwidth_bit_lshift1_ + fixedwidth_bit_lshift2_));
    //std::cout << moduleId << " " << moduleCode << " " << ss << std::endl;

    return ss;
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripProjective(unsigned moduleId, float r, float phi, float z, float ds) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));

    int n_phi = projective_max_nx_;
    int n_z   = projective_nz_;

    int i_phi = std::floor((phi - phiMins_.at(lay16)) / projective_phiBins_.at(lay16));
    int i_z   = std::floor((z - zMins_.at(lay16)) / projective_zBins_.at(lay16));

    i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi) ? (n_phi - 1) : i_phi;  // proper range
    i_z       = (i_z   < 0) ? 0 : (i_z   >= n_z  ) ? (n_z   - 1) : i_z;    // proper range

    unsigned ss = i_z * n_phi + i_phi;
    return ss;
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripFountain(unsigned moduleId, float r, float phi, float z, float ds) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));

    // For barrel, correct phi based on ds and dr
    //if (lay16 < 6)
    //    phi = phi + drCorrs_.at(lay16) * ds * (r - rMeans_.at(lay16));

    int n_phi = fountain_max_nx_;
    int n_z   = fountain_nz_;
    int i_phi, i_z;


    i_phi = std::floor((phi - phiMins_.at(lay16)) / fountain_phiBins_.at(lay16));
    i_z   = std::floor((z - zMins_.at(lay16)) / fountain_zBins_.at(lay16));

    i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi) ? (n_phi - 1) : i_phi;  // proper range
    i_z       = (i_z   < 0) ? 0 : (i_z   >= n_z  ) ? (n_z   - 1) : i_z;    // proper range


    //unsigned ss = i_z * n_phi + i_phi;
    unsigned ss = ((i_z & 0x7) << 9) | (i_phi & 0x1ff);  // use magic number of 512 (= 1<<9)
    ss = ss | (lay16 << 12); //encode also first 4 bit with layer information
    return ss;
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripFountainFlower(unsigned moduleId, float r, float phi, float z, float ds) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));

    // For barrel, correct phi based on ds and dr
    //if (lay16 < 6)
    //    phi = phi + drCorrs_.at(lay16) * ds * (r - rMeans_.at(lay16));

//    08/13 - HYBRID ATTEMPT. nz == nr on disks; pT = 5; phi0 = 0; + charge (phi(rmax) > phi(rmin); B=3.8T
//    phibins is phi width hardcoded * sf
//    if FLOWER= false, just export the barrel SS definition, with hardcoded phi widths.
//    ssId are from 0 to n_phi, n_z for barrel; from (n_phi+1) elsewhere

    int n_phi = fountain_max_nx_  ; //+2 for underflow and overflow tower
    int n_phi_disk = flower_max_nx_;
    int n_z   = fountain_nz_;
    int n_r_disk = flower_nr_;
    int i_phi, i_z;

//    std::cout << "coder -- layer " << lay16 << std::endl;
    if (lay16<6){ //barrel (5-6-7-8-9-10)real aka (0,1,2,3,4,5)index lay16
    	i_phi = std::floor((phi - phiMins_.at(lay16)) / fountain_phiBins_.at(lay16));
    	i_z   = std::floor((z - zMins_.at(lay16)) / fountainflower_zrBins_.at(lay16));
//    	std::cout << " coder done, zMins" << zMins_.at(lay16)<<std::endl;
    	i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi) ? (n_phi - 1) : i_phi;  // proper range n real ss: [0 (underflow), 1,..., n, n+1 (over)]
    	i_z       = (i_z   < 0) ? 0 : (i_z   >= n_z  ) ? (n_z   - 1) : i_z;    // proper range
    }

    else { //disks (11-20) aka (7,8,...,15)
    	float phi_reference = -asin(0.3*3.8*flower_arbitrer_charge/flower_reference_pt_ / 2 * r/100) + flower_firstSS_phizero_.at(lay16);
    	i_phi = std::floor((phi - phi_reference ) / fountain_phiBins_.at(lay16));
    	i_z   = std::floor((z - zMins_.at(lay16)) / fountainflower_zrBins_.at(lay16)) ;
//    	if(i_phi < 0) std::cout << "iphi<0 -> phi" << phi << "  r " << r << " .. phi_reference " << phi_reference << std::endl;
//    	if(i_phi>n_phi)std::cout << "iphi>n_phi-> nphi: " << n_phi << " .. phi" << phi << "  r " << r << " .. phi_reference " << phi_reference << std::endl;


    	i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi_disk) ? (n_phi_disk -1) : i_phi;  // proper range
    	i_z       = (i_z   < 0) ? 0 : (i_z   >= n_r_disk  ) ? (n_r_disk  - 1) : i_z;    // proper range

    }

    //unsigned ss = i_z * n_phi + i_phi;
    unsigned ss = ((i_z & 0x7) << 9) | (i_phi & 0x1ff);  // use magic number of 512 (= 1<<9)
    ss = ss | (lay16 << 12); //encode also first 4 bit with layer information
    return ss;
}
// _____________________________________________________________________________
void SuperstripArbiter::print() {
    switch (sstype_) {
    case SuperstripType::FIXEDWIDTH:
        std::cout << "Using fixedwidth superstrip with nstrips: " << fixedwidth_nstrips_ << ", nz: " << fixedwidth_nz_
                  << ", left shifts: " << fixedwidth_bit_lshift1_ << "," << fixedwidth_bit_lshift2_
                  << ", right shifts: " << fixedwidth_bit_rshift1_ << "," << fixedwidth_bit_rshift2_
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::PROJECTIVE:
        std::cout << "Using projective superstrip with nx: " << projective_nx_ << ", nz: " << projective_nz_
                  << ", phi bins: " << projective_phiBins_.at(0) << "," << projective_phiBins_.at(1) << "," << projective_phiBins_.at(2)
                  << "," << projective_phiBins_.at(3) << "," << projective_phiBins_.at(4) << "," << projective_phiBins_.at(5)
                  << ", z bins: " << projective_zBins_.at(0) << "," << projective_zBins_.at(1) << "," << projective_zBins_.at(2)
                  << "," << projective_zBins_.at(3) << "," << projective_zBins_.at(4) << "," << projective_zBins_.at(5)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::FOUNTAIN:
        std::cout << "Using fountain superstrip [barrel] with sf: " << fountain_sf_ << ", nz: " << fountain_nz_
				  << "Using fountain superstrip  [disk]  with sf: " << flower_sf_ << ", nr: " << flower_nr_
                  << ", phi bins (barrel): " << fountain_phiBins_.at(0) << "," << fountain_phiBins_.at(1) << "," << fountain_phiBins_.at(2)
                  << "," << fountain_phiBins_.at(3) << "," << fountain_phiBins_.at(4) << "," << fountain_phiBins_.at(5)
				  << ", phi bins (disks): " << fountain_phiBins_.at(6) << "," << fountain_phiBins_.at(7) << "," << fountain_phiBins_.at(8)
				  << "," << fountain_phiBins_.at(9) << "," << fountain_phiBins_.at(10) << ","
				  << ", z bins (barrel): " << fountain_zBins_.at(0) << "," << fountain_zBins_.at(1) << "," << fountain_zBins_.at(2)
                  << "," << fountain_zBins_.at(3) << "," << fountain_zBins_.at(4) << "," << fountain_zBins_.at(5)
				  << ", z bins (disk): " << fountain_zBins_.at(6) << "," << fountain_zBins_.at(7) << "," << fountain_zBins_.at(8)
				  << "," << fountain_zBins_.at(9) << "," << fountain_zBins_.at(10)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::FOUNTAIN_FLOWER:
    	std::cout  << "Using fountain - FLOWERS ss -----> [charge]: " << flower_arbitrer_charge << "; [reference pT]: "<< flower_reference_pt_
    			<< "  [in barrel] using fountain superstrip with sf: " << fountain_sf_ << ", nz: " << fountain_nz_
				<< ", phi bins: " << fountain_phiBins_.at(0) << "," << fountain_phiBins_.at(1) << "," << fountain_phiBins_.at(2)
				<< "," << fountain_phiBins_.at(3) << "," << fountain_phiBins_.at(4) << "," << fountain_phiBins_.at(5)
				<< ", z bins: " << fountainflower_zrBins_.at(0) << "," << fountainflower_zrBins_.at(1) << "," << fountainflower_zrBins_.at(2)
				<< "," << fountainflower_zrBins_.at(3) << "," << fountainflower_zrBins_.at(4) << "," << fountainflower_zrBins_.at(5)
				<< "[in disk, only 5 positive] using flower superstrip with sf: " << flower_sf_ << ", nz: " << flower_nr_
				<< ", phi bins: " << fountain_phiBins_.at(6) << "," << fountain_phiBins_.at(7) << "," << fountain_phiBins_.at(8)
				<< "," << fountain_phiBins_.at(9) << "," << fountain_phiBins_.at(10) << ","
				<< ", R bins (if layer not present, is 0: " << fountainflower_zrBins_.at(0) << "," << fountainflower_zrBins_.at(1) << "," << fountainflower_zrBins_.at(2)
				<< "," << fountainflower_zrBins_.at(3) << "," << fountainflower_zrBins_.at(4) << "," << fountainflower_zrBins_.at(5)
				<< ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
    	break;

    case SuperstripType::FOUNTAINOPT:
        std::cout << "Using optimized pt range: " << fountainopt_pt_ << ", nz: " << fountain_nz_
                  << ", phi bins: " << fountain_phiBins_.at(0) << "," << fountain_phiBins_.at(1) << "," << fountain_phiBins_.at(2)
                  << "," << fountain_phiBins_.at(3) << "," << fountain_phiBins_.at(4) << "," << fountain_phiBins_.at(5)
                  << ", z bins: " << fountain_zBins_.at(0) << "," << fountain_zBins_.at(1) << "," << fountain_zBins_.at(2)
                  << "," << fountain_zBins_.at(3) << "," << fountain_zBins_.at(4) << "," << fountain_zBins_.at(5)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    default:
        break;
    }
}
