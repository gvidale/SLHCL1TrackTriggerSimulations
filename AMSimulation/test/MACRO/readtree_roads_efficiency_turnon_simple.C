/*
Road efficiency & turn on curves

19/07/2017

Giorgio Vidale

 */

#define readtree_roads_cxx
#include "readtree_roads.h"
#include <TH1.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include "functions.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/StubCleaner.h"

void readtree_roads::Loop(TString key)
{

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Int_t evt_not_reconstructed=0;

	Int_t roads_before_reconstruction=0;
	Float_t average_roads_before_reconstruction=0;
	Float_t average_true_mu_perevent=0;
	Float_t average_recognized_mutrack_perevent=0;
	Float_t average_recognized_mutrack=0;
	Float_t average_roads_fired=0;

	Int_t genuine_mu=0;
	Int_t all_stub=0;
	Int_t n_roads_fired;
	Int_t n_evt_noroads=0;
	Int_t n_stubs;
	Int_t trk_index;
	unsigned int stub_index;

	Int_t * layer_phi_zeta;
	Float_t trkParts_ChargeOverPt;

	set<Int_t> setOflayers;		//remove overlap in same layer
	pair < set<Int_t>::iterator , Bool_t>  insertYES;


//	Histo with mu tracking particle (id0) attributes count. Confront counts of all roads with counts from roads with at least 4 stubs form mu.

	TH1F * mu_trk_phi = new TH1F("mu_trk_phi","mu_trk_phi",20,0.60,1.75);
	mu_trk_phi->Sumw2();
	mu_trk_phi->SetTitle("'true TT25' mu-trkpart phi count; phi; parts count ");

	TH1F * mu_trk_eta = new TH1F("mu_trk_eta","mu_trk_eta", 20,-0.5,1.5);
	mu_trk_eta->Sumw2();
	mu_trk_eta->SetTitle("'true TT25' mu-trkpart eta count ; eta; parts count");

	TH1F * mu_trk_pT = new TH1F("mu_trk_pT","mu_trk_pT",200,0,200);
	mu_trk_pT->Sumw2();
	mu_trk_pT->SetTitle("'true TT25' mu-trkpart pT count; trkpart_pT; parts count");

	TH1F * mu_trk_phi_4 = new TH1F("mu_trk_phi_4","mu_trk_phi_4",20,0.60,1.75);
	mu_trk_phi_4->Sumw2();
	mu_trk_phi_4->SetTitle("road efficiency (TT25 @5 stub) ; phi; % efficiency");

	TH1F * mu_trk_eta_4 = new TH1F("mu_trk_eta_4","mu_trk_eta_4", 20,-0.5,1.5);
	mu_trk_eta_4->Sumw2();
	mu_trk_eta_4->SetTitle("road efficiency (TT25 @5 stub) ; eta; % efficiency");

	TH1F * mu_trk_pT_4 = new TH1F("mu_trk_pT_4","mu_trk_pT_4",200,0,200);
	mu_trk_pT_4->Sumw2();
	mu_trk_pT_4->SetTitle("road efficiency (TT25 @5 stub) ; pT; % efficiency");

	TH1I * road_read_before_reconstruction = new TH1I("road_read_before_reconstruction","road_read_before_reconstruction",1500,0,1500);

//	*****

	Long64_t nbytes = 0, nb = 0;

//	LOOP EVENTS

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		if (jentry % 500 == 0) cout << "@@@ event " << jentry << " ... Roads fired 'til now: " << average_roads_fired << endl;

		roads_before_reconstruction = 0;

//		APPLY PHASE SPACE CUT  (trkPart 0)
		trkParts_ChargeOverPt = float(trkParts_charge->at(0))/float(trkParts_pt->at(0));
		int aux_TT = TrackParametersToTT().get_tt(trkParts_phi->at(0), trkParts_ChargeOverPt,trkParts_eta->at(0),trkParts_vz->at(0));
		if (aux_TT != 25) continue;
//		*****


		average_true_mu_perevent++; //average number of true mu that pass the filters.

		mu_trk_phi->Fill(trkParts_phi->at(0));
		mu_trk_pT->Fill(trkParts_pt->at(0));
		mu_trk_eta->Fill(trkParts_eta->at(0));


//		*****

//		LOOP OVER PATTERN FIRED
		n_roads_fired=AMTTRoads_patternRef->size();
		if(n_roads_fired<1) {
			++n_evt_noroads;	//# of roads that has not fired
			continue; //jump to other event if this one hasn't fired.
		}

		for (Int_t patt_count=0; patt_count<n_roads_fired; ++patt_count){

			genuine_mu = 0;

			//LOOP OVER layers for each pattern
			Int_t l=(*AMTTRoads_stubRefs)[patt_count].size();
			assert(l==6);
			for (l =0; l<6; ++l){

				//LOOP OVER STUBS LIST for that layer, in that pattern, for that event.
				n_stubs=(*AMTTRoads_stubRefs)[patt_count][l].size();
				for (Int_t stub_count = 0; stub_count< n_stubs; ++stub_count){

					stub_index=(*AMTTRoads_stubRefs)[patt_count][l][stub_count];

					trk_index=TTStubs_tpId->at(stub_index);

//					if that trkpart is trk 0, go ahead
					if( trk_index!=0) continue;
//
					layer_phi_zeta = phi_zeta(TTStubs_modId->at(stub_index)); //extract layer info //a[0] is the layer
					insertYES = setOflayers.insert(layer_phi_zeta[0]);   //  do not overcount overlap;

					if(insertYES.second) {
						genuine_mu++;
						//for that particular mu-trk, count the number of stubs that leaves in the road
					}
				}
			}
//			keep trace of how many roads I need to throw out AM to get a reconstruction candidate.
			roads_before_reconstruction++;

//			ROAD EFFICIENCY 4 - 5 - 6 Stub?

			if (genuine_mu>=5){
				mu_trk_phi_4->Fill(trkParts_phi->at(0));
				mu_trk_pT_4->Fill(trkParts_pt->at(0));
				mu_trk_eta_4->Fill(trkParts_eta->at(0));

				average_recognized_mutrack_perevent++;
				road_read_before_reconstruction->Fill(roads_before_reconstruction);
				average_roads_before_reconstruction+=roads_before_reconstruction;

				break;
			}

			if(patt_count==n_roads_fired-1){
				evt_not_reconstructed++;
			}
		}

//		Erase the setOfLayer, before checking new road
		setOflayers.erase(setOflayers.begin(),setOflayers.end());


//		SOME AVERAGE ESTIMATES til now...
		average_roads_fired+=n_roads_fired;
		average_roads_before_reconstruction+=roads_before_reconstruction;




	}

	//FINAL AVERAGE QUANTITIES
	average_roads_fired=average_roads_fired/average_true_mu_perevent;
	average_roads_before_reconstruction=average_roads_before_reconstruction/nentries;
	average_recognized_mutrack_perevent=average_recognized_mutrack_perevent/average_true_mu_perevent;
	average_true_mu_perevent = average_true_mu_perevent/nentries;

	cout << endl;
	cout << "*************" << endl;
	cout << "total events read: " << nentries << endl;
	cout << "evt not reconstructed " <<evt_not_reconstructed << endl;
	cout << "... average of " << average_roads_fired << " fired roads per event" << endl;
	cout << "... average of " << average_roads_before_reconstruction << " roads before reconstruction" << endl;
	cout << "... average of " << average_true_mu_perevent << " true mu track per event" <<endl;
	cout << "... average of " << average_recognized_mutrack_perevent << " mu track recognized per event" << endl;

	TCanvas * c[3];
	for(Int_t k=0; k<4 ;++k){
		c[k]=new TCanvas();
	}


	//   Plot cross check turn on curves. Number of roads / event that has at least for mu-stub in them
	c[0]->Divide(1,2);
	c[0]->cd(1);
	mu_trk_phi->Draw();
	mu_trk_phi_4->Divide(mu_trk_phi); //ratio as per title
	c[0]->cd(2);
	mu_trk_phi_4->Draw();

	c[1]->Divide(1,2);
	c[1]->cd(1);
	mu_trk_pT->Draw();
	mu_trk_pT_4->Divide(mu_trk_pT);
	c[1]->cd(2);
	mu_trk_pT_4->Draw();

	c[2]->Divide(1,2);
	c[2]->cd(1);
	mu_trk_eta->Draw();
	mu_trk_eta_4->Divide(mu_trk_eta);
	c[2]->cd(2);
	mu_trk_eta_4->Draw();

}



//      **********************
//               MAIN
//      **********************

void readtree_roads_efficiency_turnon_simple(){


	gStyle->SetOptStat(111111);
	TString key;
	TString fName;
	TString name;



	name="roads_mu_PU200_8K_sf1_nz1_pt3_0718"; //PU200 (8400 evt)
//	name="roads_mu_PU0_TT25_50M_sf1_nz1_0719"; //only mu ((10000 evt)

	fName = "../roads/"+name+".root"; //name of file root to read, with path
	key = "_"+name;     //just a key for histos

//	HERE I DEFINE OUTPUT .ROOT WITH THE HISTOGRAMS TO CHECK

	bool isfOpen;
	TFile* f = 0; isfOpen = false;
	f = new TFile(name+"_efficiency_trk0_PU200.root","RECREATE");
	isfOpen = f->IsOpen();
	if (!isfOpen) {
		cout << "ERROR. Not able to load the confrontoBranches file. Exiting..." << endl;
		return;
	}

	readtree_roads a(fName);
	f -> cd();
	a.Loop(key);

	f -> Write();

	cout << "ho scritto" << endl;

	return;


}

