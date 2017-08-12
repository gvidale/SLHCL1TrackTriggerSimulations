#define layerSelectionStudy_cxx
#include "layerSelectionStudy.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TText.h>
#include <iostream>
#include <set>
#include <map>


//using namespace std;

void printSet(std::set<vector <UInt_t> > setLayers) {
	for (auto it = setLayers.begin(); it != setLayers.end(); ++it) {
		vector <UInt_t> v = *it;
		for (UInt_t i=0; i<v.size(); ++i) cout << v.at(i) << "\t";
		cout << endl;
	}
}

void layerSelectionStudy::Loop(UInt_t nEv)
{
//   In a ROOT session, you can do:
//      Root > .L layerSelectionStudy.C
//      Root > layerSelectionStudy t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1F* hNStubs = new TH1F("hNStubs","hNStubs",10,0,10);
   TH1F* hIConf    = new TH1F("hIConf","hIConf",100,0,100);
   TH1F* hGoodConf = new TH1F("hGoodConf","hGoodConf",100,0,100);

   Int_t iCompressedLayer[23] = {
  		 -1,-1,-1,-1,-1,
  		 0 , 1, 2, 3, 4, 5,
			 6 , 7, 8, 9,10,
			 -1,-1,
			 11,12,13,14,15
   };

   const UInt_t nLayers = 16;
   Bool_t aConfig[nLayers] = {false};

   const UInt_t nValidConfs = 15;
   Bool_t bValidConfigs[nLayers*nValidConfs] = {
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
   UInt_t aLay[nLayers*nValidConfs] = {
  		 5, 6, 7, 8, 9, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  		 5, 6, 7, 8, 9,  0, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  		 5, 6, 7, 8, 0,  0, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0,
  		 5, 6, 7, 0, 0,  0, 11, 12, 13,  0,  0,  0,  0,  0,  0,  0,
  		 5, 6, 0, 0, 0,  0, 11, 12, 13, 14,  0,  0,  0,  0,  0,  0,
  		 5, 6, 0, 0, 0,  0, 11, 12,  0, 14, 15,  0,  0,  0,  0,  0,
  		 5, 6, 0, 0, 0,  0, 11,  0, 13, 14, 15,  0,  0,  0,  0,  0,
  		 5, 0, 0, 0, 0,  0, 11, 12, 13, 14, 15,  0,  0,  0,  0,  0,
  		 5, 6, 7, 8, 9,  0,  0,  0,  0,  0,  0, 18,  0,  0,  0,  0,
  		 5, 6, 7, 8, 0,  0,  0,  0,  0,  0,  0, 18, 19,  0,  0,  0,
  		 5, 6, 7, 0, 0,  0,  0,  0,  0,  0,  0, 18, 19, 20,  0,  0,
  		 5, 6, 0, 0, 0,  0,  0,  0,  0,  0,  0, 18, 19, 20, 21,  0,
  		 5, 6, 0, 0, 0,  0,  0,  0,  0,  0,  0, 18, 19,  0, 21, 22,
  		 5, 6, 0, 0, 0,  0,  0,  0,  0,  0,  0, 18,  0, 20, 21, 22,
  		 5, 0, 0, 0, 0,  0,  0,  0,  0,  0,  0, 18, 19, 20, 21, 22,
   };

   std::set<vector <UInt_t> > setLayers;
   std::set<vector <UInt_t> > setValidLayers;
   for (UInt_t iConfs=0; iConfs<nValidConfs; ++iConfs) {
  	 vector <UInt_t> v;
     for (UInt_t iL=0; iL<nLayers; ++iL) {
    	 v.push_back(aLay[iL+nLayers*iConfs]);
     }
     setValidLayers.insert(v);
   }

   std::set<vector <UInt_t> > setFilteredLayers;
   std::map<vector <UInt_t>, UInt_t> mapLayers;
   std::map<vector <UInt_t>, UInt_t> mapFilteredLayers;
   std::map<vector <UInt_t>, UInt_t>::iterator it_map;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (jentry%100000==0) cout << "******** Processing event " << jentry << endl;
      if (ientry < 0 || jentry>nEv) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      hNStubs->Fill(TTStubs_modId->size());
//      if (genParts_eta->front() < 0) continue;
      if (genParts_pt ->front() < 3) continue;
      if (TTStubs_modId->size() < 6) continue;
      vector <UInt_t> vLayers;
      Bool_t bConfig[nLayers] = {false};
      UInt_t nStubs = TTStubs_modId->size();
      for (UInt_t iStub=0; iStub<nStubs; ++iStub) {
      	UInt_t iLayer = (UInt_t)TTStubs_modId->at(iStub)/10000;
      	bConfig[iCompressedLayer[iLayer]] = true;
      	vLayers.push_back(iLayer);
      }

//  		for (UInt_t i=0; i<vLayers.size(); ++i) cout << vLayers.at(i) << "\t";
//  		cout << endl;
//    	for (UInt_t iL=0; iL<nLayers; ++iL) cout << bConfig[iL] << "\t";
//  		cout << endl;

      Int_t iGoodConf = -1;
      for (UInt_t iConfs=0; iConfs<nValidConfs; ++iConfs) {
      	UInt_t nGoodLayers = 0;
      	for (UInt_t iL=0; iL<nLayers; ++iL) {
      		nGoodLayers += bConfig[iL]&&bValidConfigs[iL+nLayers*iConfs];
      		if (nGoodLayers>6) cout << "MHHH" << endl;
      	}
      	if (nGoodLayers==6) {
      		iGoodConf = iConfs;
      		break;
      	}
      }

      it_map = mapLayers.find(vLayers);
      if (it_map == mapLayers.end()) {
      	UInt_t nConfigs = mapLayers.size();
      	mapLayers.insert(std::pair<vector <UInt_t>, UInt_t>(vLayers,nConfigs));
      }
      it_map = mapLayers.find(vLayers);
      hIConf->Fill(it_map->second);
      if (iGoodConf>=0) hGoodConf->Fill(it_map->second);

      setLayers.insert(vLayers);
      std::set<vector <UInt_t> >::iterator it;
      it = setValidLayers.find(vLayers);
      if (it != setValidLayers.end()) setFilteredLayers.insert(vLayers);
   }

//   cout << "# layer configurations: " << setLayers.size() << endl;
//   printSet(setLayers);
//   cout << "# Valid configurations: " << setValidLayers.size() << endl;
//   printSet(setValidLayers);
//   cout << "# Filtered configurations: " << setFilteredLayers.size() << endl;
//   printSet(setFilteredLayers);

   TCanvas* cNStubs = new TCanvas("cNStubs","cNStubs",0,0,600,600);
   cNStubs->SetLogy();
   hNStubs->DrawCopy();

   TCanvas* cConfigs = new TCanvas("cConfigs","cConfigs",100,100,1600,600);
   cConfigs->SetLogy();
   hIConf->DrawCopy();
   hGoodConf->SetLineColor(2);
   hGoodConf->DrawCopy("same");
   for (it_map = mapLayers.begin(); it_map != mapLayers.end(); ++it_map) {
  	 TString sConf("");
  	 UInt_t iConf = it_map->second;
  	 vector <UInt_t>  v = it_map->first;
  	 for (UInt_t iL=0; iL<v.size(); ++iL) {
  		 sConf += TString(std::to_string(v.at(iL))+TString("  "));
  	 }
  	 cout << iConf << "\t\t" << sConf << endl;
  	 TText t;
  	 t.SetTextAngle(90);
  	 t.SetTextSize(0.03);
  	 t.DrawText(iConf+0.9,1,sConf);

   }
	 cout << "Good configurations' coverage: " << hGoodConf->Integral()/hIConf->Integral() << endl;
}

void layerSelectionStudyRun(UInt_t nEv=10000) {
	layerSelectionStudy l;
	l.Loop(nEv);
	return;
}
