#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "untuplizer.h"

void counter(){
  
  TreeReader data("partial200to400.root");
  
  TH1D* h_photon1 = new TH1D(" ", " Photon",50,0,350);
  TH1D* h_photon2 = new TH1D(" ", " ",50,0,350);
  TH1D* h_bjet    = new TH1D(" ", " Events with #gamma and b jets",50,0,350);
  TH1D* h_bbjet    = new TH1D(" ", "Events with #gamma and at least 2 b jets ",50,0,350);
   
  Int_t event1 = 0;
  Int_t event2 = 0;

  // begin of event loop                                                                                                                                    
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){
    
    // print progress                                                                                                                                       
    if ( ev % 50000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev+1, data.GetEntriesFast());
    
    data.GetEntry(ev);
    
    
    Int_t    nGenPar     = data.GetInt("nGenPar");
    Int_t*   genParId    = data.GetPtrInt("genParId");
    Int_t*   genParSt    = data.GetPtrInt("genParSt");
    Int_t*   genMo1      = data.GetPtrInt("genMo1");
    Int_t*   genMo2      = data.GetPtrInt("genMo2");
    Float_t* genParPt    = data.GetPtrFloat("genParPt");
    Float_t* genParEta   = data.GetPtrFloat("genParEta");
    Float_t* genParPhi   = data.GetPtrFloat("genParPhi");
    Float_t* genParM     = data.GetPtrFloat("genParM");
    
    TLorentzVector bquark(0,0,0,0); 
    TLorentzVector gamma(0,0,0,0);
    
    Int_t photonIndex = -1;
    vector<int> jetIndices;  

    for(int gen=0; gen < nGenPar; gen++) {
      if(genParSt[gen] == 3 && genParId[gen] == 22) photonIndex = gen;
      else if(genParSt[gen] == 3 && ( (abs(genParId[gen]) < 6 && abs(genParId[gen]) > 0) || genParId[gen] == 22)) 
	jetIndices.push_back(gen);
    }
    
    if(photonIndex < 0) continue; 
    
    Int_t phoMo1 = genMo1[photonIndex];
    Int_t phoMo2 = genMo2[photonIndex];

    gamma.SetPtEtaPhiM(genParPt[photonIndex],genParEta[photonIndex],genParPhi[photonIndex],0); 
    h_photon1->Fill(gamma.Pt());
    h_photon2->Fill(gamma.Pt());

    Int_t bcount   =  0;
    
    for(int ij=0; ij < jetIndices.size(); ij++) {
      Int_t jetIndex = jetIndices[ij];
      if(jetIndex < 0) continue;
      
      Int_t jetMo1 = genMo1[jetIndex];
      Int_t jetMo2 = genMo2[jetIndex]; 
      
      if(phoMo1 != jetMo1 || phoMo2 != jetMo2) continue;
      
      if( abs(genParId[jetIndex]) != 5) continue;
     
      bquark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      bcount++;
    
    }

      if(bcount > 0) { 
	h_bjet->Fill(bquark.Pt());
	event1++;
      }
      if(bcount > 1) { 
	h_bbjet->Fill(bquark.Pt());
	event2++;
      }
    
  } // end of event loop 
  std::cout<< "atleast 1 b jets= " << event1 << std::endl;
  std::cout<< "atleast 2 b jets= " << event2 << std::endl;

  TCanvas *c = new TCanvas("c","c",900, 900);
  c->Divide(2,2);
  c->cd(1);
  h_bjet->SetLineColor(kYellow);
  h_bjet->SetLineWidth(2);
  h_bjet->GetXaxis()->SetTitle("p_{T}^{#gamma}(GeV) ");
  h_bjet->GetYaxis()->SetTitle("Entries");
  h_bbjet->SetLineColor(kBlue);
  h_bbjet->SetLineWidth(2);
  h_bbjet->GetXaxis()->SetTitle("p_{T}^{#gamma}(GeV) ");
  h_bbjet->GetYaxis()->SetTitle("Entries");
  h_bbjet->Draw();
  h_bbjet->Draw("same");
  TLegend *legend = new TLegend(0.65,0.2,0.8,0.3);
  legend->AddEntry(h_tada,"at least 1 b jet","l");
  legend->AddEntry(h_yada,"at least 2 b jets", "l");
  legend->SetFillColor(0);
  legend->Draw();

  c->cd(2);
  h_photon1->SetLineColor(kRed);
  h_photon1->SetLineWidth(2);
  h_photon1->GetXaxis()->SetTitle("p_{T}^{#gamma}(GeV) ");
  h_photon1->GetYaxis()->SetTitle("Entries");
  h_photon1->Draw();

  c->cd(3);
  TH1D *ratio1 = new TH1D(" ", "Events with #gamma and at least 1 b jet / total number of events ",50,0,350);
  ratio1->Divide(h_bjet,h_photon1);
  ratio1->SetLineWidth(2);
  ratio1->SetLineColor(kYellow);
  ratio1->GetXaxis()->SetTitle("p_{T}^{#gamma}(GeV)");
  ratio1->GetYaxis()->SetTitle("Entries");
  ratio1->Draw();

   c->cd(4);
  TH1D *ratio2 = new TH1D(" ", "Events with #gamma and at least 2 b jets / total number of events ",50,0,350);
  ratio2->Divide(h_bbjet,h_photon2);
  ratio2->SetLineWidth(2);
  ratio2->SetLineColor(kBlue);
  ratio2->GetXaxis()->SetTitle("p_{T}^{#gamma}(GeV)");
  ratio2->GetYaxis()->SetTitle("Entries");
  ratio2->Draw();
}

