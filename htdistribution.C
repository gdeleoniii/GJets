#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TH1F.h>
#include <TH1D.h>
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
  
  TreeReader data("gjets40to100_2.root");
  
  TH1D* h_parton  = new TH1D(" ","HT distribution",16,30,110);

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
    Float_t  HT          = data.GetFloat("HT");    
    
    Int_t photonIndex = -1;
    vector<int> jetIndices;  
    
    for(int gen=0; gen < nGenPar; gen++) {
      if(genParSt[gen] == 3 && genParId[gen] == 22) photonIndex = gen;
      else if(genParSt[gen] == 3 && ( (abs(genParId[gen]) < 6 && abs(genParId[gen]) > 0) || genParId[gen] == 21))
	jetIndices.push_back(gen);
    }
    
    if(photonIndex < 0) continue; 
    
    Int_t phoMo1 = genMo1[photonIndex];
    Int_t phoMo2 = genMo2[photonIndex];
    
    Float_t sum = HT;
    for(int ij=0; ij < jetIndices.size(); ij++) {
      Int_t jetIndex = jetIndices[ij];
      if(jetIndex < 0) continue;
      
      Int_t jetMo1 = genMo1[jetIndex]; 
      Int_t jetMo2 = genMo2[jetIndex];
     
      if(phoMo1 != jetMo1 || phoMo2 != jetMo2) continue;

    }
    h_parton->Fill(sum);
    
  } // end of event loop 
  Double_t scale = 20930/h_parton->Integral(); // 20930 is the cross section of HT 40 to 100
  h_parton->Scale(scale);
  h_parton->SetLineColor(kOrange-1); 
  h_parton->SetFillColor(kOrange-1);
  h_parton->GetXaxis()->SetTitle("H_{T} (GeV)");
  h_parton->Draw();
}

