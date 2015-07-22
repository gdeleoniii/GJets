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
  
  TreeReader data("partial400toinf.root");
  
  TH1D* h_parton  = new TH1D(" ","HT 400 to Inf",50,0,900);

  Int_t event = 0;
  
  // begin of event loop                                                                                                                                    
  for (Long64_t ev = 0; ev < /*data.GetEntriesFast()*/100000; ev++){
    
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
    
    TLorentzVector dquark(0,0,0,0);
    TLorentzVector uquark(0,0,0,0);
    TLorentzVector squark(0,0,0,0);
    TLorentzVector cquark(0,0,0,0);
    TLorentzVector bquark(0,0,0,0);
    TLorentzVector gluon(0,0,0,0);
    
    Int_t photonIndex = -1;
    vector<int> jetIndices;  
    
    for(int gen=0; gen < nGenPar; gen++) {
      if(genParSt[gen] == 3 && genParId[gen] == 22) photonIndex = gen; //conditions for photon
      else if(genParSt[gen] == 3 && ( (abs(genParId[gen]) < 6 && abs(genParId[gen]) > 0) || genParId[gen] == 21)) //conditions for gluons or jets 
	jetIndices.push_back(gen);
    }
    
    if(photonIndex < 0) continue; //for protection
    
    Int_t phoMo1 = genMo1[photonIndex]; //declaring mother for photon
    Int_t phoMo2 = genMo2[photonIndex]; //declaring mother for photon
    //genMomParId is used only when 1 particle decayed, for collision and likes, we use genMo1 and genMo2    
    
    
    for(int ij=0; ij < jetIndices.size(); ij++) {
      Int_t jetIndex = jetIndices[ij];
      if(jetIndex < 0) continue; //for protection
      
      Int_t jetMo1 = genMo1[jetIndex]; //declaring mother for b jets
      Int_t jetMo2 = genMo2[jetIndex]; //declaring mother for b jets
      
      if(phoMo1 != jetMo1 || phoMo2 != jetMo2) continue; //checking for same mother particle

      if(genParId[jetIndex] == 21) gluon.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      else if(abs(genParId[jetIndex]) == 1)
	dquark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      else if(abs(genParId[jetIndex]) == 2)
        uquark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      else if(abs(genParId[jetIndex]) == 3)
        cquark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      else if(abs(genParId[jetIndex]) == 4)
        squark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
      else if(abs(genParId[jetIndex]) == 5)
        bquark.SetPtEtaPhiM(genParPt[jetIndex],genParEta[jetIndex],genParPhi[jetIndex],genParM[jetIndex]);
       
      Float_t p = sqrt(gluon.Pt()*gluon.Pt()) + sqrt(dquark.Pt()*dquark.Pt()) + sqrt(squark.Pt()*squark.Pt()) + sqrt(uquark.Pt()*uquark.Pt()) + sqrt(cquark.Pt()*cquark.Pt()) + sqrt(bquark.Pt()*bquark.Pt());
      h_parton->Fill(p);
      event++;
    }
    
    
  } // end of event loop 
  h_parton->SetLineColor(kRed-1);
  h_parton->SetFillColor(kRed-1);
  h_parton->GetXaxis()->SetTitle("H_{T} (GeV)");
  h_parton->Draw();
  std::cout << "HT = " << event << std::endl;
}

