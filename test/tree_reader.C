#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include <iomanip>

float AssignEtaBin(float eta){

  float ref_eta = -9999.;

  for (float m = -2.50; m < 2.45; m+=0.05){
    if (eta - m >= 0. && eta - m < 0.05){
      ref_eta = roundf(m*100)/100; //Round to two decimal places
    }
  }
  if (ref_eta == -9999.){
    std::cout << "eta: " << eta << "  ref_eta: " << ref_eta << std::endl;
    std::cout << "SOMETHING WENT WRONG" << std::endl;
  }
  return ref_eta;
}

void Prepare1DHisto(TH1F* histo[2][2], int B_or_E, const char* title_x, const char* title_y){

  float DATA_normalization = histo[B_or_E][0]->Integral();
  float MC_normalization = histo[B_or_E][1]->Integral();
  float norm_factor = DATA_normalization/MC_normalization;

  histo[B_or_E][1]->Scale(norm_factor);
  histo[B_or_E][1]->SetLineColor(800);
  histo[B_or_E][1]->SetFillColor(800);
  histo[B_or_E][1]->GetXaxis()->SetTitle(title_x);
  histo[B_or_E][1]->GetXaxis()->SetTitleSize(0.05);
  histo[B_or_E][1]->GetXaxis()->SetTitleOffset(0.7);
  histo[B_or_E][1]->GetYaxis()->SetTitle(title_y);
  histo[B_or_E][1]->GetYaxis()->SetTitleOffset(1.5);
  TH1F* h_MCErr = (TH1F*)histo[B_or_E][1]->Clone(); //Plotting MC uncertainties
  h_MCErr->SetFillStyle(3008);
  h_MCErr->SetMarkerStyle(1);
  h_MCErr->SetFillColor(kOrange+7);
  histo[B_or_E][0]->SetMarkerStyle(20); //Data

  histo[B_or_E][1]->Draw("hist");
  h_MCErr->Draw("SAME,E2");
  histo[B_or_E][0]->Draw("SAME,PE");

  histo[B_or_E][0]->Write();
  histo[B_or_E][1]->Write();
 
}

void Prepare2DHisto(TH2F* histo2D[2], int D_or_MC, const char* title_x, const char* title_y){

  histo2D[D_or_MC]->GetXaxis()->SetTitle(title_x);
  histo2D[D_or_MC]->GetYaxis()->SetTitle(title_y);
  histo2D[D_or_MC]->Draw("COLZ");
  histo2D[D_or_MC]->Write();

}

void tree_reader() {
    
  //TFile* fileIn_DATA = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/Run2017B.root");
  TFile* fileIn_DATA = new TFile("../ntuples/2017_DATA.root");
  TFile* fileIn_MC   = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/DY_LO.root");
  
  TTree* treeIn_DATA = (TTree*)fileIn_DATA->Get("tnpEleTrig/fitter_tree");
  TTree* treeIn_MC   = (TTree*)fileIn_MC->Get("tnpEleTrig/fitter_tree");
  
  TTree* tree_array[2] = {treeIn_DATA, treeIn_MC};
  
  TH1F* h_fBrem_B_DATA = new TH1F("h_fBrem_B_DATA","",25,0.,1.);
  TH1F* h_fBrem_E_DATA = new TH1F("h_fBrem_E_DATA","",25,0.,1.);
  TH1F* h_fBrem_B_MC   = new TH1F("h_fBrem_B_MC","",25,0.,1.);
  TH1F* h_fBrem_E_MC   = new TH1F("h_fBrem_E_MC","",25,0.,1.);
  
  TH1F* h_fBrem_Plus_B_DATA = new TH1F("h_fBrem_Plus_B_DATA","",25,0.,1.);
  TH1F* h_fBrem_Plus_E_DATA = new TH1F("h_fBrem_Plus_E_DATA","",25,0.,1.);
  TH1F* h_fBrem_Plus_B_MC   = new TH1F("h_fBrem_Plus_B_MC","",25,0.,1.);
  TH1F* h_fBrem_Plus_E_MC   = new TH1F("h_fBrem_Plus_E_MC","",25,0.,1.);
  
  TH1F* h_fBrem_Minus_B_DATA = new TH1F("h_fBrem_Minus_B_DATA","",25,0.,1.);
  TH1F* h_fBrem_Minus_E_DATA = new TH1F("h_fBrem_Minus_E_DATA","",25,0.,1.);
  TH1F* h_fBrem_Minus_B_MC   = new TH1F("h_fBrem_Minus_B_MC","",25,0.,1.);
  TH1F* h_fBrem_Minus_E_MC   = new TH1F("h_fBrem_Minus_E_MC","",25,0.,1.);
  
  TH1F* h_mass_B_DATA = new TH1F("h_mass_B_DATA","",30,60.,120.);
  TH1F* h_mass_E_DATA = new TH1F("h_mass_E_DATA","",30,60.,120.);
  TH1F* h_mass_B_MC   = new TH1F("h_mass_B_MC","",30,60.,120.);
  TH1F* h_mass_E_MC   = new TH1F("h_mass_E_MC","",30,60.,120.);
  
  TH2F* h2_fBrem_phi_DATA = new TH2F("h2_fBrem_phi_DATA","",100,-3.14,3.14,50,0.,1.);
  TH2F* h2_fBrem_phi_MC   = new TH2F("h2_fBrem_phi_MC","",100,-3.14,3.14,50,0.,1.);
  TH2F* h2_fBrem_pT_DATA  = new TH2F("h2_fBrem_pT_DATA","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_pT_MC    = new TH2F("h2_fBrem_pT_MC","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_eta_DATA = new TH2F("h2_fBrem_eta_DATA","",100,-2.5,2.5,50,0.,1.);
  TH2F* h2_fBrem_eta_MC   = new TH2F("h2_fBrem_eta_MC","",100,-2.5,2.5,50,0.,1.);

  TH1F* h_array_fBrem[2][2]       = {{h_fBrem_B_DATA,h_fBrem_B_MC},{h_fBrem_E_DATA,h_fBrem_E_MC}};
  TH1F* h_array_fBrem_Plus[2][2]  = {{h_fBrem_Plus_B_DATA,h_fBrem_Plus_B_MC},{h_fBrem_Plus_E_DATA,h_fBrem_Plus_E_MC}};
  TH1F* h_array_fBrem_Minus[2][2] = {{h_fBrem_Minus_B_DATA,h_fBrem_Minus_B_MC},{h_fBrem_Minus_E_DATA,h_fBrem_Minus_E_MC}};
  TH1F* h_array_mass[2][2]        = {{h_mass_B_DATA,h_mass_B_MC},{h_mass_E_DATA,h_mass_E_MC}};
  TH2F* h2_array_fBrem_phi[2]     = {h2_fBrem_phi_DATA,h2_fBrem_phi_MC};
  TH2F* h2_array_fBrem_pT[2]      = {h2_fBrem_pT_DATA,h2_fBrem_pT_MC};
  TH2F* h2_array_fBrem_eta[2]     = {h2_fBrem_eta_DATA,h2_fBrem_eta_MC};

  //TH1F* h_pT_DATA  = new TH1F("h_pT_DATA","",47,33.,80.);
  //TH1I* h_nPV_DATA = new TH1I("h_nPV_DATA","",60,0.,60.);
  
  float i_range = -2.50;
  float iplus_range = -2.45;
  map<std::string,TH1F> myMap_DATA;
  map<std::string,TH1F> myMap_MC;
  
  for (int d = 0; d < 100; d++){
    i_range = roundf(i_range*100)/100;
    iplus_range = roundf(iplus_range*100)/100;
    myMap_DATA[std::to_string(i_range)] = TH1F(TString::Format("DATA_h_fBrem_etaBins_%f", i_range) + TString::Format("_%f", iplus_range),"",50,0.,1.);
    myMap_MC[std::to_string(i_range)]   = TH1F(TString::Format("MC_h_fBrem_etaBins_%f", i_range) + TString::Format("_%f", iplus_range),"",50,0.,1.);
    i_range += 0.05;
    iplus_range += 0.05;
  }
  
  for (int D_or_MC = 0; D_or_MC < 2; D_or_MC++){
    tree_array[D_or_MC]->SetBranchStatus("*",0);
    tree_array[D_or_MC]->SetBranchStatus("passHltEle32DoubleEGWPTightGsf",1);
    tree_array[D_or_MC]->SetBranchStatus("passEGL1SingleEGOr",1);
    tree_array[D_or_MC]->SetBranchStatus("passingMVA94Xwp80isoV2",1);
    tree_array[D_or_MC]->SetBranchStatus("el_pt",1);
    tree_array[D_or_MC]->SetBranchStatus("el_eta",1);
    tree_array[D_or_MC]->SetBranchStatus("el_phi",1);
    tree_array[D_or_MC]->SetBranchStatus("tag_Ele_pt",1);
    tree_array[D_or_MC]->SetBranchStatus("tag_sc_eta",1);
    tree_array[D_or_MC]->SetBranchStatus("mass",1);
    tree_array[D_or_MC]->SetBranchStatus("el_fbrem",1);
    tree_array[D_or_MC]->SetBranchStatus("el_q",1);
    tree_array[D_or_MC]->SetBranchStatus("el_dz",1);
    tree_array[D_or_MC]->SetBranchStatus("event_nPV",1);
    if (D_or_MC == 1) tree_array[D_or_MC]->SetBranchStatus("totWeight",1);
    
    int v_passHltEle32;
    int v_passEGL1;
    int v_passMVA;
    int v_nPV;
    float v_el_pt;
    float v_el_eta;
    float v_el_phi;
    float v_el_fbrem;
    float v_el_q;
    float v_el_dz;
    float v_mass;
    float v_Tag_el_pt;
    float v_Tag_sc_eta;
    float v_totWeight;
    
    tree_array[D_or_MC]->SetBranchAddress("passHltEle32DoubleEGWPTightGsf",&v_passHltEle32);
    tree_array[D_or_MC]->SetBranchAddress("passEGL1SingleEGOr",&v_passEGL1);
    tree_array[D_or_MC]->SetBranchAddress("passingMVA94Xwp80isoV2",&v_passMVA);
    tree_array[D_or_MC]->SetBranchAddress("event_nPV",&v_nPV);
    tree_array[D_or_MC]->SetBranchAddress("el_pt",&v_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("el_eta",&v_el_eta);
    tree_array[D_or_MC]->SetBranchAddress("el_phi",&v_el_phi);
    tree_array[D_or_MC]->SetBranchAddress("el_fbrem",&v_el_fbrem);
    tree_array[D_or_MC]->SetBranchAddress("el_q",&v_el_q);
    tree_array[D_or_MC]->SetBranchAddress("el_dz",&v_el_dz);
    tree_array[D_or_MC]->SetBranchAddress("tag_Ele_pt",&v_Tag_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("tag_sc_eta",&v_Tag_sc_eta);
    tree_array[D_or_MC]->SetBranchAddress("mass",&v_mass);
    if (D_or_MC == 1) tree_array[D_or_MC]->SetBranchAddress("totWeight",&v_totWeight);
    
    for (Long64_t n = 0; n < tree_array[D_or_MC]->GetEntriesFast(); n++){
      tree_array[D_or_MC]->GetEntry(n);
      
      if (!(v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passMVA == 1 && v_el_pt > 33. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      //if (!(v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passMVA == 1 && v_el_pt > 28. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      //if not (entry.passHltEle27WPTightGsf and entry.passingMVA94Xwp80isoV2 and entry.el_pt > 28 and math.fabs(entry.el_eta) < \ 2.50) continue;                                                                                                                  
      if (!(v_Tag_el_pt > 35. && fabs(v_Tag_sc_eta) <= 2.1)) continue;
      
      if (!(v_mass > 60. && v_mass < 120.)) continue;
      
      if (!(v_el_fbrem > 0.)) continue;

      //if (!(v_el_pt < 43.)) continue;
      if (!(v_nPV > 23.)) continue;
      
      //Assign the right eta bin
      float eta_bin = AssignEtaBin(v_el_eta);
      if (eta_bin == -9999.) break;
      
      float EventWeight = 1.;
      if (D_or_MC == 0){
	myMap_DATA.find(std::to_string(eta_bin))->second.Fill(v_el_fbrem);
      }
      else{
	EventWeight = v_totWeight;
	myMap_MC.find(std::to_string(eta_bin))->second.Fill(v_el_fbrem,EventWeight);
      }

      int B_or_E = -1;
      if (fabs(v_el_eta) < 1.2) B_or_E = 0; //Corresponds to BARREL
      else B_or_E = 1; //Corresponds to ENDCAPS
      
      h_array_fBrem[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      h_array_mass[B_or_E][D_or_MC]->Fill(v_mass,EventWeight);
      if (v_el_q == 1) h_array_fBrem_Plus[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      if (v_el_q == -1) h_array_fBrem_Minus[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      h2_array_fBrem_phi[D_or_MC]->Fill(v_el_phi,v_el_fbrem,EventWeight);      
      h2_array_fBrem_pT[D_or_MC]->Fill(v_el_pt,v_el_fbrem,EventWeight);      
      h2_array_fBrem_eta[D_or_MC]->Fill(v_el_eta,v_el_fbrem,EventWeight);      
    }
  }
  
  gStyle->SetOptStat(0);
  TFile* fileOut_plots = new TFile("histos/cpp/plots_highPU.root","RECREATE");
  fileOut_plots->cd();

  TCanvas* c_fBrem_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem, 0, "f_{brem}", "Events/0.04"); // 0 stands for BARREL
  c_fBrem_B->SaveAs("plots/cpp/fBrem_B_highPU.pdf");
  
  TCanvas* c_fBrem_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem, 1, "f_{brem}", "Events/0.04"); // 1 stands for BARREL
  c_fBrem_E->SaveAs("plots/cpp/fBrem_E_highPU.pdf");
  
  TCanvas* c_fBrem_Plus_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Plus, 0, "f_{brem}", "Events/0.04");
  c_fBrem_Plus_B->SaveAs("plots/cpp/fBrem_Plus_B_highPU.pdf");
  
  TCanvas* c_fBrem_Minus_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Minus, 0, "f_{brem}", "Events/0.04");
  c_fBrem_Minus_B->SaveAs("plots/cpp/fBrem_Minus_B_highPU.pdf");
  
  TCanvas* c_fBrem_Plus_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Plus, 1, "f_{brem}", "Events/0.04");
  c_fBrem_Plus_E->SaveAs("plots/cpp/fBrem_Plus_E_highPU.pdf");
  
  TCanvas* c_fBrem_Minus_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Minus, 1, "f_{brem}", "Events/0.04");
  c_fBrem_Minus_E->SaveAs("plots/cpp/fBrem_Minus_E_highPU.pdf");
  
  TCanvas* c_mass_B = new TCanvas();
  Prepare1DHisto(h_array_mass, 0, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_B->SaveAs("plots/cpp/mass_B_highPU.pdf");
  
  TCanvas* c_mass_E = new TCanvas();
  Prepare1DHisto(h_array_mass, 1, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_E->SaveAs("plots/cpp/mass_E_highPU.pdf");
  
  TCanvas* c_fBrem_vs_phi_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_phi, 0, "#phi^{e}", "f_{brem}"); // 0 stands for DATA
  c_fBrem_vs_phi_DATA->SaveAs("plots/cpp/fBrem_vs_phi_DATA_highPU.pdf");
  
  TCanvas* c_fBrem_vs_phi_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_phi, 1, "#phi^{e}", "f_{brem}"); // 1 stands for MC
  c_fBrem_vs_phi_MC->SaveAs("plots/cpp/fBrem_vs_phi_MC_highPU.pdf");
  
  TCanvas* c_fBrem_vs_pT_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_pT, 0, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  c_fBrem_vs_pT_DATA->SaveAs("plots/cpp/fBrem_vs_pT_DATA_highPU.pdf");
  
  TCanvas* c_fBrem_vs_pT_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_pT, 1, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  c_fBrem_vs_pT_MC->SaveAs("plots/cpp/fBrem_vs_pT_MC_highPU.pdf");
  
  TCanvas* c_fBrem_vs_eta_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_eta, 0, "#eta^{e}", "f_{brem}");
  c_fBrem_vs_eta_DATA->SaveAs("plots/cpp/fBrem_vs_eta_DATA_highPU.pdf");
  
  TCanvas* c_fBrem_vs_eta_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_eta, 1, "#eta^{e}", "f_{brem}");
  c_fBrem_vs_eta_MC->SaveAs("plots/cpp/fBrem_vs_eta_MC_highPU.pdf");

  //h_pT_DATA->Write();
  //h_nPV_DATA->Write();
  
  TFile* fileOut = new TFile("histos/cpp/fBrem_notFolded_full2017_highPU.root","RECREATE");
  fileOut->cd();
  gDirectory->mkdir("DATA");
  fileOut->cd("DATA");
  for (float i = -2.50; i < 2.45; i+=0.05){
    //std::cout << std::to_string(i) << std::endl;
    myMap_DATA.find(std::to_string(roundf(i*100)/100))->second.Write();
  }
  // std::cout << "ciao" << std::endl;
  fileOut->cd();
  gDirectory->mkdir("MC");
  fileOut->cd("MC");
  for (float i = -2.50; i < 2.45; i+=0.05){
    myMap_MC.find(std::to_string(roundf(i*100)/100))->second.Write();
  }
  
  return;
}


