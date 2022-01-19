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

using namespace std;

float AssignEtaBin(float eta){

  float ref_eta = -9999.;

  for (float m = -2.50; m < 2.45; m+=0.05){
    if (eta - m >= 0. && eta - m < 0.05){
      ref_eta = roundf(m*100)/100; //Round to two decimal places
    }
  }
  if (ref_eta == -9999.){
    cout << "eta: " << eta << "  ref_eta: " << ref_eta << endl;
    cout << "SOMETHING WENT WRONG WITH ETA BIN ASSIGNMENT" << endl;
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

void tree_reader_2016() {
    
  // TFile* fileIn_DATA_0 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016B_L1matched.root");
  // TFile* fileIn_DATA_1 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016C_L1matched.root");
  // TFile* fileIn_DATA_2 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016D_L1matched.root");
  // TFile* fileIn_DATA_3 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016E_L1matched.root");
  // TFile* fileIn_DATA_4 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016F_L1matched.root");
  TFile* fileIn_DATA_5 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged/Run2016F_L1merged.root");
  TFile* fileIn_DATA_6 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged/Run2016G_L1matched.root");
  TFile* fileIn_DATA_7 = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged/Run2016H_L1matched.root");
  
  // TFile* fileIn_MC_0   = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/DY_NLO_L1matched.root");
  TFile* fileIn_MC_1   = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged/DY_NLO_L1matched.root");

  TFile* file_PVz_reweight_DATA = new TFile("PVz_weights_DATA_2016postVFP.root");
  TFile* file_PVz_reweight_MC   = new TFile("PVz_weights_MC_2016postVFP.root");
  TH1F* h_PVz_reweight = new TH1F();

  // TTree* treeIn_DATA_0 = (TTree*)fileIn_DATA_0->Get("tnpEleTrig/fitter_tree");
  // TTree* treeIn_DATA_1 = (TTree*)fileIn_DATA_1->Get("tnpEleTrig/fitter_tree");
  // TTree* treeIn_DATA_2 = (TTree*)fileIn_DATA_2->Get("tnpEleTrig/fitter_tree");
  // TTree* treeIn_DATA_3 = (TTree*)fileIn_DATA_3->Get("tnpEleTrig/fitter_tree");
  // TTree* treeIn_DATA_4 = (TTree*)fileIn_DATA_4->Get("tnpEleTrig/fitter_tree");
  TTree* treeIn_DATA_5 = (TTree*)fileIn_DATA_5->Get("tnpEleTrig/fitter_tree");
  TTree* treeIn_DATA_6 = (TTree*)fileIn_DATA_6->Get("tnpEleTrig/fitter_tree");
  TTree* treeIn_DATA_7 = (TTree*)fileIn_DATA_7->Get("tnpEleTrig/fitter_tree");

  // TTree* treeIn_MC_0   = (TTree*)fileIn_MC_0->Get("tnpEleTrig/fitter_tree");
  TTree* treeIn_MC_1   = (TTree*)fileIn_MC_1->Get("tnpEleTrig/fitter_tree");
  
  //TTree* tree_array[10] = {treeIn_DATA_0,treeIn_DATA_1,treeIn_DATA_2,treeIn_DATA_3,treeIn_DATA_4,treeIn_DATA_5,treeIn_DATA_6,treeIn_DATA_7,treeIn_MC_0,treeIn_MC_1};
  //TTree* tree_array[6] = {treeIn_DATA_0,treeIn_DATA_1,treeIn_DATA_2,treeIn_DATA_3,treeIn_DATA_4,treeIn_MC_0}; //Pre-VFP
  TTree* tree_array[4] = {treeIn_DATA_5,treeIn_DATA_6,treeIn_DATA_7,treeIn_MC_1}; //Post-VFP
  std::cout << std::size(tree_array) << std::endl; 

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
  
  TH1F* h_mass_B_DATA = new TH1F("h_mass_B_DATA","",120,60.,120.);
  TH1F* h_mass_E_DATA = new TH1F("h_mass_E_DATA","",120,60.,120.);
  TH1F* h_mass_B_MC   = new TH1F("h_mass_B_MC","",120,60.,120.);
  TH1F* h_mass_E_MC   = new TH1F("h_mass_E_MC","",120,60.,120.);

  TH1F* h_mass_negFbrem_B_DATA = new TH1F("h_mass_negFbrem_B_DATA","",120,60.,120.);
  TH1F* h_mass_negFbrem_E_DATA = new TH1F("h_mass_negFbrem_E_DATA","",120,60.,120.);
  TH1F* h_mass_negFbrem_B_MC   = new TH1F("h_mass_negFbrem_B_MC","",120,60.,120.);
  TH1F* h_mass_negFbrem_E_MC   = new TH1F("h_mass_negFbrem_E_MC","",120,60.,120.);

  TH1F* h_el_energy_DATA = new TH1F("h_el_energy_DATA","",270,30.,300.);
  TH1F* h_el_energy_MC   = new TH1F("h_el_energy_MC","",270,30.,300.);

  TH1F* h_el_pT_DATA = new TH1F("h_el_pT_DATA","",117,33.,150.);
  TH1F* h_el_pT_MC   = new TH1F("h_el_pT_MC","",117,33.,150.);

  TH1F* h_PVz_DATA = new TH1F("h_PVz_DATA","",60,-15.,15.);
  TH1F* h_PVz_MC   = new TH1F("h_PVz_MC","",60,-15.,15.);

  TH2F* h2_fBrem_phi_DATA   = new TH2F("h2_fBrem_phi_DATA","",100,-3.14,3.14,50,0.,1.);
  TH2F* h2_fBrem_phi_MC     = new TH2F("h2_fBrem_phi_MC","",100,-3.14,3.14,50,0.,1.);
  TH2F* h2_fBrem_pT_DATA    = new TH2F("h2_fBrem_pT_DATA","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_pT_MC      = new TH2F("h2_fBrem_pT_MC","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_eta_DATA   = new TH2F("h2_fBrem_eta_DATA","",100,-2.5,2.5,50,0.,1.);
  TH2F* h2_fBrem_eta_MC     = new TH2F("h2_fBrem_eta_MC","",100,-2.5,2.5,50,0.,1.);
  TH2F* h2_etaDiff_PVz_DATA = new TH2F("h2_etaDiff_PVz_DATA","",30,-15.,15.,40,-0.2,0.2);
  TH2F* h2_etaDiff_PVz_MC   = new TH2F("h2_etaDiff_PVz_MC","",30,-15.,15.,40,-0.2,0.2);

  TH1F* h_array_fBrem[2][2]         = {{h_fBrem_B_DATA,h_fBrem_B_MC},{h_fBrem_E_DATA,h_fBrem_E_MC}};
  TH1F* h_array_fBrem_Plus[2][2]    = {{h_fBrem_Plus_B_DATA,h_fBrem_Plus_B_MC},{h_fBrem_Plus_E_DATA,h_fBrem_Plus_E_MC}};
  TH1F* h_array_fBrem_Minus[2][2]   = {{h_fBrem_Minus_B_DATA,h_fBrem_Minus_B_MC},{h_fBrem_Minus_E_DATA,h_fBrem_Minus_E_MC}};
  TH1F* h_array_mass[2][2]          = {{h_mass_B_DATA,h_mass_B_MC},{h_mass_E_DATA,h_mass_E_MC}};
  TH1F* h_array_mass_negFbrem[2][2] = {{h_mass_negFbrem_B_DATA,h_mass_negFbrem_B_MC},{h_mass_negFbrem_E_DATA,h_mass_negFbrem_E_MC}};
  TH1F* h_array_el_energy[2][2]     = {{h_el_energy_DATA,h_el_energy_MC},{h_el_energy_DATA,h_el_energy_MC}}; //Same histo for barrel and endcaps
  TH1F* h_array_el_pT[2][2]         = {{h_el_pT_DATA,h_el_pT_MC},{h_el_pT_DATA,h_el_pT_MC}}; //Same histo for barrel and endcaps
  TH1F* h_array_PVz[2][2]           = {{h_PVz_DATA,h_PVz_MC},{h_PVz_DATA,h_PVz_MC}}; //Same histo for barrel and endcaps
  TH2F* h2_array_fBrem_phi[2]       = {h2_fBrem_phi_DATA,h2_fBrem_phi_MC};
  TH2F* h2_array_fBrem_pT[2]        = {h2_fBrem_pT_DATA,h2_fBrem_pT_MC};
  TH2F* h2_array_fBrem_eta[2]       = {h2_fBrem_eta_DATA,h2_fBrem_eta_MC};
  TH2F* h2_array_etaDiff_PVz[2]     = {h2_etaDiff_PVz_DATA,h2_etaDiff_PVz_MC};

  //TH1F* h_pT_DATA  = new TH1F("h_pT_DATA","",47,33.,80.);
  //TH1I* h_nPV_DATA = new TH1I("h_nPV_DATA","",60,0.,60.);
  
  //Indexes for range in eta
  float e_range = -2.50; 
  float e_step = 0.05;
  float eplus_range = e_range + e_step;

  int local_D_or_MC = 999;

  map<string, pair<TH1F,TH1F> > myMap_fBrem;
  map<string, pair<TH1F,TH1F> > myMap_energy;
  map<string, pair<TH1F,TH1F> > myMap_xOverX0;
  map<string, pair<TH2F,TH2F> > myMap_fBrem_vs_PVz;

  for (int e = 0; e < 100; e++){
    e_range = roundf(e_range*100)/100;
    eplus_range = roundf(eplus_range*100)/100;

    string map_index = to_string(e_range);
    TString n_fBrem_DATA = TString::Format("DATA_h_fBrem_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);
    TString n_fBrem_MC   = TString::Format("MC_h_fBrem_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);

    TString n_energy_DATA = TString::Format("DATA_h_energy_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);
    TString n_energy_MC   = TString::Format("MC_h_energy_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);

    TString n_xOverX0_DATA = TString::Format("DATA_h_xOverX0_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);
    TString n_xOverX0_MC   = TString::Format("MC_h_xOverX0_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);

    TString n_fBrem_vs_PVz_DATA = TString::Format("DATA_h2_fBrem_vs_PVz_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);
    TString n_fBrem_vs_PVz_MC   = TString::Format("MC_h2_fBrem_vs_PVz_etaBins_%f_", e_range) + TString::Format("%f", eplus_range);

    myMap_fBrem.insert(make_pair(map_index, make_pair( TH1F(n_fBrem_DATA,"",200,-1.,1.), TH1F(n_fBrem_MC,"",200,-1.,1.) )));
    
    myMap_energy.insert(make_pair(map_index, make_pair( TH1F(n_energy_DATA,"",217,33.,250.), TH1F(n_energy_MC,"",217,33.,250.) )));
    
    myMap_xOverX0.insert(make_pair(map_index, make_pair( TH1F(n_xOverX0_DATA,"",200,0.,5.), TH1F(n_xOverX0_MC,"",200,0.,5.) )));

    myMap_fBrem_vs_PVz.insert(make_pair(map_index, make_pair( TH2F(n_fBrem_vs_PVz_DATA,"",30,-15.,15.,50,-1.,1.), TH2F(n_fBrem_vs_PVz_MC,"",30,-15.,15.,50,-1.,1.) )));

    e_range += 0.05;
    eplus_range += 0.05;
  }
  
  for (int D_or_MC = 0; D_or_MC < std::size(tree_array); D_or_MC++){
    tree_array[D_or_MC]->SetBranchStatus("*",0);
    tree_array[D_or_MC]->SetBranchStatus("passHltEle27WPTightGsf",1);
    tree_array[D_or_MC]->SetBranchStatus("passingCutBasedTight94XV2",1);
    tree_array[D_or_MC]->SetBranchStatus("el_relPfLepIso03",1);
    tree_array[D_or_MC]->SetBranchStatus("el_pt",1);
    tree_array[D_or_MC]->SetBranchStatus("el_eta",1);
    tree_array[D_or_MC]->SetBranchStatus("el_sc_eta",1);
    tree_array[D_or_MC]->SetBranchStatus("el_phi",1);
    tree_array[D_or_MC]->SetBranchStatus("el_e",1);
    tree_array[D_or_MC]->SetBranchStatus("tag_Ele_pt",1);
    tree_array[D_or_MC]->SetBranchStatus("tag_sc_eta",1);
    tree_array[D_or_MC]->SetBranchStatus("mass",1);
    tree_array[D_or_MC]->SetBranchStatus("el_fbrem",1);
    tree_array[D_or_MC]->SetBranchStatus("el_q",1);
    tree_array[D_or_MC]->SetBranchStatus("el_dz",1);
    tree_array[D_or_MC]->SetBranchStatus("event_nPV",1);
    tree_array[D_or_MC]->SetBranchStatus("event_PrimaryVertex_z",1);
    //Turn on the weight branch for MC
    //if (D_or_MC > 7) tree_array[D_or_MC]->SetBranchStatus("totWeight",1);
    //if (D_or_MC > 4) tree_array[D_or_MC]->SetBranchStatus("totWeight",1); //PreVFP
    if (D_or_MC > 2) tree_array[D_or_MC]->SetBranchStatus("totWeight",1); //PostVFP
   
    int v_passHltEle27;
    int v_passCutBasedID;
    int v_nPV;
    float v_PVz;
    float v_el_PFRelIso;
    float v_el_pt;
    float v_el_eta;
    float v_el_sc_eta;
    float v_el_phi;
    float v_el_energy;
    float v_el_fbrem;
    float v_el_q;
    float v_el_dz;
    float v_mass;
    float v_Tag_el_pt;
    float v_Tag_sc_eta;
    float v_totWeight;
    

    tree_array[D_or_MC]->SetBranchAddress("passHltEle27WPTightGsf",&v_passHltEle27);
    tree_array[D_or_MC]->SetBranchAddress("passingCutBasedTight94XV2",&v_passCutBasedID);
    tree_array[D_or_MC]->SetBranchAddress("el_relPfLepIso03",&v_el_PFRelIso);
    tree_array[D_or_MC]->SetBranchAddress("event_nPV",&v_nPV);
    tree_array[D_or_MC]->SetBranchAddress("event_PrimaryVertex_z",&v_PVz);
    tree_array[D_or_MC]->SetBranchAddress("el_pt",&v_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("el_eta",&v_el_eta);
    tree_array[D_or_MC]->SetBranchAddress("el_sc_eta",&v_el_sc_eta);
    tree_array[D_or_MC]->SetBranchAddress("el_phi",&v_el_phi);
    tree_array[D_or_MC]->SetBranchAddress("el_e",&v_el_energy);
    tree_array[D_or_MC]->SetBranchAddress("el_fbrem",&v_el_fbrem);
    tree_array[D_or_MC]->SetBranchAddress("el_q",&v_el_q);
    tree_array[D_or_MC]->SetBranchAddress("el_dz",&v_el_dz);
    tree_array[D_or_MC]->SetBranchAddress("tag_Ele_pt",&v_Tag_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("tag_sc_eta",&v_Tag_sc_eta);
    tree_array[D_or_MC]->SetBranchAddress("mass",&v_mass);
    //Set the address of the weight branch for MC
    //if (D_or_MC > 7) tree_array[D_or_MC]->SetBranchAddress("totWeight",&v_totWeight);
    //if (D_or_MC > 4) tree_array[D_or_MC]->SetBranchAddress("totWeight",&v_totWeight); //PreVFP
    if (D_or_MC > 2) tree_array[D_or_MC]->SetBranchAddress("totWeight",&v_totWeight); //PostVFP

    //if (D_or_MC <= 7) h_PVz_reweight = (TH1F*)file_PVz_reweight_DATA->Get("h_reweight_PVz");
    //if (D_or_MC <= 4) h_PVz_reweight = (TH1F*)file_PVz_reweight_DATA->Get("h_reweight_PVz"); //PreVFT
    if (D_or_MC <= 2) h_PVz_reweight = (TH1F*)file_PVz_reweight_DATA->Get("h_reweight_PVz"); //PostVFT
    else h_PVz_reweight = (TH1F*)file_PVz_reweight_MC->Get("h_reweight_PVz");

    std::cout << "loop #" << D_or_MC << std::endl;
    
    for (Long64_t n = 0; n < tree_array[D_or_MC]->GetEntriesFast(); n++){
      tree_array[D_or_MC]->GetEntry(n);
      
      if (!(v_passHltEle27 == 1 && v_passCutBasedID == 1 && v_el_PFRelIso < 0.25  && v_el_pt > 33. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;

      if (!(v_Tag_el_pt > 35. && fabs(v_Tag_sc_eta) < 2.1)) continue;
      
      if (!(v_mass > 60. && v_mass < 120.)) continue;
      
      //if (!(v_el_fbrem > 0.)) continue;
      //if(!(v_PVz > 0.)) continue;
      //if (!(v_el_pt < 43.)) continue;
      //if (!(v_nPV > 23.)) continue;
      
      //Assign the right eta bin
      float eta_bin = AssignEtaBin(v_el_eta);
      if (eta_bin == -9999.) break;
     
      //Fill the maps
      string map_index = to_string(eta_bin);
      float EventWeight = 1.;
      
      //PVz_reweighting
      EventWeight = EventWeight*(h_PVz_reweight->GetBinContent(h_PVz_reweight->FindBin(v_PVz)));

      //if (D_or_MC <= 7){
      //if (D_or_MC <= 4){ //PreVFP
      if (D_or_MC <= 2){ //PostVFP
	myMap_fBrem.find(map_index)->second.first.Fill(v_el_fbrem,EventWeight);
	myMap_energy.find(map_index)->second.first.Fill(v_el_energy,EventWeight);
	myMap_xOverX0.find(map_index)->second.first.Fill(-(log(1 - v_el_fbrem)),EventWeight);
	myMap_fBrem_vs_PVz.find(map_index)->second.first.Fill(v_PVz,v_el_fbrem,EventWeight);
      }
      else{
	EventWeight *= v_totWeight;
	myMap_fBrem.find(map_index)->second.second.Fill(v_el_fbrem,EventWeight);
	myMap_energy.find(map_index)->second.second.Fill(v_el_energy,EventWeight);
	myMap_xOverX0.find(map_index)->second.second.Fill(-(log(1 - v_el_fbrem)),EventWeight);
	myMap_fBrem_vs_PVz.find(map_index)->second.second.Fill(v_PVz,v_el_fbrem,EventWeight);
      }

      int B_or_E = -1;
      if (fabs(v_el_eta) < 1.2) B_or_E = 0; //Corresponds to BARREL
      else B_or_E = 1; //Corresponds to ENDCAPS
      
      //Fill the histograms stored in arrays
      //if (D_or_MC <= 7) local_D_or_MC = 0;
      //if (D_or_MC <= 4) local_D_or_MC = 0; //PreVFP
      if (D_or_MC <= 2) local_D_or_MC = 0; //PostVFP
      else local_D_or_MC = 1;
      h_array_fBrem[B_or_E][local_D_or_MC]->Fill(v_el_fbrem,EventWeight);
      h_array_el_energy[B_or_E][local_D_or_MC]->Fill(v_el_energy,EventWeight);
      h_array_el_pT[B_or_E][local_D_or_MC]->Fill(v_el_pt,EventWeight);
      h_array_PVz[B_or_E][local_D_or_MC]->Fill(v_PVz,EventWeight);
      if (v_el_fbrem >= 0.) h_array_mass[B_or_E][local_D_or_MC]->Fill(v_mass,EventWeight);
      if (v_el_fbrem < 0.) h_array_mass_negFbrem[B_or_E][local_D_or_MC]->Fill(v_mass,EventWeight);
      if (v_el_q == 1) h_array_fBrem_Plus[B_or_E][local_D_or_MC]->Fill(v_el_fbrem,EventWeight);
      if (v_el_q == -1) h_array_fBrem_Minus[B_or_E][local_D_or_MC]->Fill(v_el_fbrem,EventWeight);
      h2_array_fBrem_phi[local_D_or_MC]->Fill(v_el_phi,v_el_fbrem,EventWeight);      
      h2_array_fBrem_pT[local_D_or_MC]->Fill(v_el_pt,v_el_fbrem,EventWeight);      
      h2_array_fBrem_eta[local_D_or_MC]->Fill(v_el_eta,v_el_fbrem,EventWeight);      
      h2_array_etaDiff_PVz[local_D_or_MC]->Fill(v_PVz,(v_el_eta-v_el_sc_eta),EventWeight);      
    }
  }
  
  // gStyle->SetOptStat(0);
  TFile* fileOut_plots = new TFile("histos/09_12_2021/plots_CutBased_DYNLO_2016postVFP_PVzReweighted.root","RECREATE");
  fileOut_plots->cd();

  TCanvas* c_fBrem_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem, 0, "f_{brem}", "Events/0.04"); // 0 stands for BARREL
  c_fBrem_B->SaveAs("plots/09_12_2021/fBrem_B_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem, 1, "f_{brem}", "Events/0.04"); // 1 stands for ENDCAP
  c_fBrem_E->SaveAs("plots/09_12_2021/fBrem_E_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_Plus_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Plus, 0, "f_{brem}", "Events/0.04");
  c_fBrem_Plus_B->SaveAs("plots/09_12_2021/fBrem_Plus_B_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_Minus_B = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Minus, 0, "f_{brem}", "Events/0.04");
  c_fBrem_Minus_B->SaveAs("plots/09_12_2021/fBrem_Minus_B_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_Plus_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Plus, 1, "f_{brem}", "Events/0.04");
  c_fBrem_Plus_E->SaveAs("plots/09_12_2021/fBrem_Plus_E_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_Minus_E = new TCanvas();
  Prepare1DHisto(h_array_fBrem_Minus, 1, "f_{brem}", "Events/0.04");
  c_fBrem_Minus_E->SaveAs("plots/09_12_2021/fBrem_Minus_E_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_mass_B = new TCanvas();
  Prepare1DHisto(h_array_mass, 0, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_B->SaveAs("plots/09_12_2021/mass_B_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_mass_E = new TCanvas();
  Prepare1DHisto(h_array_mass, 1, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_E->SaveAs("plots/09_12_2021/mass_E_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_mass_negFbrem_B = new TCanvas();
  Prepare1DHisto(h_array_mass_negFbrem, 0, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_negFbrem_B->SaveAs("plots/09_12_2021/mass_negFbrem_B_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_mass_negFbrem_E = new TCanvas();
  Prepare1DHisto(h_array_mass_negFbrem, 1, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  c_mass_negFbrem_E->SaveAs("plots/09_12_2021/mass_negFbrem_E_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_el_energy = new TCanvas();
  Prepare1DHisto(h_array_el_energy, 0, "E^{e} (GeV)", "Events/1.0 GeV");
  c_el_energy->SaveAs("plots/09_12_2021/el_energy_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_el_pT = new TCanvas();
  Prepare1DHisto(h_array_el_pT, 0, "#it{p}_{T}^{e} (GeV)", "Events/1.0 GeV");
  c_el_pT->SaveAs("plots/09_12_2021/el_pT_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_PVz = new TCanvas();
  Prepare1DHisto(h_array_PVz, 0, "PV z (cm)", "Events/0.5 cm");
  c_PVz->SaveAs("plots/09_12_2021/PVz_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_phi_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_phi, 0, "#phi^{e}", "f_{brem}"); // 0 stands for DATA
  c_fBrem_vs_phi_DATA->SaveAs("plots/09_12_2021/fBrem_vs_phi_DATA_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_phi_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_phi, 1, "#phi^{e}", "f_{brem}"); // 1 stands for MC
  c_fBrem_vs_phi_MC->SaveAs("plots/09_12_2021/fBrem_vs_phi_MC_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_pT_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_pT, 0, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  c_fBrem_vs_pT_DATA->SaveAs("plots/09_12_2021/fBrem_vs_pT_DATA_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_pT_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_pT, 1, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  c_fBrem_vs_pT_MC->SaveAs("plots/09_12_2021/fBrem_vs_pT_MC_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_eta_DATA = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_eta, 0, "#eta^{e}", "f_{brem}");
  c_fBrem_vs_eta_DATA->SaveAs("plots/09_12_2021/fBrem_vs_eta_DATA_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");
  
  TCanvas* c_fBrem_vs_eta_MC = new TCanvas();
  Prepare2DHisto(h2_array_fBrem_eta, 1, "#eta^{e}", "f_{brem}");
  c_fBrem_vs_eta_MC->SaveAs("plots/09_12_2021/fBrem_vs_eta_MC_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_etaDiff_vs_PVz_DATA = new TCanvas();
  Prepare2DHisto(h2_array_etaDiff_PVz, 0, "PV_{z}","#eta^{e}-#eta^{e}_{SC}");
  c_etaDiff_vs_PVz_DATA->SaveAs("plots/09_12_2021/etaDiff_vs_PVz_DATA_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  TCanvas* c_etaDiff_vs_PVz_MC = new TCanvas();
  Prepare2DHisto(h2_array_etaDiff_PVz, 1, "PV_{z}","#eta^{e}-#eta^{e}_{SC}");
  c_etaDiff_vs_PVz_MC->SaveAs("plots/09_12_2021/etaDiff_vs_PVz_MC_CutBased_DYNLO_2016postVFP_PVzReweighted.pdf");

  //h_pT_DATA->Write();
  //h_nPV_DATA->Write();
  
  TFile* fileOut = new TFile("histos/09_12_2021/fBrem_notFolded_2016postVFP_PVzReweighted.root","RECREATE");
  fileOut->cd();
  gDirectory->mkdir("DATA");
  fileOut->cd("DATA");
  for (float i = -2.50; i < 2.45; i+=0.05){
    myMap_fBrem.find(to_string(roundf(i*100)/100))->second.first.Write();
    myMap_energy.find(to_string(roundf(i*100)/100))->second.first.Write();
    myMap_xOverX0.find(to_string(roundf(i*100)/100))->second.first.Write();
    myMap_fBrem_vs_PVz.find(to_string(roundf(i*100)/100))->second.first.Write();
  }
  fileOut->cd();
  gDirectory->mkdir("MC");
  fileOut->cd("MC");
  for (float i = -2.50; i < 2.45; i+=0.05){
    myMap_fBrem.find(to_string(roundf(i*100)/100))->second.second.Write();
    myMap_energy.find(to_string(roundf(i*100)/100))->second.second.Write();
    myMap_xOverX0.find(to_string(roundf(i*100)/100))->second.second.Write();
    myMap_fBrem_vs_PVz.find(to_string(roundf(i*100)/100))->second.second.Write();
  }
  
  return;
  }


