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

bool EtaFolding = true;

float AssignEtaBin(float eta, bool EtaFolding){

  float ref_eta = -9999.;

  if(eta > 0 && eta < 0.000001) eta = 0.01;
  if(eta < 0 && eta > -0.00001) eta = -0.01;

  for (float m = -2.50; m < 2.45; m+=0.05){
    if ((eta - m >= 0. && eta - m < 0.05) && !EtaFolding){
      ref_eta = roundf(m*100)/100; //Round to two decimal places
    }
    if ((fabs(eta) + m <= 0. && fabs(eta) + m > -0.05) && EtaFolding){
      ref_eta = fabs(roundf(m*100)/100);
    }
  }
  if (ref_eta == -9999.){
    cout << "eta: " << eta << "  ref_eta: " << ref_eta << endl;
    cout << "SOMETHING WENT WRONG WITH ETA BIN ASSIGNMENT" << endl;
  }
  return ref_eta;
}

float AssignPhiBin(float phi){

  float ref_phi = -9999.;

  for (float m = -M_PI; m < 3.13; m+=M_PI){ //Divide into two emispheres in phi
    //cout << "m: " << m <<endl;
    if (phi - m >= 0. && phi - m <= M_PI){
      ref_phi = roundf(m*100)/100; //Round to two decimal places
    }
  }
  if (ref_phi == -9999.){
    cout << "phi: " << phi << "  ref_phi: " << ref_phi << endl;
    cout << "SOMETHING WENT WRONG WITH PHI BIN ASSIGNMENT" << endl;
  }
  //cout << ref_phi << endl;
  return ref_phi;
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

void tree_reader_phi() {
    
  //TFile* fileIn_DATA = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/Run2017B.root");
  TFile* fileIn_DATA = new TFile("/eos/home-r/rselvati/ntuples_fBrem/2017_DATA.root");
  TFile* fileIn_MC   = new TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/DY_NLO.root");

  TFile* file_PVz_reweight_DATA = new TFile("PVz_weights_DATA.root");
  TFile* file_PVz_reweight_MC   = new TFile("PVz_weights_MC.root");
  TH1F* h_PVz_reweight = new TH1F();
  
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

  TH1F* h_PVz_DATA = new TH1F("h_PVz_DATA","",60,-15.,15.);
  TH1F* h_PVz_MC   = new TH1F("h_PVz_MC","",60,-15.,15.);

  TH2F* h2_fBrem_phi_DATA = new TH2F("h2_fBrem_phi_DATA","",100,-M_PI,M_PI,50,0.,1.);
  TH2F* h2_fBrem_phi_MC   = new TH2F("h2_fBrem_phi_MC","",100,-M_PI,M_PI,50,0.,1.);
  TH2F* h2_fBrem_pT_DATA  = new TH2F("h2_fBrem_pT_DATA","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_pT_MC    = new TH2F("h2_fBrem_pT_MC","",50,33.,80.,50,0.,1.);
  TH2F* h2_fBrem_eta_DATA = new TH2F("h2_fBrem_eta_DATA","",100,-2.5,2.5,50,0.,1.);
  TH2F* h2_fBrem_eta_MC   = new TH2F("h2_fBrem_eta_MC","",100,-2.5,2.5,50,0.,1.);

  TH1F* h_array_fBrem[2][2]         = {{h_fBrem_B_DATA,h_fBrem_B_MC},{h_fBrem_E_DATA,h_fBrem_E_MC}};
  TH1F* h_array_fBrem_Plus[2][2]    = {{h_fBrem_Plus_B_DATA,h_fBrem_Plus_B_MC},{h_fBrem_Plus_E_DATA,h_fBrem_Plus_E_MC}};
  TH1F* h_array_fBrem_Minus[2][2]   = {{h_fBrem_Minus_B_DATA,h_fBrem_Minus_B_MC},{h_fBrem_Minus_E_DATA,h_fBrem_Minus_E_MC}};
  TH1F* h_array_mass[2][2]          = {{h_mass_B_DATA,h_mass_B_MC},{h_mass_E_DATA,h_mass_E_MC}};
  TH1F* h_array_mass_negFbrem[2][2] = {{h_mass_negFbrem_B_DATA,h_mass_negFbrem_B_MC},{h_mass_negFbrem_E_DATA,h_mass_negFbrem_E_MC}};
  TH1F* h_array_el_energy[2][2]     = {{h_el_energy_DATA,h_el_energy_MC},{h_el_energy_DATA,h_el_energy_MC}}; //Same histo for barrel and endcaps
  TH1F* h_array_PVz[2][2]           = {{h_PVz_DATA,h_PVz_MC},{h_PVz_DATA,h_PVz_MC}}; //Same histo for barrel and endcaps
  TH2F* h2_array_fBrem_phi[2]       = {h2_fBrem_phi_DATA,h2_fBrem_phi_MC};
  TH2F* h2_array_fBrem_pT[2]        = {h2_fBrem_pT_DATA,h2_fBrem_pT_MC};
  TH2F* h2_array_fBrem_eta[2]       = {h2_fBrem_eta_DATA,h2_fBrem_eta_MC};

  //TH1F* h_pT_DATA  = new TH1F("h_pT_DATA","",47,33.,80.);
  //TH1I* h_nPV_DATA = new TH1I("h_nPV_DATA","",60,0.,60.);
  
  //Indexes for range in eta
  float e_range = -2.50; 
  float e_step = 0.05;
  float eplus_range = e_range + e_step;

  map<string, pair<TH1F,TH1F> > myMap_fBrem;
  map<string, pair<TH1F,TH1F> > myMap_energy;
  map<string, pair<TH1F,TH1F> > myMap_xOverX0;

  string map_index;

  //for (int d = 0; d < 100; d++){
  for (int e = 0; e < 50; e++){
    e_range = roundf(e_range*100)/100;
    eplus_range = roundf(eplus_range*100)/100;

    //Indexes fot range in phi
    float p_range = -3.14;
    float p_step = 3.14;
    float pplus_range = p_range + p_step;

    for (int p = 0; p < 2; p++){

      p_range = roundf(p_range*100)/100;
      pplus_range = roundf(pplus_range*100)/100;
      
      if(EtaFolding) map_index = to_string(fabs(e_range)) + "_" + to_string(p_range);
      else map_index = to_string(e_range) + "_" + to_string(p_range);
      //cout << "il map index: " << map_index << endl;
      TString n_fBrem_DATA = TString::Format("DATA_h_fBrem_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);
      TString n_fBrem_MC   = TString::Format("MC_h_fBrem_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);

      TString n_energy_DATA = TString::Format("DATA_h_energy_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);
      TString n_energy_MC   = TString::Format("MC_h_energy_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);

      TString n_xOverX0_DATA = TString::Format("DATA_h_xOverX0_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);
      TString n_xOverX0_MC   = TString::Format("MC_h_xOverX0_etaBin_%f", e_range) + TString::Format("_%f", eplus_range) + TString::Format("_phiBin_%f", p_range) + TString::Format("_%f", pplus_range);
       

      myMap_fBrem.insert(make_pair(map_index, make_pair( TH1F(n_fBrem_DATA,"",200,-1.,1.), TH1F(n_fBrem_MC,"",200,-1.,1.) )));

      myMap_energy.insert(make_pair(map_index, make_pair( TH1F(n_energy_DATA,"",217,33.,250.), TH1F(n_energy_MC,"",217,33.,250.) )));
      
      myMap_xOverX0.insert(make_pair(map_index, make_pair( TH1F(n_xOverX0_DATA,"",200,0.,5.), TH1F(n_xOverX0_MC,"",200,0.,5.) )));

      p_range += 3.14;
      pplus_range += 3.14;
    }

    e_range += 0.05;
    eplus_range += 0.05;
  }
  
  for (int D_or_MC = 0; D_or_MC < 2; D_or_MC++){
    tree_array[D_or_MC]->SetBranchStatus("*",0);
    tree_array[D_or_MC]->SetBranchStatus("passHltEle32DoubleEGWPTightGsf",1);
    tree_array[D_or_MC]->SetBranchStatus("passEGL1SingleEGOr",1);
    //tree_array[D_or_MC]->SetBranchStatus("passingMVA94Xwp80isoV2",1);
    tree_array[D_or_MC]->SetBranchStatus("passingCutBasedTight94XV2",1);
    //tree_array[D_or_MC]->SetBranchStatus("passingCutBasedLoose94XV2",1);
    tree_array[D_or_MC]->SetBranchStatus("el_relPfLepIso03",1);
    tree_array[D_or_MC]->SetBranchStatus("el_pt",1);
    tree_array[D_or_MC]->SetBranchStatus("el_eta",1);
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
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleFull5x5SigmaIEtaIEtaCut",1);
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleDEtaInSeedCut",1);
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleHadronicOverEMEnergyScaledCut",1);
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleRelPFIsoScaledCut",1);
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleMissingHitsCut",1);
    // tree_array[D_or_MC]->SetBranchStatus("passingCutBasedVeto94XV2GsfEleConversionVetoCut",1);
    if (D_or_MC == 1) tree_array[D_or_MC]->SetBranchStatus("totWeight",1);
    
    // int v_passVeto5x5sigma;
    // int v_passVetoDEta;
    // int v_passVetoHOverE;
    // int v_passVetoPFIso;
    // int v_passVetoMissingHits;
    // int v_passVetoConversion;

    int v_passHltEle32;
    int v_passEGL1;
    int v_passMVA;
    int v_passCutBasedID;
    int v_nPV;
    float v_PVz;
    float v_el_PFRelIso;
    float v_el_pt;
    float v_el_eta;
    float v_el_phi;
    float v_el_energy;
    float v_el_fbrem;
    float v_el_q;
    float v_el_dz;
    float v_mass;
    float v_Tag_el_pt;
    float v_Tag_sc_eta;
    float v_totWeight;
    
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleFull5x5SigmaIEtaIEtaCut",&v_passVeto5x5sigma);
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleDEtaInSeedCut",&v_passVetoDEta);
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleHadronicOverEMEnergyScaledCut",&v_passVetoHOverE);
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleRelPFIsoScaledCut",&v_passVetoPFIso);
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleMissingHitsCut",&v_passVetoMissingHits);
    // tree_array[D_or_MC]->SetBranchAddress("passingCutBasedVeto94XV2GsfEleConversionVetoCut",&v_passVetoConversion);

    tree_array[D_or_MC]->SetBranchAddress("passHltEle32DoubleEGWPTightGsf",&v_passHltEle32);
    tree_array[D_or_MC]->SetBranchAddress("passEGL1SingleEGOr",&v_passEGL1);
    //tree_array[D_or_MC]->SetBranchAddress("passingMVA94Xwp80isoV2",&v_passMVA);
    tree_array[D_or_MC]->SetBranchAddress("passingCutBasedTight94XV2",&v_passCutBasedID);
    //tree_array[D_or_MC]->SetBranchAddress("passingCutBasedLoose94XV2",&v_passCutBasedID);
    tree_array[D_or_MC]->SetBranchAddress("el_relPfLepIso03",&v_el_PFRelIso);
    tree_array[D_or_MC]->SetBranchAddress("event_nPV",&v_nPV);
    tree_array[D_or_MC]->SetBranchAddress("event_PrimaryVertex_z",&v_PVz);
    tree_array[D_or_MC]->SetBranchAddress("el_pt",&v_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("el_eta",&v_el_eta);
    tree_array[D_or_MC]->SetBranchAddress("el_phi",&v_el_phi);
    tree_array[D_or_MC]->SetBranchAddress("el_e",&v_el_energy);
    tree_array[D_or_MC]->SetBranchAddress("el_fbrem",&v_el_fbrem);
    tree_array[D_or_MC]->SetBranchAddress("el_q",&v_el_q);
    tree_array[D_or_MC]->SetBranchAddress("el_dz",&v_el_dz);
    tree_array[D_or_MC]->SetBranchAddress("tag_Ele_pt",&v_Tag_el_pt);
    tree_array[D_or_MC]->SetBranchAddress("tag_sc_eta",&v_Tag_sc_eta);
    tree_array[D_or_MC]->SetBranchAddress("mass",&v_mass);
    if (D_or_MC == 1) tree_array[D_or_MC]->SetBranchAddress("totWeight",&v_totWeight);

    if (D_or_MC == 0) h_PVz_reweight = (TH1F*)file_PVz_reweight_DATA->Get("h_reweight_PVz");
    else h_PVz_reweight = (TH1F*)file_PVz_reweight_MC->Get("h_reweight_PVz");
    
    for (Long64_t n = 0; n < tree_array[D_or_MC]->GetEntriesFast(); n++){
      tree_array[D_or_MC]->GetEntry(n);
      
      //if (!( v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passVeto5x5sigma == 1 && v_passVetoDEta == 1 && v_passVetoHOverE == 1 && v_passVetoPFIso == 1 && v_passVetoMissingHits == 1 && v_passVetoConversion == 1 && v_el_PFRelIso < 0.25 && v_el_pt > 33. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      if (!(v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passCutBasedID == 1 && v_el_PFRelIso < 0.25  && v_el_pt > 33. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      //if (!(v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passMVA == 1 && v_el_pt > 33. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      //if (!(v_passHltEle32 == 1 && v_passEGL1 == 1 && v_passMVA == 1 && v_el_pt > 28. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;
      //if (!(v_passHltEle27 == 1 && v_passCutBasedID == 1 && v_el_PFRelIso < 0.25  && v_el_pt > 28. && fabs(v_el_eta) < 2.50 && fabs(v_el_dz) < 0.01)) continue;

      if (!(v_Tag_el_pt > 35. && fabs(v_Tag_sc_eta) < 2.1)) continue;
      
      if (!(v_mass > 60. && v_mass < 120.)) continue;
      
      //if (!(v_el_fbrem > 0.)) continue;

      //if (!(v_el_pt < 43.)) continue;
      //if (!(v_nPV > 23.)) continue;
      
      //Assign the right eta bin
      float eta_bin = AssignEtaBin(v_el_eta, EtaFolding);
      if (eta_bin == -9999.) break;

      //Assign the right phi bin
      //cout << v_el_phi << endl;
      //if(v_el_phi > 3.14) v_el_phi = 3.1399;
      //if(v_el_phi < -3.14) v_el_phi = -3.1399;
      float phi_bin = AssignPhiBin(v_el_phi);

      if (phi_bin == -9999.) break;

      string map_index = to_string(eta_bin) + "_" + to_string(phi_bin);
      float EventWeight = 1.;

      //PVz_reweighting
      EventWeight = EventWeight*(h_PVz_reweight->GetBinContent(h_PVz_reweight->FindBin(v_PVz)));

      if (D_or_MC == 0){
      	myMap_fBrem.find(map_index)->second.first.Fill(v_el_fbrem,EventWeight);
      	myMap_energy.find(map_index)->second.first.Fill(v_el_energy,EventWeight);
      	myMap_xOverX0.find(map_index)->second.first.Fill(-(log(1 - v_el_fbrem)),EventWeight);
      }
      else{
      	EventWeight = v_totWeight;
	EventWeight = EventWeight*(h_PVz_reweight->GetBinContent(h_PVz_reweight->FindBin(v_PVz)));
      	myMap_fBrem.find(map_index)->second.second.Fill(v_el_fbrem,EventWeight);
      	myMap_energy.find(map_index)->second.second.Fill(v_el_energy,EventWeight);
      	myMap_xOverX0.find(map_index)->second.second.Fill(-(log(1 - v_el_fbrem)),EventWeight);
      }

      int B_or_E = -1;
      if (fabs(v_el_eta) < 1.2) B_or_E = 0; //Corresponds to BARREL
      else B_or_E = 1; //Corresponds to ENDCAPS
   
      // h_array_fBrem[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      // h_array_el_energy[B_or_E][D_or_MC]->Fill(v_el_energy,EventWeight);
      // h_array_PVz[B_or_E][D_or_MC]->Fill(v_PVz,EventWeight);
      // if (v_el_fbrem >= 0.) h_array_mass[B_or_E][D_or_MC]->Fill(v_mass,EventWeight);
      // if (v_el_fbrem < 0.) h_array_mass_negFbrem[B_or_E][D_or_MC]->Fill(v_mass,EventWeight);
      // if (v_el_q == 1) h_array_fBrem_Plus[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      // if (v_el_q == -1) h_array_fBrem_Minus[B_or_E][D_or_MC]->Fill(v_el_fbrem,EventWeight);
      // h2_array_fBrem_phi[D_or_MC]->Fill(v_el_phi,v_el_fbrem,EventWeight);      
      // h2_array_fBrem_pT[D_or_MC]->Fill(v_el_pt,v_el_fbrem,EventWeight);      
      // h2_array_fBrem_eta[D_or_MC]->Fill(v_el_eta,v_el_fbrem,EventWeight);      
    }
  }
  
  // gStyle->SetOptStat(0);
  // TFile* fileOut_plots = new TFile("histos/02_06_2021/plots_CutBased_DYNLO_custom.root","RECREATE");
  // fileOut_plots->cd();

  // TCanvas* c_fBrem_B = new TCanvas();
  // Prepare1DHisto(h_array_fBrem, 0, "f_{brem}", "Events/0.04"); // 0 stands for BARREL
  // c_fBrem_B->SaveAs("plots/02_06_2021/fBrem_B_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_E = new TCanvas();
  // Prepare1DHisto(h_array_fBrem, 1, "f_{brem}", "Events/0.04"); // 1 stands for ENDCAP
  // c_fBrem_E->SaveAs("plots/02_06_2021/fBrem_E_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_Plus_B = new TCanvas();
  // Prepare1DHisto(h_array_fBrem_Plus, 0, "f_{brem}", "Events/0.04");
  // c_fBrem_Plus_B->SaveAs("plots/02_06_2021/fBrem_Plus_B_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_Minus_B = new TCanvas();
  // Prepare1DHisto(h_array_fBrem_Minus, 0, "f_{brem}", "Events/0.04");
  // c_fBrem_Minus_B->SaveAs("plots/02_06_2021/fBrem_Minus_B_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_Plus_E = new TCanvas();
  // Prepare1DHisto(h_array_fBrem_Plus, 1, "f_{brem}", "Events/0.04");
  // c_fBrem_Plus_E->SaveAs("plots/02_06_2021/fBrem_Plus_E_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_Minus_E = new TCanvas();
  // Prepare1DHisto(h_array_fBrem_Minus, 1, "f_{brem}", "Events/0.04");
  // c_fBrem_Minus_E->SaveAs("plots/02_06_2021/fBrem_Minus_E_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_mass_B = new TCanvas();
  // Prepare1DHisto(h_array_mass, 0, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  // c_mass_B->SaveAs("plots/02_06_2021/mass_B_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_mass_E = new TCanvas();
  // Prepare1DHisto(h_array_mass, 1, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  // c_mass_E->SaveAs("plots/02_06_2021/mass_E_CutBased_DYNLO_custom.pdf");

  // TCanvas* c_mass_negFbrem_B = new TCanvas();
  // Prepare1DHisto(h_array_mass_negFbrem, 0, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  // c_mass_negFbrem_B->SaveAs("plots/02_06_2021/mass_negFbrem_B_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_mass_negFbrem_E = new TCanvas();
  // Prepare1DHisto(h_array_mass_negFbrem, 1, "#it{m}_{ee} (GeV)", "Events/2.0 GeV");
  // c_mass_negFbrem_E->SaveAs("plots/02_06_2021/mass_negFbrem_E_CutBased_DYNLO_custom.pdf");

  // TCanvas* c_el_energy = new TCanvas();
  // Prepare1DHisto(h_array_el_energy, 0, "E^{e} (GeV)", "Events/1.0 GeV");
  // c_el_energy->SaveAs("plots/02_06_2021/el_energy_CutBased_DYNLO_custom.pdf");

  // TCanvas* c_PVz = new TCanvas();
  // Prepare1DHisto(h_array_PVz, 0, "PV z (cm)", "Events/0.5 cm");
  // c_PVz->SaveAs("plots/02_06_2021/PVz_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_phi_DATA = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_phi, 0, "#phi^{e}", "f_{brem}"); // 0 stands for DATA
  // c_fBrem_vs_phi_DATA->SaveAs("plots/02_06_2021/fBrem_vs_phi_DATA_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_phi_MC = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_phi, 1, "#phi^{e}", "f_{brem}"); // 1 stands for MC
  // c_fBrem_vs_phi_MC->SaveAs("plots/02_06_2021/fBrem_vs_phi_MC_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_pT_DATA = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_pT, 0, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  // c_fBrem_vs_pT_DATA->SaveAs("plots/02_06_2021/fBrem_vs_pT_DATA_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_pT_MC = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_pT, 1, "#it{p}_{T}^{e} (GeV)", "f_{brem}");
  // c_fBrem_vs_pT_MC->SaveAs("plots/02_06_2021/fBrem_vs_pT_MC_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_eta_DATA = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_eta, 0, "#eta^{e}", "f_{brem}");
  // c_fBrem_vs_eta_DATA->SaveAs("plots/02_06_2021/fBrem_vs_eta_DATA_CutBased_DYNLO_custom.pdf");
  
  // TCanvas* c_fBrem_vs_eta_MC = new TCanvas();
  // Prepare2DHisto(h2_array_fBrem_eta, 1, "#eta^{e}", "f_{brem}");
  // c_fBrem_vs_eta_MC->SaveAs("plots/02_06_2021/fBrem_vs_eta_MC_CutBased_DYNLO_custom.pdf");

  //h_pT_DATA->Write();
  //h_nPV_DATA->Write();
  
  TFile* fileOut = new TFile("histos/01_07_2021/fBrem_etaFolded_EtaPhiBins_PVz_reweighted.root","RECREATE");
  fileOut->cd();
  gDirectory->mkdir("DATA");
  fileOut->cd("DATA");
  for (float e = 0.05; e < 2.50; e+=0.05){
    for (float p = -3.14; p < 3.13; p+=3.14){
      string map_index = to_string(roundf(e*100)/100) + "_" + to_string(roundf(p*100)/100);
      //cout << "map_index_DATA:" << map_index;
      
      myMap_fBrem.find(map_index)->second.first.Write();
      myMap_energy.find(map_index)->second.first.Write();
      myMap_xOverX0.find(map_index)->second.first.Write();
    }
  }
  fileOut->cd();
  gDirectory->mkdir("MC");
  fileOut->cd("MC");
  for (float e = 0.05; e < 2.50; e+=0.05){
    for (float p = -3.14; p < 3.13; p+=3.14){
      string map_index = to_string(roundf(e*100)/100) + "_" + to_string(roundf(p*100)/100);
      //cout << "map_index_DATA:" << map_index;
      myMap_fBrem.find(map_index)->second.second.Write();
      myMap_energy.find(map_index)->second.second.Write();
      myMap_xOverX0.find(map_index)->second.second.Write();
    }
  }
  
  return;
}


