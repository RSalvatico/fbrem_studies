from ROOT import TFile, TH1D, TH2D, TCanvas, TGraphErrors, gDirectory, gStyle, kOrange
import math
import numpy as np
import copy
from array import array

#fileIn_DATA = TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/Run2017B.root")
fileIn_DATA = TFile("../ntuples/2017_DATA.root")
#fileIn_DATA = TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-06-09/2016/merged/Run2016G.root")
fileIn_MC   = TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/DY_LO.root")
#fileIn_MC   = TFile("/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-06-09/2016/merged/DY_LO.root")
treeIn_DATA = fileIn_DATA.Get("tnpEleTrig/fitter_tree")
treeIn_MC   = fileIn_MC.Get("tnpEleTrig/fitter_tree")

tree_list = [treeIn_DATA,treeIn_MC]

h_fBrem_B_DATA = TH1D("h_fBrem_B_DATA","",25,0.,1.)
h_fBrem_E_DATA = TH1D("h_fBrem_E_DATA","",25,0.,1.)
h_fBrem_B_MC   = TH1D("h_fBrem_B_MC","",25,0.,1.)
h_fBrem_E_MC   = TH1D("h_fBrem_E_MC","",25,0.,1.)

h_fBrem_Plus_B_DATA = TH1D("h_fBrem_Plus_B_DATA","",25,0.,1.)
h_fBrem_Plus_E_DATA = TH1D("h_fBrem_Plus_E_DATA","",25,0.,1.)
h_fBrem_Plus_B_MC   = TH1D("h_fBrem_Plus_B_MC","",25,0.,1.)
h_fBrem_Plus_E_MC   = TH1D("h_fBrem_Plus_E_MC","",25,0.,1.)

h_fBrem_Minus_B_DATA = TH1D("h_fBrem_Minus_B_DATA","",25,0.,1.)
h_fBrem_Minus_E_DATA = TH1D("h_fBrem_Minus_E_DATA","",25,0.,1.)
h_fBrem_Minus_B_MC   = TH1D("h_fBrem_Minus_B_MC","",25,0.,1.)
h_fBrem_Minus_E_MC   = TH1D("h_fBrem_Minus_E_MC","",25,0.,1.)

h_mass_B_DATA = TH1D("h_mass_B_DATA","",30,60.,120.)
h_mass_E_DATA = TH1D("h_mass_E_DATA","",30,60.,120.)
h_mass_B_MC   = TH1D("h_mass_B_MC","",30,60.,120.)
h_mass_E_MC   = TH1D("h_mass_E_MC","",30,60.,120.)

h2_fBrem_phi_DATA = TH2D("h2_fBrem_phi_DATA","",100,-3.14,3.14,50,0.,1.)
h2_fBrem_phi_MC   = TH2D("h2_fBrem_phi_MC","",100,-3.14,3.14,50,0.,1.)
h2_fBrem_pT_DATA  = TH2D("h2_fBrem_pT_DATA","",50,33.,80.,50,0.,1.)
h2_fBrem_pT_MC    = TH2D("h2_fBrem_pT_MC","",50,33.,80.,50,0.,1.)
h2_fBrem_eta_DATA = TH2D("h2_fBrem_eta_DATA","",100,-2.5,2.5,50,0.,1.)
h2_fBrem_eta_MC   = TH2D("h2_fBrem_eta_MC","",100,-2.5,2.5,50,0.,1.)

h_dict_DATA = dict()
h_dict_MC   = dict()

for i in np.arange(-2.50,2.50,0.05):
    str_i = '%.2f' % i #Generate a string with two decimal places
    str_iplus = '%.2f' %(i+0.05)
    str_i = str_i.replace(".","_")
    str_iplus = str_iplus.replace(".","_")
    #h_dict_DATA[i] = TH1D("DATA_h_xOverX0_etaBins_" + str_i + "_" + str_iplus,"",50,0.,4.) 
    #h_dict_MC[i]   = TH1D("MC_h_xOverX0_etaBins_" + str_i + "_" + str_iplus,"",50,0.,4.) 

    h_dict_DATA[i] = TH1D("DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus,"",50,0.,1.) 
    h_dict_MC[i]   = TH1D("MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus,"",50,0.,1.) 


for tree in tree_list:
    #ciao = 0
    #Only activate desired branches to speed up
    tree.SetBranchStatus("*",0)
    tree.SetBranchStatus("passHltEle32DoubleEGWPTightGsf",1)
    tree.SetBranchStatus("passEGL1SingleEGOr",1)
    tree.SetBranchStatus("passingMVA94Xwp80isoV2",1)
    tree.SetBranchStatus("el_pt",1)
    tree.SetBranchStatus("el_eta",1)
    tree.SetBranchStatus("el_phi",1)
    tree.SetBranchStatus("tag_Ele_pt",1)
    tree.SetBranchStatus("tag_sc_eta",1)
    tree.SetBranchStatus("mass",1)
    tree.SetBranchStatus("el_fbrem",1)
    tree.SetBranchStatus("el_q",1)
    if tree == treeIn_MC:
        tree.SetBranchStatus("totWeight",1)

    v_passHltEle32 = array('i',[0])
    v_passEGL1     = array('i',[0])
    v_passMVA      = array('i',[0])
    v_el_pt        = array('f',[0])
    v_el_eta       = array('f',[0])
    v_el_phi       = array('f',[0])
    v_Tag_el_pt    = array('f',[0])
    v_Tag_sc_eta   = array('f',[0])
    v_mass         = array('f',[0])
    v_el_fbrem     = array('f',[0])
    v_el_q         = array('f',[0])
    v_totWeight    = array('f',[0])
    
    tree.SetBranchAddress("passHltEle32DoubleEGWPTightGsf",v_passHltEle32)
    tree.SetBranchAddress("passEGL1SingleEGOr",v_passEGL1)
    tree.SetBranchAddress("passingMVA94Xwp80isoV2",v_passMVA)
    tree.SetBranchAddress("el_pt",v_el_pt)
    tree.SetBranchAddress("el_eta",v_el_eta)
    tree.SetBranchAddress("el_phi",v_el_phi)
    tree.SetBranchAddress("el_fbrem",v_el_fbrem)
    tree.SetBranchAddress("el_q",v_el_q)
    tree.SetBranchAddress("tag_Ele_pt",v_Tag_el_pt)
    tree.SetBranchAddress("tag_sc_eta",v_Tag_sc_eta)
    tree.SetBranchAddress("mass",v_mass)
    if tree == treeIn_MC:
        tree.SetBranchAddress("totWeight",v_totWeight)

    #i = 0
    for entry in xrange(tree.GetEntriesFast()):
    #while tree.GetEntry(i):
        tree.GetEntry(entry)
        #i += 1
        #ciao += 1
        #if ciao == 20000:
            #break
        #print "ele32: ", v_passHltEle32[0], "  EGL1: ", v_passEGL1[0], "  MVA: ", v_passMVA[0], "  pt: ", v_el_pt[0], "  eta: ", v_el_eta[0]
        if not (v_passHltEle32[0] == 1 and v_passEGL1[0] == 1 and v_passMVA[0] == 1 and v_el_pt[0] > 33. and math.fabs(v_el_eta[0]) < 2.50):
        #if not (entry.passHltEle27WPTightGsf and entry.passingMVA94Xwp80isoV2 and entry.el_pt > 28 and math.fabs(entry.el_eta) < 2.50):
            continue
        #print v_el_pt[0]

        if not (v_Tag_el_pt[0] > 35. and math.fabs(v_Tag_sc_eta[0]) <= 2.1):
            continue

        if not (v_mass[0] > 60. and v_mass[0] < 120.):
            continue

        #el_eta = entry.el_eta
        #abs_el_eta = math.fabs(el_eta)
        #fBrem = entry.el_fbrem
        if not v_el_fbrem[0] > 0.:
            continue
        #xOverX0 = -(math.log(1-fBrem))
        #rounded_abs_el_eta = round(abs_el_eta,2) #Round to 2 decimal places
        #rounded_el_eta = round(v_el_eta[0],2) #Round to 2 decimal places
        
        index = -9999.
        for i in np.arange(-2.50,2.50,0.05):
            if (v_el_eta[0] - i) >= 0. and (v_el_eta[0] - i) < 0.05:# or (rounded_abs_el_eta - i) == 0.:
                index = i
        if index == -9999.:
            print rounded_el_eta
            print "SOMETHING WENT WRONG"
            break
        
        if tree == treeIn_DATA:
            #h_dict_DATA[index].Fill(xOverX0)
            h_dict_DATA[index].Fill(v_el_fbrem[0])
        else:
            #h_dict_MC[index].Fill(xOverX0,entry.totWeight)
            h_dict_MC[index].Fill(v_el_fbrem[0],v_totWeight[0])

        if math.fabs(v_el_eta[0]) < 1.2:
            if tree == treeIn_DATA:
                h_fBrem_B_DATA.Fill(v_el_fbrem[0])
                h_mass_B_DATA.Fill(v_mass[0])
                if v_el_q[0] == 1:
                    h_fBrem_Plus_B_DATA.Fill(v_el_fbrem[0])
                else:
                    h_fBrem_Minus_B_DATA.Fill(v_el_fbrem[0])
            else:
                h_fBrem_B_MC.Fill(v_el_fbrem[0],v_totWeight[0])
                h_mass_B_MC.Fill(v_mass[0],v_totWeight[0])
                if v_el_q[0] == 1:
                    h_fBrem_Plus_B_MC.Fill(v_el_fbrem[0],v_totWeight[0])
                else:
                    h_fBrem_Minus_B_MC.Fill(v_el_fbrem[0],v_totWeight[0])

        if v_el_eta[0] <= -1.2 or v_el_eta[0] >= 1.2:
            if tree == treeIn_DATA:
                h_fBrem_E_DATA.Fill(v_el_fbrem[0])
                h_mass_E_DATA.Fill(v_mass[0])
                if v_el_q[0] == 1:
                    h_fBrem_Plus_E_DATA.Fill(v_el_fbrem[0])
                else:
                    h_fBrem_Minus_E_DATA.Fill(v_el_fbrem[0])
            else:
                h_fBrem_E_MC.Fill(v_el_fbrem[0],v_totWeight[0])
                h_mass_E_MC.Fill(v_mass[0],v_totWeight[0])
                if v_el_q[0] == 1:
                    h_fBrem_Plus_E_MC.Fill(v_el_fbrem[0],v_totWeight[0])
                else:
                    h_fBrem_Minus_E_MC.Fill(v_el_fbrem[0],v_totWeight[0])
        
        if tree == treeIn_DATA:
            h2_fBrem_phi_DATA.Fill(v_el_phi[0],v_el_fbrem[0])
            h2_fBrem_pT_DATA.Fill(v_el_pt[0],v_el_fbrem[0])
            h2_fBrem_eta_DATA.Fill(v_el_eta[0],v_el_fbrem[0])
        else:
            h2_fBrem_phi_MC.Fill(v_el_phi[0],v_el_fbrem[0],v_totWeight[0])
            h2_fBrem_pT_MC.Fill(v_el_pt[0],v_el_fbrem[0],v_totWeight[0])
            h2_fBrem_eta_MC.Fill(v_el_eta[0],v_el_fbrem[0],v_totWeight[0])


gStyle.SetOptStat(0)
fileOut_plots = TFile("histos/plots.root","RECREATE")
fileOut_plots.cd()
c_fBrem_B = TCanvas()
DATA_normalization_B = h_fBrem_B_DATA.Integral()
MC_normalization_B = h_fBrem_B_MC.Integral()
norm_factor_B = DATA_normalization_B/MC_normalization_B
h_fBrem_B_MC.Scale(norm_factor_B)
h_fBrem_B_MC.SetLineColor(800)
h_fBrem_B_MC.SetFillColor(800)
h_fBrem_B_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_B_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_B_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_B_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_B_MC.Draw("hist")
h_MCErr_fBrem_B = copy.deepcopy(h_fBrem_B_MC)
h_MCErr_fBrem_B.SetFillStyle(3008)
h_MCErr_fBrem_B.SetMarkerStyle(1)
h_MCErr_fBrem_B.SetFillColor(kOrange+7)
h_MCErr_fBrem_B.Draw("sameE2")
h_fBrem_B_DATA.SetMarkerStyle(20)
h_fBrem_B_DATA.Draw("SAME,PE")
c_fBrem_B.SaveAs("plots/fBrem_B.pdf")
h_fBrem_B_DATA.Write()
h_fBrem_B_MC.Write()

c_fBrem_E = TCanvas()
DATA_normalization_E = h_fBrem_E_DATA.Integral()
MC_normalization_E = h_fBrem_E_MC.Integral()
norm_factor_E = DATA_normalization_E/MC_normalization_E
h_fBrem_E_MC.Scale(norm_factor_E)
h_fBrem_E_MC.SetLineColor(800)
h_fBrem_E_MC.SetFillColor(800)
h_fBrem_E_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_E_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_E_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_E_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_E_MC.Draw("hist")
h_MCErr_fBrem_E = copy.deepcopy(h_fBrem_E_MC)
h_MCErr_fBrem_E.SetFillStyle(3008)
h_MCErr_fBrem_E.SetMarkerStyle(1)
h_MCErr_fBrem_E.SetFillColor(kOrange+7)
h_MCErr_fBrem_E.Draw("sameE2")
h_fBrem_E_DATA.SetMarkerStyle(20)
h_fBrem_E_DATA.Draw("SAME,PE")
c_fBrem_E.SaveAs("plots/fBrem_E.pdf")
h_fBrem_B_DATA.Write()
h_fBrem_B_MC.Write()

c_fBrem_Plus_B = TCanvas()
DATA_normalization_Plus_B = h_fBrem_Plus_B_DATA.Integral()
MC_normalization_Plus_B = h_fBrem_Plus_B_MC.Integral()
norm_factor_Plus_B = DATA_normalization_Plus_B/MC_normalization_Plus_B
h_fBrem_Plus_B_MC.Scale(norm_factor_Plus_B)
h_fBrem_Plus_B_MC.SetLineColor(800)
h_fBrem_Plus_B_MC.SetFillColor(800)
h_fBrem_Plus_B_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_Plus_B_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_Plus_B_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_Plus_B_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_Plus_B_MC.Draw("hist")
h_MCErr_fBrem_Plus_B = copy.deepcopy(h_fBrem_Plus_B_MC)
h_MCErr_fBrem_Plus_B.SetFillStyle(3008)
h_MCErr_fBrem_Plus_B.SetMarkerStyle(1)
h_MCErr_fBrem_Plus_B.SetFillColor(kOrange+7)
h_MCErr_fBrem_Plus_B.Draw("sameE2")
h_fBrem_Plus_B_DATA.SetMarkerStyle(20)
h_fBrem_Plus_B_DATA.Draw("SAME,PE")
c_fBrem_Plus_B.SaveAs("plots/fBrem_Plus_B.pdf")
h_fBrem_Plus_B_DATA.Write()
h_fBrem_Plus_B_MC.Write()

c_fBrem_Minus_B = TCanvas()
DATA_normalization_Minus_B = h_fBrem_Minus_B_DATA.Integral()
MC_normalization_Minus_B = h_fBrem_Minus_B_MC.Integral()
norm_factor_Minus_B = DATA_normalization_Minus_B/MC_normalization_Minus_B
h_fBrem_Minus_B_MC.Scale(norm_factor_Minus_B)
h_fBrem_Minus_B_MC.SetLineColor(800)
h_fBrem_Minus_B_MC.SetFillColor(800)
h_fBrem_Minus_B_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_Minus_B_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_Minus_B_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_Minus_B_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_Minus_B_MC.Draw("hist")
h_MCErr_fBrem_Minus_B = copy.deepcopy(h_fBrem_Minus_B_MC)
h_MCErr_fBrem_Minus_B.SetFillStyle(3008)
h_MCErr_fBrem_Minus_B.SetMarkerStyle(1)
h_MCErr_fBrem_Minus_B.SetFillColor(kOrange+7)
h_MCErr_fBrem_Minus_B.Draw("sameE2")
h_fBrem_Minus_B_DATA.SetMarkerStyle(20)
h_fBrem_Minus_B_DATA.Draw("SAME,PE")
c_fBrem_Minus_B.SaveAs("plots/fBrem_Minus_B.pdf")
h_fBrem_Minus_B_DATA.Write()
h_fBrem_Minus_B_MC.Write()

c_fBrem_Plus_E = TCanvas()
DATA_normalization_Plus_E = h_fBrem_Plus_E_DATA.Integral()
MC_normalization_Plus_E = h_fBrem_Plus_E_MC.Integral()
norm_factor_Plus_E = DATA_normalization_Plus_E/MC_normalization_Plus_E
h_fBrem_Plus_E_MC.Scale(norm_factor_Plus_E)
h_fBrem_Plus_E_MC.SetLineColor(800)
h_fBrem_Plus_E_MC.SetFillColor(800)
h_fBrem_Plus_E_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_Plus_E_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_Plus_E_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_Plus_E_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_Plus_E_MC.Draw("hist")
h_MCErr_fBrem_Plus_E = copy.deepcopy(h_fBrem_Plus_E_MC)
h_MCErr_fBrem_Plus_E.SetFillStyle(3008)
h_MCErr_fBrem_Plus_E.SetMarkerStyle(1)
h_MCErr_fBrem_Plus_E.SetFillColor(kOrange+7)
h_MCErr_fBrem_Plus_E.Draw("sameE2")
h_fBrem_Plus_E_DATA.SetMarkerStyle(20)
h_fBrem_Plus_E_DATA.Draw("SAME,PE")
c_fBrem_Plus_E.SaveAs("plots/fBrem_Plus_E.pdf")
h_fBrem_Plus_E_DATA.Write()
h_fBrem_Plus_E_MC.Write()

c_fBrem_Minus_E = TCanvas()
DATA_normalization_Minus_E = h_fBrem_Minus_E_DATA.Integral()
MC_normalization_Minus_E = h_fBrem_Minus_E_MC.Integral()
norm_factor_Minus_E = DATA_normalization_Minus_E/MC_normalization_Minus_E
h_fBrem_Minus_E_MC.Scale(norm_factor_Minus_E)
h_fBrem_Minus_E_MC.SetLineColor(800)
h_fBrem_Minus_E_MC.SetFillColor(800)
h_fBrem_Minus_E_MC.GetXaxis().SetTitle("f_{brem}")
h_fBrem_Minus_E_MC.GetXaxis().SetTitleSize(0.05)
h_fBrem_Minus_E_MC.GetXaxis().SetTitleOffset(0.7)
h_fBrem_Minus_E_MC.GetYaxis().SetTitle("Events/0.04")
h_fBrem_Minus_E_MC.Draw("hist")
h_MCErr_fBrem_Minus_E = copy.deepcopy(h_fBrem_Minus_E_MC)
h_MCErr_fBrem_Minus_E.SetFillStyle(3008)
h_MCErr_fBrem_Minus_E.SetMarkerStyle(1)
h_MCErr_fBrem_Minus_E.SetFillColor(kOrange+7)
h_MCErr_fBrem_Minus_E.Draw("sameE2")
h_fBrem_Minus_E_DATA.SetMarkerStyle(20)
h_fBrem_Minus_E_DATA.Draw("SAME,PE")
c_fBrem_Minus_E.SaveAs("plots/fBrem_Minus_E.pdf")
h_fBrem_Minus_E_DATA.Write()
h_fBrem_Minus_E_MC.Write()

c_mass_B = TCanvas()
DATA_normalization_mass_B = h_mass_B_DATA.Integral()
MC_normalization_mass_B = h_mass_B_MC.Integral()
norm_factor_mass_B = DATA_normalization_mass_B/MC_normalization_mass_B
h_mass_B_MC.Scale(norm_factor_mass_B)
h_mass_B_MC.SetLineColor(800)
h_mass_B_MC.SetFillColor(800)
h_mass_B_MC.GetXaxis().SetTitle("#it{m}_{ee}")
h_mass_B_MC.GetXaxis().SetTitleSize(0.05)
h_mass_B_MC.GetXaxis().SetTitleOffset(0.7)
h_mass_B_MC.GetYaxis().SetTitle("Events/2.0 GeV")
h_mass_B_MC.Draw("hist")
h_MCErr_mass_B = copy.deepcopy(h_mass_B_MC)
h_MCErr_mass_B.SetFillStyle(3008)
h_MCErr_mass_B.SetMarkerStyle(1)
h_MCErr_mass_B.SetFillColor(kOrange+7)
h_MCErr_mass_B.Draw("sameE2")
h_mass_B_DATA.SetMarkerStyle(20)
h_mass_B_DATA.Draw("SAME,PE")
c_mass_B.SaveAs("plots/mass_B.pdf")
h_mass_B_DATA.Write()
h_mass_B_MC.Write()

c_mass_E = TCanvas()
DATA_normalization_mass_E = h_mass_E_DATA.Integral()
MC_normalization_mass_E = h_mass_E_MC.Integral()
norm_factor_mass_E = DATA_normalization_mass_E/MC_normalization_mass_E
h_mass_E_MC.Scale(norm_factor_mass_E)
h_mass_E_MC.SetLineColor(800)
h_mass_E_MC.SetFillColor(800)
h_mass_E_MC.GetXaxis().SetTitle("#it{m}_{ee}")
h_mass_E_MC.GetXaxis().SetTitleSize(0.05)
h_mass_E_MC.GetXaxis().SetTitleOffset(0.7)
h_mass_E_MC.GetYaxis().SetTitle("Events/2.0 GeV")
h_mass_E_MC.Draw("hist")
h_MCErr_mass_E = copy.deepcopy(h_mass_E_MC)
h_MCErr_mass_E.SetFillStyle(3008)
h_MCErr_mass_E.SetMarkerStyle(1)
h_MCErr_mass_E.SetFillColor(kOrange+7)
h_MCErr_mass_E.Draw("sameE2")
h_mass_E_DATA.SetMarkerStyle(20)
h_mass_E_DATA.Draw("SAME,PE")
c_mass_E.SaveAs("plots/mass_E.pdf")
h_mass_B_DATA.Write()
h_mass_B_MC.Write()

c_fBrem_vs_phi_DATA = TCanvas()
h2_fBrem_phi_DATA.GetXaxis().SetTitle("#phi^{e}")
h2_fBrem_phi_DATA.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_phi_DATA.Draw("COLZ")
c_fBrem_vs_phi_DATA.SaveAs("plots/fBrem_vs_phi_DATA.pdf")
h2_fBrem_phi_DATA.Write()

c_fBrem_vs_phi_MC = TCanvas()
h2_fBrem_phi_MC.GetXaxis().SetTitle("#phi^{e}")
h2_fBrem_phi_MC.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_phi_MC.Draw("COLZ")
c_fBrem_vs_phi_MC.SaveAs("plots/fBrem_vs_phi_MC.pdf")
h2_fBrem_phi_MC.Write()

c_fBrem_vs_pT_DATA = TCanvas()
h2_fBrem_pT_DATA.GetXaxis().SetTitle("#it{p}_{T}^{e} (GeV)")
h2_fBrem_pT_DATA.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_pT_DATA.Draw("COLZ")
c_fBrem_vs_pT_DATA.SaveAs("plots/fBrem_vs_pT_DATA.pdf")
h2_fBrem_pT_DATA.Write()

c_fBrem_vs_pT_MC = TCanvas()
h2_fBrem_pT_MC.GetXaxis().SetTitle("#it{p}_{T}^{e} (GeV)")
h2_fBrem_pT_MC.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_pT_MC.Draw("COLZ")
c_fBrem_vs_pT_MC.SaveAs("plots/fBrem_vs_pT_MC.pdf")
h2_fBrem_pT_MC.Write()

c_fBrem_vs_eta_DATA = TCanvas()
h2_fBrem_eta_DATA.GetXaxis().SetTitle("#eta^{e}")
h2_fBrem_eta_DATA.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_eta_DATA.Draw("COLZ")
c_fBrem_vs_eta_DATA.SaveAs("plots/fBrem_vs_eta_DATA.pdf")
h2_fBrem_eta_DATA.Write()

c_fBrem_vs_eta_MC = TCanvas()
h2_fBrem_eta_MC.GetXaxis().SetTitle("#eta^{e}")
h2_fBrem_eta_MC.GetYaxis().SetTitle("f_{brem}")
h2_fBrem_eta_MC.Draw("COLZ")
c_fBrem_vs_eta_MC.SaveAs("plots/fBrem_vs_eta_MC.pdf")
h2_fBrem_eta_MC.Write()

# fileOut = TFile("histos/xOverX0_CORR.root","RECREATE")
# fileOut.cd()
# gDirectory.mkdir("DATA")
# fileOut.cd("DATA")
# for i in np.arange(0.,2.55,0.05):
#     h_dict_DATA[i].Write()

# fileOut.cd()
# gDirectory.mkdir("MC")
# fileOut.cd("MC")
# for i in np.arange(0.,2.55,0.05):
#     h_dict_MC[i].Write()

fileOut = TFile("histos/fBrem_notFolded_full2017.root","RECREATE")
fileOut.cd()
gDirectory.mkdir("DATA")
fileOut.cd("DATA")
for i in np.arange(-2.50,2.50,0.05):
    h_dict_DATA[i].Write()

fileOut.cd()
gDirectory.mkdir("MC")
fileOut.cd("MC")
for i in np.arange(-2.50,2.50,0.05):
    h_dict_MC[i].Write()

raw_input()
