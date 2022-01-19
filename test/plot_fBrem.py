from ROOT import TFile, TH1F, TCanvas, TGraphErrors, TLegend, TPad, TLine, kOrange, kGreen, gStyle
import math
import numpy as np
from array import array
import sys
import copy

#fileIn = TFile("histos/02_06_2021/fBrem_notFolded_full2017_CutBased_DYNLO.root")
#fileIn = TFile("histos/28_06_2021/fBrem_notFolded_vs_PVz.root")
fileIn = TFile("histos/09_12_2021/fBrem_notFolded_postVFP_2016.root")

h_dict_DATA = dict()
h_dict_MC = dict()
c_dict = dict()

h_dict_xOverX0_DATA = dict()
h_dict_xOverX0_MC = dict()
c_dict_xOverX0 = dict()

h2_dict_fBrem_vs_PVz_DATA = dict()
h2_dict_fBrem_vs_PVz_MC = dict()
c_dict_fBrem_vs_PVz_DATA = dict()
c_dict_fBrem_vs_PVz_MC = dict()

def plot_fBrem():
    
    for i in np.arange(-2.50,2.50,0.05):
        str_i = '%.6f' % i
        str_iplus = '%.6f' %(i+0.05)
        #str_i = str_i.replace(".","_")
        #str_iplus = str_iplus.replace(".","_")
        
        h_dict_DATA[i] = fileIn.Get("DATA/DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
        h_dict_MC[i]   = fileIn.Get("MC/MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
        
        #Remove some unpleasant zeros from the histo names
        str_i = str_i.replace("0000","")
        str_iplus = str_iplus.replace("0000","")
        
        c_dict[i] = TCanvas()
        DATA_normalization = h_dict_DATA[i].Integral()
        MC_normalization = h_dict_MC[i].Integral()
        norm_factor = DATA_normalization/MC_normalization
        h_dict_MC[i].Scale(norm_factor)
        h_dict_MC[i].SetLineColor(800)
        h_dict_MC[i].SetFillColor(800)
        h_dict_MC[i].GetXaxis().SetTitle("f_{brem}")
        h_dict_MC[i].GetXaxis().SetTitleSize(0.05)
        h_dict_MC[i].GetXaxis().SetTitleOffset(0.7)
        h_dict_MC[i].GetYaxis().SetTitle("Events/0.01")
        h_dict_MC[i].GetYaxis().SetTitleOffset(1.5)
        h_dict_MC[i].Draw("hist")
        h_MCErr = copy.deepcopy(h_dict_MC[i])
        h_MCErr.SetFillStyle(3008)
        h_MCErr.SetMarkerStyle(1)
        h_MCErr.SetFillColor(kOrange-7)
        h_MCErr.Draw("sameE2")
        h_dict_DATA[i].SetMarkerStyle(20)
        h_dict_DATA[i].Draw("SAME,PE")
        
        c_dict[i].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/fBrem_etaBins/full2016/prePVz_recalibration/postVFP/h_fBrem_etaBins_" + str_i + "_" + str_iplus + ".png")

def plot_fBrem_vs_PVz():

    for i in np.arange(-2.50,2.50,0.05):
        str_i = '%.6f' % i
        str_iplus = '%.6f' %(i+0.05)
        
        h2_dict_fBrem_vs_PVz_DATA[i] = fileIn.Get("DATA/DATA_h2_fBrem_vs_PVz_etaBins_" + str_i + "_" + str_iplus)
        h2_dict_fBrem_vs_PVz_MC[i]   = fileIn.Get("MC/MC_h2_fBrem_vs_PVz_etaBins_" + str_i + "_" + str_iplus)
        
        #Remove some unpleasant zeros from the histo names
        str_i = str_i.replace("0000","")
        str_iplus = str_iplus.replace("0000","")


        corr_factor_DATA = h2_dict_fBrem_vs_PVz_DATA[i].GetCorrelationFactor()
        corr_factor_MC   = h2_dict_fBrem_vs_PVz_MC[i].GetCorrelationFactor()
        print str_i,"-",str_iplus
        print "cf DATA: ", corr_factor_DATA, "  cf MC: ", corr_factor_MC

        #gStyle.SetOptStat(0)
        c_dict_fBrem_vs_PVz_DATA[i] = TCanvas()
        h2_dict_fBrem_vs_PVz_DATA[i].GetXaxis().SetTitle("PV z (cm)")
        h2_dict_fBrem_vs_PVz_DATA[i].GetXaxis().SetTitleSize(0.05)
        h2_dict_fBrem_vs_PVz_DATA[i].GetXaxis().SetTitleOffset(0.7)
        h2_dict_fBrem_vs_PVz_DATA[i].GetYaxis().SetTitle("f_{brem}")
        h2_dict_fBrem_vs_PVz_DATA[i].GetYaxis().SetTitleOffset(1.0)
        h2_dict_fBrem_vs_PVz_DATA[i].GetYaxis().SetTitleSize(0.05)
        h2_dict_fBrem_vs_PVz_DATA[i].Draw("COLZ")
        #c_dict_fBrem_vs_PVz_DATA[i].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/fBrem_vs_PVz_etaBins/DATA/h_fBrem_vs_PVz_etaBins_" + str_i + "_" + str_iplus + ".png")
        
        c_dict_fBrem_vs_PVz_MC[i] = TCanvas()
        h2_dict_fBrem_vs_PVz_MC[i].GetXaxis().SetTitle("PV z (cm)")
        h2_dict_fBrem_vs_PVz_MC[i].GetXaxis().SetTitleSize(0.05)
        h2_dict_fBrem_vs_PVz_MC[i].GetXaxis().SetTitleOffset(0.7)
        h2_dict_fBrem_vs_PVz_MC[i].GetYaxis().SetTitle("f_{brem}")
        h2_dict_fBrem_vs_PVz_MC[i].GetYaxis().SetTitleOffset(1.0)
        h2_dict_fBrem_vs_PVz_MC[i].GetYaxis().SetTitleSize(0.05)
        h2_dict_fBrem_vs_PVz_MC[i].Draw("COLZ")
        #c_dict_fBrem_vs_PVz_MC[i].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/fBrem_vs_PVz_etaBins/MC/h_fBrem_vs_PVz_etaBins_" + str_i + "_" + str_iplus + ".png")

def plot_xOverX0():

    for i in np.arange(-2.50,2.50,0.05):
        str_i = '%.6f' % i
        str_iplus = '%.6f' %(i+0.05)
        #str_i = str_i.replace(".","_")
        #str_iplus = str_iplus.replace(".","_")
        
        h_dict_DATA[i] = fileIn.Get("DATA/DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
        h_dict_MC[i]   = fileIn.Get("MC/MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
        h_dict_xOverX0_DATA[i] = fileIn.Get("DATA/DATA_h_xOverX0_etaBins_" + str_i + "_" + str_iplus)
        h_dict_xOverX0_MC[i]   = fileIn.Get("MC/MC_h_xOverX0_etaBins_" + str_i + "_" + str_iplus)
        
        #Remove some unpleasant zeros from the histo names
        str_i = str_i.replace("0000","")
        str_iplus = str_iplus.replace("0000","")
        
        mean_DATA = h_dict_DATA[i].GetMean()
        mean_MC   = h_dict_MC[i].GetMean()
        uncert_mean_DATA = h_dict_DATA[i].GetStdDev()/math.sqrt(h_dict_DATA[i].GetEntries()) #std.dev/sqrt(N)
        uncert_mean_MC   = h_dict_MC[i].GetStdDev()/math.sqrt(h_dict_MC[i].GetEntries()) #std.dev/sqrt(N)
        
        xOverX0_DATA = -(math.log(1 - mean_DATA))
        xOverX0_MC   = -(math.log(1 - mean_MC))
        uncert_xOverX0_DATA = uncert_mean_DATA/(1. - mean_DATA) #error propagation
        uncert_xOverX0_MC   = uncert_mean_MC/(1. - mean_MC) #error propagation
        
        c_dict_xOverX0[i] = TCanvas()
        DATA_xOverX0_normalization = h_dict_xOverX0_DATA[i].Integral()
        MC_xOverX0_normalization = h_dict_xOverX0_MC[i].Integral()
        norm_factor_xOverX0 = DATA_xOverX0_normalization/MC_xOverX0_normalization
        h_dict_xOverX0_MC[i].Scale(norm_factor_xOverX0)
        h_dict_xOverX0_MC[i].SetLineColor(800)
        h_dict_xOverX0_MC[i].SetFillColor(800)
        h_dict_xOverX0_MC[i].GetXaxis().SetTitle("x/X0")
        h_dict_xOverX0_MC[i].GetXaxis().SetTitleSize(0.05)
        h_dict_xOverX0_MC[i].GetXaxis().SetTitleOffset(0.7)
        h_dict_xOverX0_MC[i].GetYaxis().SetTitle("Events/0.02")
        h_dict_xOverX0_MC[i].GetYaxis().SetTitleOffset(1.5)
        h_dict_xOverX0_MC[i].Draw("hist")
        h_MCErr_xOverX0 = copy.deepcopy(h_dict_xOverX0_MC[i])
        h_MCErr_xOverX0.SetFillStyle(3008)
        h_MCErr_xOverX0.SetMarkerStyle(1)
        h_MCErr_xOverX0.SetFillColor(kOrange-7)
        h_MCErr_xOverX0.Draw("sameE2")
        h_dict_xOverX0_DATA[i].SetMarkerStyle(20)
        h_dict_xOverX0_DATA[i].Draw("SAME,PE")

        print str_i, " ", str_iplus
        print "DATA: ", xOverX0_DATA, "\pm", uncert_xOverX0_DATA, "  MC: ", xOverX0_MC, "\pm", uncert_xOverX0_MC
        print
        h_line_DATA = TLine(xOverX0_DATA,0,xOverX0_DATA,h_dict_xOverX0_DATA[i].GetMaximum())
        h_line_MC   = TLine(xOverX0_MC,0,xOverX0_MC,h_dict_xOverX0_DATA[i].GetMaximum())
        h_line_DATA.SetLineColor(2)
        h_line_MC.SetLineColor(kGreen+1)
        h_line_DATA.Draw("SAME")
        h_line_MC.Draw("SAME")
        
        
        #c_dict_xOverX0[i].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/xOverX0_etaBins/full2017/02_01_2021/h_xOverX0_etaBins_" + str_i + "_" + str_iplus + ".png")

if __name__ == "__main__":
    
    plot_fBrem()
    #plot_fBrem_vs_PVz()
    #plot_xOverX0()