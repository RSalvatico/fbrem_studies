from ROOT import TFile, TH1F, TCanvas, TGraphErrors, TLegend, TPad, TLine, kOrange, kGreen
import math
import numpy as np
from array import array
import sys
import copy

fileIn = TFile("histos/17_06_2021/fBrem_etaFolded_EtaPhiBins.root")

h_dict_DATA = dict()
h_dict_MC = dict()
c_dict = dict()

h_dict_xOverX0_DATA = dict()
h_dict_xOverX0_MC = dict()
c_dict_xOverX0 = dict()

for i in np.arange(-2.50,0.00,0.05):
    str_i = '%.6f' % i
    str_iplus = '%.6f' %(i+0.05)
    #str_i = str_i.replace(".","_")
    #str_iplus = str_iplus.replace(".","_")
    
    for p in np.arange(-3.14,3.14,3.14):
    
        str_p = '%.6f' % p
        str_pplus = '%.6f' %(p+3.14)

        if str_p == "0.000000": str_p = "-0.000000"

        h_dict_DATA[str_i+str_p] = fileIn.Get("DATA/DATA_h_fBrem_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)
        h_dict_MC[str_i+str_p]   = fileIn.Get("MC/MC_h_fBrem_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)
        
        h_dict_xOverX0_DATA[str_i+str_p] = fileIn.Get("DATA/DATA_h_xOverX0_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)
        h_dict_xOverX0_MC[str_i+str_p]   = fileIn.Get("MC/MC_h_xOverX0_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)

        print "DATA/DATA_h_fBrem_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus

        #c_dict[str_i+str_p] = TCanvas()
        DATA_normalization = h_dict_DATA[str_i+str_p].Integral()
        MC_normalization = h_dict_MC[str_i+str_p].Integral()
        norm_factor = DATA_normalization/MC_normalization
        h_dict_MC[str_i+str_p].Scale(norm_factor)
        h_dict_MC[str_i+str_p].SetLineColor(800)
        h_dict_MC[str_i+str_p].SetFillColor(800)
        h_dict_MC[str_i+str_p].GetXaxis().SetTitle("f_{brem}")
        h_dict_MC[str_i+str_p].GetXaxis().SetTitleSize(0.05)
        h_dict_MC[str_i+str_p].GetXaxis().SetTitleOffset(0.7)
        h_dict_MC[str_i+str_p].GetYaxis().SetTitle("Events/0.01")
        h_dict_MC[str_i+str_p].GetYaxis().SetTitleOffset(1.5)
        h_dict_MC[str_i+str_p].Draw("hist")
        h_MCErr = copy.deepcopy(h_dict_MC[str_i+str_p])
        h_MCErr.SetFillStyle(3008)
        h_MCErr.SetMarkerStyle(1)
        h_MCErr.SetFillColor(kOrange-7)
        h_MCErr.Draw("sameE2")
        h_dict_DATA[str_i+str_p].SetMarkerStyle(20)
        h_dict_DATA[str_i+str_p].Draw("SAME,PE")
        
        #c_dict[str_i+str_p].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/fBrem_etaBins/full2017/CutBased_custom/h_fBrem_etaBins_" + str_i + "_" + str_iplus + ".png")
        
        mean_DATA = h_dict_DATA[str_i+str_p].GetMean()
        mean_MC   = h_dict_MC[str_i+str_p].GetMean()
        uncert_mean_DATA = h_dict_DATA[str_i+str_p].GetStdDev()/math.sqrt(h_dict_DATA[str_i+str_p].GetEntries()) #std.dev/sqrt(N)
        uncert_mean_MC   = h_dict_MC[str_i+str_p].GetStdDev()/math.sqrt(h_dict_MC[str_i+str_p].GetEntries()) #std.dev/sqrt(N)
        
        xOverX0_DATA = -(math.log(1 - mean_DATA))
        xOverX0_MC   = -(math.log(1 - mean_MC))
        uncert_xOverX0_DATA = uncert_mean_DATA/(1. - mean_DATA) #error propagation
        uncert_xOverX0_MC   = uncert_mean_MC/(1. - mean_MC) #error propagation
        
        c_dict_xOverX0[str_i+str_p] = TCanvas()
        DATA_xOverX0_normalization = h_dict_xOverX0_DATA[str_i+str_p].Integral()
        MC_xOverX0_normalization = h_dict_xOverX0_MC[str_i+str_p].Integral()
        norm_factor_xOverX0 = DATA_xOverX0_normalization/MC_xOverX0_normalization
        h_dict_xOverX0_MC[str_i+str_p].Scale(norm_factor_xOverX0)
        h_dict_xOverX0_MC[str_i+str_p].SetLineColor(800)
        h_dict_xOverX0_MC[str_i+str_p].SetFillColor(800)
        h_dict_xOverX0_MC[str_i+str_p].GetXaxis().SetTitle("x/X0")
        h_dict_xOverX0_MC[str_i+str_p].GetXaxis().SetTitleSize(0.05)
        h_dict_xOverX0_MC[str_i+str_p].GetXaxis().SetTitleOffset(0.7)
        h_dict_xOverX0_MC[str_i+str_p].GetYaxis().SetTitle("Events/0.025")
        h_dict_xOverX0_MC[str_i+str_p].GetYaxis().SetTitleOffset(1.5)
        h_dict_xOverX0_MC[str_i+str_p].Draw("hist")
        h_MCErr_xOverX0 = copy.deepcopy(h_dict_xOverX0_MC[str_i+str_p])
        h_MCErr_xOverX0.SetFillStyle(3008)
        h_MCErr_xOverX0.SetMarkerStyle(1)
        h_MCErr_xOverX0.SetFillColor(kOrange-7)
        h_MCErr_xOverX0.Draw("sameE2")
        h_dict_xOverX0_DATA[str_i+str_p].SetMarkerStyle(20)
        h_dict_xOverX0_DATA[str_i+str_p].Draw("SAME,PE")
        
        h_line_DATA = TLine(xOverX0_DATA,0,xOverX0_DATA,h_dict_xOverX0_DATA[str_i+str_p].GetMaximum())
        h_line_MC   = TLine(xOverX0_MC,0,xOverX0_MC,h_dict_xOverX0_DATA[str_i+str_p].GetMaximum())
        h_line_DATA.SetLineColor(2)
        h_line_MC.SetLineColor(kGreen+1)
        h_line_DATA.Draw("SAME")
        h_line_MC.Draw("SAME")
        
        #Remove some unpleasant zeros from the histo names
        n_str_i = str_i.replace("0000","")
        n_str_iplus = str_iplus.replace("0000","")
        n_str_p = str_p.replace("0000","")
        n_str_pplus = str_pplus.replace("0000","")
        
        c_dict_xOverX0[str_i+str_p].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/xOverX0_EtaPhiBins/17_06_2021/h_xOverX0_etaBin_" + n_str_i + "_" + n_str_iplus + "_phiBin_" + n_str_p + "_" + n_str_pplus + ".png")




