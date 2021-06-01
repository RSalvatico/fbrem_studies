from ROOT import TFile, TH1D, TH2D, TCanvas, TGraphErrors, TLegend, TPad, TLine, kOrange
import math
import numpy as np
from array import array
import sys
import copy

fileIn = TFile("histos/cpp/fBrem_notFolded_full2017_CutBased_DYNLO_custom.root")
h_dict_DATA = dict()
h_dict_MC = dict()
c_dict = dict()


for i in np.arange(-2.50,2.50,0.05):
    str_i = '%.6f' % i
    str_iplus = '%.6f' %(i+0.05)
    #str_i = str_i.replace(".","_")
    #str_iplus = str_iplus.replace(".","_")
    h_dict_DATA[i] = fileIn.Get("DATA/DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
    h_dict_MC[i]   = fileIn.Get("MC/MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
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
    h_dict_MC[i].GetYaxis().SetTitle("Events/0.02")
    h_dict_MC[i].GetYaxis().SetTitleOffset(1.5)
    h_dict_MC[i].Draw("hist")
    h_MCErr = copy.deepcopy(h_dict_MC[i])
    h_MCErr.SetFillStyle(3008)
    h_MCErr.SetMarkerStyle(1)
    h_MCErr.SetFillColor(kOrange-7)
    h_MCErr.Draw("sameE2")
    h_dict_DATA[i].SetMarkerStyle(20)
    h_dict_DATA[i].Draw("SAME,PE")

    str_i = str_i.replace("0000","")
    str_iplus = str_iplus.replace("0000","")

    c_dict[i].SaveAs("/afs/cern.ch/user/r/rselvati/www/fBrem/fBrem_etaBins/full2017/CutBased_custom/h_fBrem_etaBins_" + str_i + "_" + str_iplus + ".png")




