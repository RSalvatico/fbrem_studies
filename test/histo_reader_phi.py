from ROOT import TFile, TH1F, TCanvas, TGraph, TGraphErrors, TLegend, TPad, TLine, kGreen
import math
import numpy as np
from array import array
import sys

fileIn = TFile("histos/17_06_2021/fBrem_etaFolded_EtaPhiBins.root")
h_dict_DATA = dict()
h_dict_MC = dict()
h_dict_energy_DATA = dict()
h_dict_energy_MC = dict()

#Some arrays for TGraphErrors                                                                                                                         
x = array( 'f',[])
x_err = array( 'f',[])
y_DATA = array( 'f',[])
y_err_DATA = array( 'f',[])
y_MC = array( 'f',[])
y_err_MC = array( 'f',[])
y_diff = array( 'f',[])
y_err_diff = array( 'f',[])

y_fBrem_asymmetry_DATA = array( 'f',[])
y_fBrem_asymmetry_MC   = array( 'f',[])
y_fBrem_asymmetry_err_DATA = array( 'f',[])
y_fBrem_asymmetry_err_MC = array( 'f',[])
x_fBrem_asymmetry = array( 'f',[])
x_fBrem_asymmetry_err = array( 'f',[])
dict_fBrem_DATA = dict()
dict_fBrem_MC   = dict()
dict_fBrem_err_DATA = dict()
dict_fBrem_err_MC   = dict()


calibration_file = TFile("Calibration_t1.root")


for i in np.arange(-2.50,0.00,0.05):

    str_i = '%.6f' % i
    str_iplus = '%.6f' %(i+0.05)

    for p in np.arange(-3.14,3.14,3.14):
    
        str_p = '%.6f' % p
        str_pplus = '%.6f' %(p+3.14)

        if str_p == "0.000000": str_p = "-0.000000"

        h_dict_DATA[str_i+str_p] = fileIn.Get("DATA/DATA_h_fBrem_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)
        h_dict_MC[str_i+str_p]   = fileIn.Get("MC/MC_h_fBrem_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)

        h_dict_energy_DATA[i] = fileIn.Get("DATA/DATA_h_energy_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)
        h_dict_energy_MC[i]   = fileIn.Get("MC/MC_h_energy_etaBin_" + str_i + "_" + str_iplus + "_phiBin_" + str_p + "_" + str_pplus)

        mean_DATA = h_dict_DATA[str_i+str_p].GetMean()
        mean_MC   = h_dict_MC[str_i+str_p].GetMean()
        uncert_mean_DATA = h_dict_DATA[str_i+str_p].GetStdDev()/math.sqrt(h_dict_DATA[str_i+str_p].GetEntries()) #std.dev/sqrt(N)
        uncert_mean_MC   = h_dict_MC[str_i+str_p].GetStdDev()/math.sqrt(h_dict_MC[str_i+str_p].GetEntries()) #std.dev/sqrt(N)

        # xOverX0_DATA = -(math.log(1 - mean_DATA))
        # xOverX0_MC   = -(math.log(1 - mean_MC))
        # uncert_xOverX0_DATA = uncert_mean_DATA/(1. - mean_DATA) #error propagation
        # uncert_xOverX0_MC   = uncert_mean_MC/(1. - mean_MC) #error propagation

        # mean_energy_DATA = h_dict_energy_DATA[i].GetMean()
        # mean_energy_MC   = h_dict_energy_MC[i].GetMean()
        
        # #Choose the most appropriate calibration file
        # if math.fabs(mean_energy_DATA - 25.) < math.fabs(mean_energy_DATA - 50.):
        #     gr_t1 = calibration_file.Get("Calibration_graph_t1_25_GeV")
        # elif math.fabs(mean_energy_DATA - 25.) > math.fabs(mean_energy_DATA - 50.) and math.fabs(mean_energy_DATA - 50.) < math.fabs(mean_energy_DATA - 100.):
        #     gr_t1 = calibration_file.Get("Calibration_graph_t1_50_GeV")
        # elif math.fabs(mean_energy_DATA - 50.) > math.fabs(mean_energy_DATA - 100.) and math.fabs(mean_energy_DATA - 100.) < math.fabs(mean_energy_DATA - 200.):
        #     gr_t1 = calibration_file.Get("Calibration_graph_t1_100_GeV")
        # elif math.fabs(mean_energy_DATA - 100.) > math.fabs(mean_energy_DATA - 200.):
        #     gr_t1 = calibration_file.Get("Calibration_graph_t1_200_GeV")

        # dict_calib_xOverX0_DATA[str_i+str_p] = gr_t1.Eval(xOverX0_DATA)
        # dict_calib_xOverX0_MC[str_i+str_p]   = gr_t1.Eval(xOverX0_MC)
        
        if str_p == "-3.140000":
            dict_fBrem_DATA[str_i+str_p] = mean_DATA 
            dict_fBrem_MC[str_i+str_p] = mean_MC
            dict_fBrem_err_DATA[str_i+str_p] = uncert_mean_DATA
            dict_fBrem_err_MC[str_i+str_p] = uncert_mean_MC

        if str_p == "-0.000000":
            dict_fBrem_DATA[str_i+"-3.140000"] = (mean_DATA - dict_fBrem_DATA[str_i+"-3.140000"])/(mean_DATA + dict_fBrem_DATA[str_i+"-3.140000"])       
            dict_fBrem_MC[str_i+"-3.140000"] = (mean_MC - dict_fBrem_MC[str_i+"-3.140000"])/(mean_MC + dict_fBrem_MC[str_i+"-3.140000"])

            dict_fBrem_err_DATA[str_i+"-3.140000"] = (2./(dict_fBrem_DATA[str_i+"-3.140000"] + mean_DATA)**2)*math.sqrt((mean_DATA*dict_fBrem_err_DATA[str_i+"-3.140000"])**2 + ( dict_fBrem_DATA[str_i+"-3.140000"]*uncert_mean_DATA)**2)
            dict_fBrem_err_MC[str_i+"-3.140000"] = (2./(dict_fBrem_MC[str_i+"-3.140000"] + mean_MC)**2)*math.sqrt((mean_MC*dict_fBrem_err_MC[str_i+"-3.140000"])**2 + ( dict_fBrem_MC[str_i+"-3.140000"]*uncert_mean_MC)**2)



for l in np.arange(-2.50,0.00,0.05):
    str_i = '%.6f' % l
    #str_i = str_i.replace("0000","")
    #str_i = str_i.replace("-","")
    x_fBrem_asymmetry.append(-l-0.025)
    x_fBrem_asymmetry_err.append(0.)
    y_fBrem_asymmetry_DATA.append(dict_fBrem_DATA[str_i+"-3.140000"])
    y_fBrem_asymmetry_MC.append(dict_fBrem_MC[str_i+"-3.140000"])
    y_fBrem_asymmetry_err_DATA.append(dict_fBrem_err_DATA[str_i+"-3.140000"])
    y_fBrem_asymmetry_err_MC.append(dict_fBrem_err_MC[str_i+"-3.140000"])
    
    #print "dict_fBrem_DATA[i]: ", dict_fBrem_DATA[str_i], "  dict_fBrem_MC[i]: ", dict_fBrem_MC[str_i]
    y_diff.append(dict_fBrem_DATA[str_i+"-3.140000"] - dict_fBrem_MC[str_i+"-3.140000"])
    y_err_diff.append(math.sqrt(dict_fBrem_err_DATA[str_i+"-3.140000"]**2 + dict_fBrem_err_MC[str_i+"-3.140000"]**2))

    
c_asymmetry = TCanvas()
c_asymmetry.SetGrid()

pad1 = TPad("pad_values","",0,0.28,1,1.)
pad2 = TPad("pad_ratio",'',0,0,1,0.25)
pad1.SetBottomMargin(0.02)
pad1.SetBorderMode(0)
pad1.SetBorderSize(0)
pad1.SetFrameBorderSize(0)
pad2.SetBorderSize(0)
pad2.SetFrameBorderSize(0)
pad2.SetBottomMargin(0.3)
pad2.SetBorderMode(0)
pad1.SetGrid()
pad1.Draw()
pad2.Draw()

pad1.cd()

g_fBrem_asymmetry_DATA = TGraphErrors(len(x_fBrem_asymmetry), x_fBrem_asymmetry, y_fBrem_asymmetry_DATA, x_fBrem_asymmetry_err, y_fBrem_asymmetry_err_DATA)
g_fBrem_asymmetry_DATA.SetMarkerColor(2)
g_fBrem_asymmetry_DATA.SetMarkerStyle(20)
g_fBrem_asymmetry_DATA.SetLineColor(2)
g_fBrem_asymmetry_DATA.SetTitle("")
g_fBrem_asymmetry_DATA.GetXaxis().SetTitle("|#eta^{e}|")
g_fBrem_asymmetry_DATA.GetXaxis().SetLabelSize(0)
g_fBrem_asymmetry_DATA.GetYaxis().SetTitle("(f_{brem}^{+} - f_{brem}^{-}) / (f_{brem}^{+} + f_{brem}^{-})")
g_fBrem_asymmetry_DATA.GetYaxis().SetTitleOffset(1.0)
g_fBrem_asymmetry_DATA.GetYaxis().SetTitleSize(0.05)
g_fBrem_asymmetry_DATA.GetYaxis().SetRangeUser(-0.05,0.05)
g_fBrem_asymmetry_DATA.Draw("AP")

g_fBrem_asymmetry_MC = TGraphErrors(len(x_fBrem_asymmetry), x_fBrem_asymmetry, y_fBrem_asymmetry_MC, x_fBrem_asymmetry_err, y_fBrem_asymmetry_err_MC)
g_fBrem_asymmetry_MC.SetMarkerColor(kGreen+1)
g_fBrem_asymmetry_MC.SetMarkerStyle(24)
g_fBrem_asymmetry_MC.SetLineColor(kGreen+1)
g_fBrem_asymmetry_MC.SetTitle("")
g_fBrem_asymmetry_MC.GetXaxis().SetTitle("|#eta^{e}|")
g_fBrem_asymmetry_MC.GetXaxis().SetLabelSize(0)
g_fBrem_asymmetry_MC.GetYaxis().SetTitle("(f_{brem}^{+} - f_{brem}^{-}) / (f_{brem}^{+} + f_{brem}^{-})")
g_fBrem_asymmetry_MC.GetYaxis().SetTitleOffset(1.0)
g_fBrem_asymmetry_MC.GetYaxis().SetTitleSize(0.05)
g_fBrem_asymmetry_MC.GetYaxis().SetRangeUser(-0.05,0.05)
g_fBrem_asymmetry_MC.Draw("SAME,P")

legend_asymmetry = TLegend(0.7,0.7,0.85,0.95)
legend_asymmetry.SetHeader(" ")
legend_asymmetry.SetFillColor(0)
legend_asymmetry.SetBorderSize(0)
legend_asymmetry.SetLineColor(1)
legend_asymmetry.SetLineStyle(1)
legend_asymmetry.SetLineWidth(1)
legend_asymmetry.SetFillStyle(0)
legend_asymmetry.AddEntry(g_fBrem_asymmetry_DATA,"Data","lep")
legend_asymmetry.AddEntry(g_fBrem_asymmetry_MC,"MC","lep")
legend_asymmetry.Draw("SAME")

pad2.cd()
pad2.SetTopMargin(0.03)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)

g_diff = TGraphErrors(len(x_fBrem_asymmetry), x_fBrem_asymmetry, y_diff, x_fBrem_asymmetry_err, y_err_diff)
g_diff.SetMarkerStyle(20)
g_diff.SetMarkerSize(0.5)
g_diff.SetTitle("")
g_diff.GetXaxis().SetTitle("|#eta^{e}|")
g_diff.GetYaxis().SetTitle("Data-MC")
#g_diff.GetYaxis().SetRangeUser(-0.2,0.2)
g_diff.GetYaxis().SetRangeUser(-0.05,0.05)
g_diff.GetXaxis().SetLabelSize(0.13)
g_diff.GetXaxis().SetTitleSize(0.15)
g_diff.GetXaxis().SetTitleOffset(0.7)
g_diff.GetYaxis().SetTitleSize(0.14)
g_diff.GetYaxis().SetTitleOffset(0.35)
g_diff.GetYaxis().SetLabelSize(0.10)
g_diff.GetYaxis().SetNdivisions(502,False)
line_on_zero = TLine(g_fBrem_asymmetry_DATA.GetXaxis().GetXmin(),0.,g_fBrem_asymmetry_DATA.GetXaxis().GetXmax(),0.)
line_on_zero.SetLineColor(4)
line_on_zero.SetLineStyle(2)
g_diff.Draw("AP")
line_on_zero.Draw("SAME")


c_asymmetry.SaveAs("plots/fBrem_phi_asymmetry.pdf")

raw_input()
