from ROOT import TFile, TH1F, TCanvas, TGraph, TGraphErrors, TLegend, TPad, TLine, kGreen
import math
import numpy as np
from array import array
import sys

if sys.argv[1] == "median":
    useMedian = True
else:
    useMedian = False

#fileIn = TFile("histos/xOverX0_2016.root")
#fileIn = TFile("histos/01_07_2021/fBrem_notFolded_PVz_reweighted.root")
#fileIn = TFile("histos/09_12_2021/fBrem_notFolded_postVFP_2016.root")
#fileIn = TFile("histos/09_12_2021/fBrem_notFolded_2016preVFP_PVzReweighted.root")
fileIn = TFile("histos/19_01_2022/highPU/fBrem_notFolded_2017_highPU.root")
#fileIn2 = TFile("histos/19_01_2022/lowPU/fBrem_notFolded_2017_highPU.root")

#fileIn = TFile("histos/01_07_2021/fBrem_notFolded_etaSC.root")
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
y_ratio = array( 'f',[])
y_err_ratio = array( 'f',[])
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
dict_calib_xOverX0_DATA = dict()
dict_calib_xOverX0_MC = dict()
dict_calib_xOverX0_err_DATA = dict()
dict_calib_xOverX0_err_MC = dict()

if useMedian:
    calibration_file = TFile("Calibration_t3.root")
    gr_t3 = calibration_file.Get("Calibration_graph_t3")
else:
    calibration_file = TFile("Calibration_t1.root")


for i in np.arange(-2.50,2.50,0.05):

    str_i = '%.6f' % i
    str_iplus = '%.6f' %(i+0.05)
    #str_i = str_i.replace(".","_")
    #str_iplus = str_iplus.replace(".","_")
    h_dict_DATA[i] = fileIn.Get("DATA/DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
    h_dict_MC[i]   = fileIn.Get("MC/MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
    h_dict_energy_DATA[i] = fileIn.Get("DATA/DATA_h_energy_etaBins_" + str_i + "_" + str_iplus)
    h_dict_energy_MC[i]   = fileIn.Get("MC/MC_h_energy_etaBins_" + str_i + "_" + str_iplus)

    str_i = str_i.replace("0000","")
    str_iplus = str_iplus.replace("0000","")

    #print str_i, str_iplus

    if useMedian:
        median_DATA = np.zeros(1, dtype=float)
        median_MC = np.zeros(1, dtype=float)
        prob = np.zeros(1, dtype=float)
        prob[0] = 0.5
        h_dict_DATA[i].ComputeIntegral()
        h_dict_MC[i].ComputeIntegral()
        h_dict_DATA[i].GetQuantiles(1,median_DATA,prob)
        h_dict_MC[i].GetQuantiles(1,median_MC,prob)
        xOverX0_DATA = -(math.log(1 - median_DATA))
        xOverX0_MC   = -(math.log(1 - median_MC))
        uncert_median_DATA = 1.2533*(h_dict_DATA[i].GetStdDev()/math.sqrt(h_dict_DATA[i].GetEntries())) #1.2533 *std.dev/sqrt(N)
        uncert_median_MC   = 1.2533*(h_dict_MC[i].GetStdDev()/math.sqrt(h_dict_MC[i].GetEntries())) #1.2533 *std.dev/sqrt(N)
        uncert_xOverX0_DATA = uncert_median_DATA/(1. - median_DATA) #error propagation
        uncert_xOverX0_MC   = uncert_median_MC/(1. - median_MC) #error propagation

        calib_xOverX0_DATA = gr_t3.Eval(xOverX0_DATA)
        calib_xOverX0_MC   = gr_t3.Eval(xOverX0_MC)

        y_DATA.append(calib_xOverX0_DATA)
        y_MC.append(calib_xOverX0_MC)
        y_err_DATA.append(uncert_median_DATA)
        y_err_MC.append(uncert_median_MC)

    else:
        mean_DATA = h_dict_DATA[i].GetMean()
        mean_MC   = h_dict_MC[i].GetMean()
        uncert_mean_DATA = h_dict_DATA[i].GetStdDev()/math.sqrt(h_dict_DATA[i].GetEntries()) #std.dev/sqrt(N)
        uncert_mean_MC   = h_dict_MC[i].GetStdDev()/math.sqrt(h_dict_MC[i].GetEntries()) #std.dev/sqrt(N)
        
        xOverX0_DATA = -(math.log(1 - mean_DATA))
        xOverX0_MC   = -(math.log(1 - mean_MC))
        uncert_xOverX0_DATA = uncert_mean_DATA/(1. - mean_DATA) #error propagation
        uncert_xOverX0_MC   = uncert_mean_MC/(1. - mean_MC) #error propagation

        mean_energy_DATA = h_dict_energy_DATA[i].GetMean()
        mean_energy_MC   = h_dict_energy_MC[i].GetMean()

        #Choose the most appropriate calibration file
        if math.fabs(mean_energy_DATA - 25.) < math.fabs(mean_energy_DATA - 50.):
            gr_t1 = calibration_file.Get("Calibration_graph_t1_25_GeV")
        elif math.fabs(mean_energy_DATA - 25.) > math.fabs(mean_energy_DATA - 50.) and math.fabs(mean_energy_DATA - 50.) < math.fabs(mean_energy_DATA - 100.):
            gr_t1 = calibration_file.Get("Calibration_graph_t1_50_GeV")
        elif math.fabs(mean_energy_DATA - 50.) > math.fabs(mean_energy_DATA - 100.) and math.fabs(mean_energy_DATA - 100.) < math.fabs(mean_energy_DATA - 200.):
            gr_t1 = calibration_file.Get("Calibration_graph_t1_100_GeV")
        elif math.fabs(mean_energy_DATA - 100.) > math.fabs(mean_energy_DATA - 200.):
            gr_t1 = calibration_file.Get("Calibration_graph_t1_200_GeV")

        calib_xOverX0_DATA = gr_t1.Eval(xOverX0_DATA)
        calib_xOverX0_MC   = gr_t1.Eval(xOverX0_MC)

        print str_i,"_",str_iplus

        if str_i == "-0.00":
            #Calculation of uncertainties must happen before calculation of the central values
            dict_calib_xOverX0_err_DATA["0.05"] = (2./(dict_calib_xOverX0_DATA["0.05"] + gr_t1.Eval(xOverX0_DATA))**2)*math.sqrt((gr_t1.Eval(xOverX0_DATA)*dict_calib_xOverX0_err_DATA["0.05"])**2 + ( dict_calib_xOverX0_DATA["0.05"]*uncert_xOverX0_DATA)**2)
            dict_calib_xOverX0_err_MC["0.05"] = (2./(dict_calib_xOverX0_MC["0.05"] + gr_t1.Eval(xOverX0_MC))**2)*math.sqrt((gr_t1.Eval(xOverX0_MC)*dict_calib_xOverX0_err_MC["0.05"])**2 + ( dict_calib_xOverX0_MC["0.05"]*uncert_xOverX0_MC)**2)

            dict_calib_xOverX0_DATA["0.05"] = (gr_t1.Eval(xOverX0_DATA) - dict_calib_xOverX0_DATA["0.05"])/(gr_t1.Eval(xOverX0_DATA) + dict_calib_xOverX0_DATA["0.05"])  
            dict_calib_xOverX0_MC["0.05"] = (gr_t1.Eval(xOverX0_MC) - dict_calib_xOverX0_MC["0.05"])/(gr_t1.Eval(xOverX0_MC) + dict_calib_xOverX0_MC["0.05"])

        if not str_i.startswith("-"):

            #Calculation of uncertainties must happen before calculation of the central values
            dict_calib_xOverX0_err_DATA[str_iplus] = (2./(dict_calib_xOverX0_DATA[str_iplus] + gr_t1.Eval(xOverX0_DATA))**2)*math.sqrt((gr_t1.Eval(xOverX0_DATA)*dict_calib_xOverX0_err_DATA[str_iplus])**2 + ( dict_calib_xOverX0_DATA[str_iplus]*uncert_xOverX0_DATA)**2)
            dict_calib_xOverX0_err_MC[str_iplus] = (2./(dict_calib_xOverX0_MC[str_iplus] + gr_t1.Eval(xOverX0_MC))**2)*math.sqrt((gr_t1.Eval(xOverX0_MC)*dict_calib_xOverX0_err_MC[str_iplus])**2 + ( dict_calib_xOverX0_MC[str_iplus]*uncert_xOverX0_MC)**2) 

            dict_calib_xOverX0_DATA[str_iplus] = (gr_t1.Eval(xOverX0_DATA) - dict_calib_xOverX0_DATA[str_iplus])/(gr_t1.Eval(xOverX0_DATA) + dict_calib_xOverX0_DATA[str_iplus])
            dict_calib_xOverX0_MC[str_iplus] = (gr_t1.Eval(xOverX0_MC) - dict_calib_xOverX0_MC[str_iplus])/(gr_t1.Eval(xOverX0_MC) + dict_calib_xOverX0_MC[str_iplus])


        if str_i.startswith("-") and not str_i == "-0.00": #Do not change the order of these if conditions
            str_i = str_i.replace("-","")
            dict_calib_xOverX0_DATA[str_i] = gr_t1.Eval(xOverX0_DATA)
            dict_calib_xOverX0_MC[str_i]   = gr_t1.Eval(xOverX0_MC)
            dict_calib_xOverX0_err_DATA[str_i] = uncert_xOverX0_DATA
            dict_calib_xOverX0_err_MC[str_i]   = uncert_xOverX0_MC

        y_DATA.append(calib_xOverX0_DATA)
        y_MC.append(calib_xOverX0_MC)

        y_err_DATA.append(uncert_xOverX0_DATA)
        y_err_MC.append(uncert_xOverX0_MC)

    x.append(i+0.025)
    x_err.append(0.)

c_graph = TCanvas()
c_graph.SetGrid()

pad1 = TPad("pad_values","",0,0.28,1,1.)
pad2 = TPad("pad_ratio",'',0,0,1,0.25)
pad1.SetBottomMargin(0.01)
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
g_xOverX0_DATA = TGraphErrors(len(x), x, y_DATA, x_err, y_err_DATA)
g_xOverX0_DATA.SetMarkerColor(2)
g_xOverX0_DATA.SetMarkerStyle(20)
g_xOverX0_DATA.SetTitle("")
g_xOverX0_DATA.GetXaxis().SetTitle("#eta^{e}")
g_xOverX0_DATA.GetXaxis().SetLabelSize(0)
if useMedian:
    g_xOverX0_DATA.GetYaxis().SetTitle("-ln [1 - #it{med}(fbrem)]")
else:
    g_xOverX0_DATA.GetYaxis().SetTitle("-ln (1 - <fbrem>)")
g_xOverX0_DATA.GetYaxis().SetTitleSize(0.05)
g_xOverX0_DATA.GetYaxis().SetRangeUser(0.,1.7)
g_xOverX0_DATA.Draw("AP")

g_xOverX0_MC = TGraphErrors(len(x), x, y_MC, x_err, y_err_MC)
g_xOverX0_MC.SetMarkerColor(kGreen+1)
g_xOverX0_MC.SetMarkerStyle(24)
g_xOverX0_MC.SetTitle("")
g_xOverX0_MC.GetXaxis().SetTitle("#eta^{e}")
g_xOverX0_MC.GetXaxis().SetLabelSize(0)
if useMedian:
    g_xOverX0_MC.GetYaxis().SetTitle("-ln [1 - #it{med}(fbrem)]")
else:
    g_xOverX0_MC.GetYaxis().SetTitle("-ln (1 - <fbrem>)")
g_xOverX0_MC.GetYaxis().SetTitleSize(0.05)
g_xOverX0_MC.GetYaxis().SetRangeUser(0.,1.7)
g_xOverX0_MC.Draw("SAME,P")

legend = TLegend(0.7,0.7,0.85,0.95)
legend.SetHeader(" ")
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.SetLineColor(1)
legend.SetLineStyle(1)
legend.SetLineWidth(1)
legend.SetFillStyle(0)
legend.AddEntry(g_xOverX0_DATA,"Data","lep")
legend.AddEntry(g_xOverX0_MC,"MC","lep")
legend.Draw("SAME")

pad2.cd()
pad2.SetTopMargin(0.03)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)

#for a_DATA, a_MC in zip(y_DATA,y_MC):
#    y_diff.append(a_DATA-a_MC)
#    y_err_diff.append(0.)

for a_DATA, a_MC, a_err_DATA, a_err_MC in zip(y_DATA,y_MC,y_err_DATA,y_err_MC):
    y_ratio.append(a_DATA/a_MC)
    y_err_ratio.append(math.sqrt(a_err_DATA**2 / a_MC**2 + a_DATA**2 * a_err_MC**2 / a_MC**4))

g_ratio = TGraphErrors(len(x), x, y_ratio, x_err, y_err_ratio)
g_ratio.SetMarkerStyle(20)
g_ratio.SetMarkerSize(0.5)
g_ratio.SetTitle("")
g_ratio.GetXaxis().SetTitle("#eta^{e}")
g_ratio.GetYaxis().SetTitle("Data/MC")
#g_ratio.GetYaxis().SetRangeUser(-0.2,0.2)
g_ratio.GetYaxis().SetRangeUser(0.9,1.1)
g_ratio.GetXaxis().SetLabelSize(0.15)
g_ratio.GetXaxis().SetTitleSize(0.15)
g_ratio.GetXaxis().SetTitleOffset(0.7)
g_ratio.GetYaxis().SetTitleSize(0.14)
g_ratio.GetYaxis().SetTitleOffset(0.35)
g_ratio.GetYaxis().SetLabelSize(0.10)
g_ratio.GetYaxis().SetNdivisions(502,False)
line_on_one = TLine(g_xOverX0_DATA.GetXaxis().GetXmin(),1.,g_xOverX0_DATA.GetXaxis().GetXmax(),1.)
line_on_one.SetLineColor(4)
line_on_one.SetLineStyle(2)
g_ratio.Draw("AP")
line_on_one.Draw("SAME")

c_ratio = TCanvas()
h_ratio = TH1F("h_ratio","", 18,0.92,1.1)
for it in xrange(13,87): #that is for |eta| ~< 1.8
    x_value = np.zeros(1, dtype=float)
    y_value = np.zeros(1, dtype=float)
    g_ratio.GetPoint(it, x_value, y_value)
    #print x_value, y_value
    h_ratio.Fill(y_value)
h_ratio.Draw()

if useMedian:
    c_graph.SaveAs("plots/xOverX0_median_ratio_CutBased_calibrated.pdf")
    c_ratio.SaveAs("plots/ratio_spread_median_CutBased_calibrated.pdf")
else:
    c_graph.SaveAs("plots/xOverX0_ratio_CutBased_2016preVFP_PVzReweighted.pdf")
    c_ratio.SaveAs("plots/ratio_spread_CutBased_2016preVFP_PVzReweighted.pdf")

    for l in np.arange(0.05,2.55,0.05):
        str_i = '%.6f' % l
        str_i = str_i.replace("0000","")
        str_i = str_i.replace("-","")
        x_fBrem_asymmetry.append(l-0.025)
        x_fBrem_asymmetry_err.append(0.)
        y_fBrem_asymmetry_DATA.append(dict_calib_xOverX0_DATA[str_i])
        y_fBrem_asymmetry_MC.append(dict_calib_xOverX0_MC[str_i])
        y_fBrem_asymmetry_err_DATA.append(dict_calib_xOverX0_err_DATA[str_i])
        y_fBrem_asymmetry_err_MC.append(dict_calib_xOverX0_err_MC[str_i])
        #print str_i, ": ", dict_calib_xOverX0_err_MC[str_i]

        y_diff.append(dict_calib_xOverX0_DATA[str_i] - dict_calib_xOverX0_MC[str_i])
        y_err_diff.append(math.sqrt(dict_calib_xOverX0_err_DATA[str_i]**2 + dict_calib_xOverX0_err_MC[str_i]**2))
        

        #print "dict_fBrem_DATA[i]: ", dict_fBrem_DATA[str_i], "  dict_fBrem_MC[i]: ", dict_fBrem_MC[str_i]
        
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
    #g_fBrem_asymmetry_DATA.GetYaxis().SetTitle("(f_{brem}^{+} - f_{brem}^{-}) / (f_{brem}^{+} + f_{brem}^{-})")
    g_fBrem_asymmetry_DATA.GetYaxis().SetTitle("(t_{1}^{+} - t_{1}^{-}) / (t_{1}^{+} + t_{1}^{-})")
    g_fBrem_asymmetry_DATA.GetYaxis().SetTitleOffset(1.0)
    g_fBrem_asymmetry_DATA.GetYaxis().SetTitleSize(0.05)
    g_fBrem_asymmetry_DATA.GetYaxis().SetRangeUser(-0.1,0.1)
    g_fBrem_asymmetry_DATA.Draw("AP")
    
    g_fBrem_asymmetry_MC = TGraphErrors(len(x_fBrem_asymmetry), x_fBrem_asymmetry, y_fBrem_asymmetry_MC, x_fBrem_asymmetry_err, y_fBrem_asymmetry_err_MC)
    g_fBrem_asymmetry_MC.SetMarkerColor(kGreen+1)
    g_fBrem_asymmetry_MC.SetMarkerStyle(24)
    g_fBrem_asymmetry_MC.SetLineColor(kGreen+1)
    g_fBrem_asymmetry_MC.SetTitle("")
    g_fBrem_asymmetry_MC.GetXaxis().SetTitle("|#eta^{e}|")
    g_fBrem_asymmetry_MC.GetXaxis().SetLabelSize(0)
    #g_fBrem_asymmetry_MC.GetYaxis().SetTitle("(f_{brem}^{+} - f_{brem}^{-}) / (f_{brem}^{+} + f_{brem}^{-})")
    g_fBrem_asymmetry_MC.GetYaxis().SetTitle("(t_{1}^{+} - t_{1}^{-}) / (t_{1}^{+} + t_{1}^{-})")
    g_fBrem_asymmetry_MC.GetYaxis().SetTitleOffset(1.0)
    g_fBrem_asymmetry_MC.GetYaxis().SetTitleSize(0.05)
    g_fBrem_asymmetry_MC.GetYaxis().SetRangeUser(-0.1,0.1)
    g_fBrem_asymmetry_MC.Draw("SAME,P")

    fOut_graph = TFile("asymmetry_MC.root","RECREATE")
    fOut_graph.cd()
    g_fBrem_asymmetry_MC.Write()
    
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

    c_asymmetry.SaveAs("plots/t1_asymmetry_2016preVFP_PVzReweighted.pdf")

raw_input()
