from ROOT import TFile, TH1D, TH2D, TCanvas, TGraphErrors, TLegend, TPad, TLine, kGreen
import math
import numpy as np
from array import array
import sys

if sys.argv[1] == "median":
    useMedian = True
else:
    useMedian = False

#fileIn = TFile("histos/xOverX0_2016.root")
#fileIn = TFile("histos/fBrem_notFolded_2016.root")
fileIn = TFile("histos/cpp/fBrem_notFolded_full2017_highPU.root")
h_dict_DATA = dict()
h_dict_MC = dict()

#Some arrays for TGraphErrors                                                                                                                         
x = array( 'f',[])
x_err = array( 'f',[])
y_DATA = array( 'f',[])
y_err_DATA = array( 'f',[])
y_MC = array( 'f',[])
y_err_MC = array( 'f',[])
y_diff = array( 'f',[])
y_err_diff = array( 'f',[])

for i in np.arange(-2.50,2.50,0.05):
    str_i = '%.6f' % i
    str_iplus = '%.6f' %(i+0.05)
    #str_i = str_i.replace(".","_")
    #str_iplus = str_iplus.replace(".","_")
    h_dict_DATA[i] = fileIn.Get("DATA/DATA_h_fBrem_etaBins_" + str_i + "_" + str_iplus)
    h_dict_MC[i]   = fileIn.Get("MC/MC_h_fBrem_etaBins_" + str_i + "_" + str_iplus)

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
        y_DATA.append(xOverX0_DATA)
        y_MC.append(xOverX0_MC)
        uncert_median_DATA = 1.2533*(h_dict_DATA[i].GetStdDev()/math.sqrt(h_dict_DATA[i].GetEntries())) #1.2533 *std.dev/sqrt(N)
        uncert_median_MC   = 1.2533*(h_dict_MC[i].GetStdDev()/math.sqrt(h_dict_MC[i].GetEntries())) #1.2533 *std.dev/sqrt(N)
        uncert_xOverX0_DATA = uncert_median_DATA/(1. - median_DATA) #error propagation
        uncert_xOverX0_MC   = uncert_median_MC/(1. - median_MC) #error propagation
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
        y_DATA.append(xOverX0_DATA)
        y_MC.append(xOverX0_MC)

        y_err_DATA.append(uncert_xOverX0_DATA)
        y_err_MC.append(uncert_xOverX0_MC)

    x.append(i+0.025)
    x_err.append(0.)

c_graph = TCanvas()
c_graph.SetGrid()

pad1 = TPad("pad_values","",0,0.28,1,1.)
pad2 = TPad("pad_diff",'',0,0,1,0.25)
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
g_xOverX0_MC.SetMarkerStyle(20)
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

for a_DATA, a_MC in zip(y_DATA,y_MC):
    y_diff.append(a_DATA/a_MC)
    y_err_diff.append(0.)

g_diff = TGraphErrors(len(x), x, y_diff, x_err, y_err_diff)
g_diff.SetMarkerStyle(20)
g_diff.SetMarkerSize(0.5)
g_diff.SetTitle("")
g_diff.GetXaxis().SetTitle("#eta^{e}")
g_diff.GetYaxis().SetTitle("Data/MC")
#g_diff.GetYaxis().SetRangeUser(-0.2,0.2)
g_diff.GetYaxis().SetRangeUser(0.7,1.3)
g_diff.GetXaxis().SetLabelSize(0.15)
g_diff.GetXaxis().SetTitleSize(0.15)
g_diff.GetXaxis().SetTitleOffset(0.7)
g_diff.GetYaxis().SetTitleSize(0.14)
g_diff.GetYaxis().SetTitleOffset(0.35)
g_diff.GetYaxis().SetLabelSize(0.10)
g_diff.GetYaxis().SetNdivisions(502,False)
line_on_zero = TLine(g_xOverX0_DATA.GetXaxis().GetXmin(),1.,g_xOverX0_DATA.GetXaxis().GetXmax(),1.)
line_on_zero.SetLineColor(4)
line_on_zero.SetLineStyle(2)
g_diff.Draw("AP")
line_on_zero.Draw("SAME")

if useMedian:
    c_graph.SaveAs("plots/xOverX0_median_ratio_highPU.pdf")
else:
    c_graph.SaveAs("plots/xOverX0_ratio_highPU.pdf")

raw_input()
