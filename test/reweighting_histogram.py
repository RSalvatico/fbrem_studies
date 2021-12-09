from ROOT import TFile, TH1F, TCanvas
import sys

isData = True
if sys.argv[1] == "MC":
    isData = False

fIn_DATA_MC = TFile("histos/cpp/plots_CutBased_DYNLO.root")
#fIn_MC = TFile()
fIn_Gauss = TFile("PVz_reweighting.root")

if isData:
    h_den = fIn_DATA_MC.Get("h_PVz_DATA")
else:
    h_den = fIn_DATA_MC.Get("h_PVz_MC")

h_num = fIn_Gauss.Get("h_reweight_PVz")

#canvas_1 = TCanvas()
# h_num.Draw("hist")
# canvas_2 = TCanvas()
# h_den.Draw("hist")
#h_num.Scale(1./h_num.Integral()) 
#h_den.Scale(1./h_den.Integral()) 
h_num.Divide(h_den)
#h_num.Scale(1./h_num.Integral())

canvas = TCanvas()
h_num.GetXaxis().SetTitle("PV z (cm)")
h_num.GetYaxis().SetTitle("Events/0.5")
h_num.Draw("hist")

if isData:
    fOut = TFile("PVz_weights_DATA.root","recreate")
else:
    fOut = TFile("PVz_weights_MC.root","recreate")

fOut.cd()
h_num.Write()
fOut.Close()

raw_input()

