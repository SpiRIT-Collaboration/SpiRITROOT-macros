void drawMultiplicity() {
  gStyle -> SetOptStat(0);
  gStyle -> SetPadRightMargin(0.04);
  gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadBottomMargin(0.15);
  gStyle -> SetPadLeftMargin(0.18);
  gStyle -> SetTitleOffset(1.2, "x");
  gStyle -> SetTitleOffset(1.6, "y");
  gStyle -> SetTitleSize(0.06, "x");
  gStyle -> SetTitleSize(0.06, "y");
  gStyle -> SetLabelSize(0.06, "x");
  gStyle -> SetLabelSize(0.06, "y");
  gStyle -> SetLegendTextSize(0.06);

  Int_t distCut = 5;

  auto chain = new TChain("mult");
  chain -> AddFile(Form("multSn132_%d.root", distCut));

  auto c1 = new TCanvas();

  auto legend = new TLegend(0.70, 0.93, 0.90, 0.54);
  legend -> SetMargin(0.5);
  legend -> SetLineStyle(0);
  legend -> SetLineColor(0);
  legend -> SetLineWidth(0);
  Int_t colors[8] = {1, 2, 8, 6, 9, 28, 3, 4};
  TString legendName[8] = {"No Cut", "1.0#sigma", "1.5#sigma", "2.0#sigma (1)", "12", "12", "124", "1245"};

  TH1D **hist = new TH1D*[8];
  for (auto iHist = 0; iHist < 8; iHist++) {
    hist[iHist] = new TH1D(Form("hist%d", iHist), "", 100, 0, 100);
    hist[iHist] -> SetLineWidth(2);

    chain -> Project(Form("hist%d", iHist), Form("mult%d", iHist), Form("mult%d!=0", iHist));
    if (iHist == 0) {
      hist[iHist] -> GetXaxis() -> SetTitle("Multiplicity");
      hist[iHist] -> GetXaxis() -> CenterTitle();
      hist[iHist] -> GetYaxis() -> SetTitle("Event counts");
      hist[iHist] -> GetYaxis() -> CenterTitle();
      hist[iHist] -> GetYaxis() -> SetRangeUser(0, 6e4);
      hist[iHist] -> SetLineColor(colors[iHist]);
      hist[iHist] -> SetLineColor(colors[0]);
      hist[iHist] -> Draw();
    } else {
      hist[iHist] -> SetLineColor(colors[iHist]);
      hist[iHist] -> Draw("same");
    }

    legend -> AddEntry(hist[iHist], legendName[iHist], "L");
  }

  legend -> Draw();

  auto text = new TLatex();
  text -> SetNDC();
  text -> DrawLatex(0.73, 0.47, "2 - vertex");
  text -> DrawLatex(0.73, 0.42, Form("4 - dist < %d", distCut));
  text -> DrawLatex(0.73, 0.32, "5 - NDF > 30");
}
