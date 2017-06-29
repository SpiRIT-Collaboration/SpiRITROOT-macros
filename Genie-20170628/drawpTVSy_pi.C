auto mpi = 139.57018;
auto mp = 938.2720813;
auto mn = 939.565346;
auto md = 1875.612762;
auto mt = 2808.921112;
auto mhe3 = 2808.39132;
auto mal = 3727.379378;
auto me = 0.5109989461;

void drawpTVSy_pi() {
  gStyle -> SetOptStat(0);
  gStyle -> SetPadRightMargin(gStyle -> GetPadRightMargin());

  auto cutFile = new TFile("pionCut.root");
  auto pim = (TCutG *) cutFile -> FindObjectAny("pim");
  auto pip = (TCutG *) cutFile -> FindObjectAny("pip");

  auto tree = new TChain("dedx");
  Bool_t sigma10, sigma15, sigma20;
  Double_t vx, vy, vz;
  Double_t dedx, mom, ndf, dx, dy, dz, dist, pt, rapCM, rapL, pzCM, KECM;
  Int_t charge, pid, run, eventid, vid, parentvid;
  tree -> AddFile("dedxSn132-LayerCut90.root");

  tree -> SetBranchAddress("run", &run);
  tree -> SetBranchAddress("eventid", &eventid);
  tree -> SetBranchAddress("dedx", &dedx);
  tree -> SetBranchAddress("mom", &mom);
  tree -> SetBranchAddress("charge", &charge);
  tree -> SetBranchAddress("pid", &pid);
  tree -> SetBranchAddress("ndf", &ndf);
  tree -> SetBranchAddress("dx", &dx); // mom direction
  tree -> SetBranchAddress("dy", &dy);
  tree -> SetBranchAddress("dz", &dz);
  tree -> SetBranchAddress("vid", &vid); // vertex id 
  tree -> SetBranchAddress("parentvid", &parentvid); // track's parent vertex id 
  tree -> SetBranchAddress("vx", &vx); // vertex position
  tree -> SetBranchAddress("vy", &vy);
  tree -> SetBranchAddress("vz", &vz);
  tree -> SetBranchAddress("dist", &dist); // poca and vertex distance
  tree -> SetBranchAddress("sigma10", &sigma10);
  tree -> SetBranchAddress("sigma15", &sigma15);
  tree -> SetBranchAddress("sigma20", &sigma20);
  tree -> SetBranchAddress("pt", &pt);
  tree -> SetBranchAddress("pzCM", &pzCM);
  tree -> SetBranchAddress("KECM", &KECM);
  tree -> SetBranchAddress("rapL", &rapL);
  tree -> SetBranchAddress("rapCM", &rapCM);

  auto theta = [&](Bool_t deg = kTRUE) -> Double_t {
    return TMath::ATan2(TMath::Sqrt(dx*dx + dy*dy), dz)*(deg ? 180./TMath::Pi() : 1.);
  };

  auto phi = [&](Bool_t deg = kTRUE) -> Double_t {
    Double_t value = TMath::ATan2(dy, dx)*(deg ? 180./TMath::Pi() : 1.);

    return value + (value < 0 ? 360. : 0);
  };

  auto hpim = new TH1D("hpim", "p - pim; p (MeV/c);", 100, 0, 500);
  auto hpip = new TH1D("hpip", "p - pip; p (MeV/c);", 100, 0, 500);
  hpip -> SetLineColor(2);
  auto hKEpim = new TH1D("hKEpim", "KE - pim; E_{CM} (MeV);", 100, 0, 310);
  auto hKEpip = new TH1D("hKEpip", "KE - pip; E_{CM} (MeV);", 100, 0, 310);
  hKEpip -> SetLineColor(2);
  auto hpim2 = new TH2D("hpim2", "pT vs y - pim; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);
  auto hpip2 = new TH2D("hpip2", "pT vs y - pip; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);
  auto hpi2 = new TH2D("hpi2", "pT vs y - pip; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);

  auto hpimAngle = new TH1D("hpimAngle", "p - pim - #phi cut; p (MeV/c)", 100, 0, 500);
  auto hpipAngle = new TH1D("hpipAngle", "p - pip - #phi cut; p (MeV/c)", 100, 0, 500);
  hpipAngle -> SetLineColor(2);
  auto hKEpimAngle = new TH1D("hKEpimAngle", "KE - pim - #phi cut; E_{CM} (MeV);", 60, 0, 310);
  auto hKEpipAngle = new TH1D("hKEpipAngle", "KE - pip - #phi cut; E_{CM} (MeV);", 60, 0, 310);
  hKEpipAngle -> SetLineColor(2);
  auto hpim2Angle = new TH2D("hpim2Angle", "pT vs y - pim - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);
  auto hpip2Angle = new TH2D("hpip2Angle", "pT vs y - pip - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);
  auto hpi2Angle = new TH2D("hpi2Angle", "pT vs y - pi - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 30, -1.5, 4, 30, 0, 450);

  auto legend = new TLegend(0.6, 0.6, 0.8, 0.8);
  legend -> SetLineColor(0);
  legend -> AddEntry(hpim, "#pi^{-}", "L");
  legend -> AddEntry(hpip, "#pi^{+}", "L");

  auto numEntries = tree -> GetEntries();
  for (Int_t iEntry = 0; iEntry < numEntries; iEntry++) {
    if (iEntry%1000000 == 0)
      cout << "Entry: " << iEntry << "/" << numEntries << " (" << (Double_t)iEntry/numEntries*100. << "%)" << endl;

    tree -> GetEntry(iEntry);

    if (!sigma20)
      continue;

    if (!(vz < -9.49569 && vz > -12.80121))
      continue;

    if (!(vx > -15 && vx < 15 && vy < -206.06 && vy > -246.06))
      continue;

    if (dist > 5)
      continue;

    if (ndf < 30)
      continue;

    if (pim -> IsInside(mom/charge, dedx)) {
      hpim -> Fill(abs(mom));
      hKEpim -> Fill(KECM);
      hpim2 -> Fill(rapCM, pt);
      hpi2 -> Fill(rapCM, pt);

      if (!((phi() > 0 && phi() < 40) || (phi() > 150 && phi() < 230) || (phi() > 310 && phi() < 360)))
        continue;

      hpimAngle -> Fill(abs(mom));
      hKEpimAngle -> Fill(KECM);
      hpim2Angle -> Fill(rapCM, pt);
      hpi2Angle -> Fill(rapCM, pt);
    } else if (pip -> IsInside(mom/charge, dedx)) {
      hpip -> Fill(abs(mom));
      hKEpip -> Fill(KECM);
      hpip2 -> Fill(rapCM, pt);
      hpi2 -> Fill(rapCM, pt);

      if (!((phi() > 0 && phi() < 40) || (phi() > 150 && phi() < 230) || (phi() > 310 && phi() < 360)))
        continue;

      hpipAngle -> Fill(abs(mom));
      hKEpipAngle -> Fill(KECM);
      hpip2Angle -> Fill(rapCM, pt);
      hpi2Angle -> Fill(rapCM, pt);
    } 
  }

  auto c1 = new TCanvas();
  hpim -> Draw();
  hpip -> Draw("same");

  legend -> Draw();

  auto c2 = new TCanvas();
  hpimAngle -> Draw();
  hpipAngle -> Draw("same");

  legend -> Draw();

  auto c3 = new TCanvas();
  hKEpim -> Draw();
  hKEpip -> Draw("same");

  legend -> Draw();

  auto c4 = new TCanvas();
  hKEpimAngle -> Draw();
  hKEpipAngle -> Draw("same");

  legend -> Draw();

  auto c5 = new TCanvas();
  hpim2 -> Draw("colz");

  auto c6 = new TCanvas();
  hpip2 -> Draw("colz");

  auto c7 = new TCanvas();
  hpim2Angle -> Draw("colz");

  auto c8 = new TCanvas();
  hpip2Angle -> Draw("colz");

  auto c9 = new TCanvas();
  auto hratio = (TH1D *) hKEpim -> Clone("hratio");
  hratio -> SetTitle("KE ratio");
  hratio -> GetYaxis() -> SetTitle("#pi^{-}/#pi^{+}");
  hratio -> Sumw2();
  hratio -> Rebin(8);
  auto temp = (TH1D *) hKEpip -> Clone("hKEpipTemp");
  temp -> Sumw2();
  temp -> Rebin(8);
  hratio -> Divide(temp);
  hratio -> Draw();

  auto c10 = new TCanvas();
  auto hratioAngle = (TH1D *) hKEpimAngle -> Clone("hratioAngle");
  hratioAngle -> SetTitle("KE ratio - #phi cut");
  hratioAngle -> GetYaxis() -> SetTitle("#pi^{-}/#pi^{+}");
  hratioAngle -> Sumw2();
  hratioAngle -> Rebin(8);
  temp = (TH1D *) hKEpipAngle -> Clone("hKEpipAngleTemp");
  temp -> Sumw2();
  temp -> Rebin(8);
  hratioAngle -> Divide(temp);
  hratioAngle -> Draw();

  auto c11 = new TCanvas();
  auto hratioCor = (TH1D *) hKEpim -> Clone("hratioCor");
  hratioCor -> SetTitle("KE ratio (BG subtracted)");
  hratioCor -> GetYaxis() -> SetTitle("#pi^{-}/#pi^{+}");
  hratioCor -> Sumw2();
  hratioCor -> Rebin(4);
  temp = (TH1D *) hKEpipCor -> Clone("hKEpipCorTemp");
  temp -> Sumw2();
  temp -> Rebin(4);
  hratioCor -> Divide(temp);
  hratioCor -> Draw();

  auto c12 = new TCanvas();
  auto hratioAngleCor = (TH1D *) hKEpimAngle -> Clone("hKEratioAngleCor");
  hratioAngleCor -> SetTitle("KE ratio - #phi cut (BG subtracted)");
  hratioAngleCor -> GetYaxis() -> SetTitle("#pi^{-}/#pi^{+}");
  hratioAngleCor -> Sumw2();
  hratioAngleCor -> Rebin(4);
  temp = (TH1D *) hKEpipAngleCor -> Clone("hKEpipAngleCorTemp");
  temp -> Sumw2();
  temp -> Rebin(4);
  hratioAngleCor -> Divide(temp);
  hratioAngleCor -> Draw();

  auto c13 = new TCanvas();
  hpi2 -> Draw("colz");

  auto c14 = new TCanvas();
  hpi2Angle -> Draw("colz");

/*
  c1 -> SaveAs("pDist.png");
  c2 -> SaveAs("pDistAngle.png");
  c3 -> SaveAs("KEDist.png");
  c4 -> SaveAs("KEDistAngle.png");
  c5 -> SaveAs("pTVSy_pim.png");
  c6 -> SaveAs("pTVSy_pip.png");
  c7 -> SaveAs("pTVSy_pim_Angle.png");
  c8 -> SaveAs("pTVSy_pip_Angle.png");
  c9 -> SaveAs("KEratio_pi.png");
  c10 -> SaveAs("KEratioAngle_pi.png");
  c11 -> SaveAs("KEratioCor_pi.png");
  c12 -> SaveAs("KEratioAngleCor_pi.png");
  c13 -> SaveAs("pTVSy_pi.png");
  c14 -> SaveAs("pTVSy_pi_Angle.png");
  */
}
