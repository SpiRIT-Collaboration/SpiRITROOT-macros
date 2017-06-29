auto mpi = 139.57018;
auto mp = 938.2720813;
auto mn = 939.565346;
auto md = 1875.612762;
auto mt = 2808.921112;
auto mhe3 = 2808.39132;
auto mal = 3727.379378;
auto me = 0.5109989461;

void drawpTVSy_proton_multCut() {
//  gStyle -> SetOptStat(0);
  gStyle -> SetPadRightMargin(gStyle -> GetPadRightMargin());

  Int_t multrun, multeventid, mult4;
  auto multFile = new TChain("mult");
  multFile -> AddFile("multSn132_5.root");
  multFile -> SetBranchAddress("run", &multrun);
  multFile -> SetBranchAddress("eventid", &multeventid);
  multFile -> SetBranchAddress("mult4", &mult4);

  auto cutFile = new TFile("pdFuncCut.root");
  auto pFuncCut = (TCutG *) cutFile -> FindObjectAny("pFuncCut");

  auto tree = new TChain("dedx");
  tree -> AddFile("../dedxSn132-LayerCut90.root");
  Bool_t sigma10, sigma15, sigma20;
  Double_t vx, vy, vz;
  Double_t dedx, mom, ndf, dx, dy, dz, dist, pt, rapCM, rapL, pzCM, KECM;
  Int_t charge, pid, run, eventid, vid, parentvid;

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

/*
  auto hpim = new TH1D("hpim", "p - pim; p (MeV/c);", 100, 0, 500);
  auto hpip = new TH1D("hpip", "p - pip; p (MeV/c);", 100, 0, 500);
  auto hpipBg = new TH1D("hpipBg", "p - pip; p (MeV/c);", 100, 0, 500);
  hpip -> SetLineColor(2);
  hpipBg -> SetLineColor(8);
  auto hKEpim = new TH1D("hKEpim", "KE - pim; E_{CM} (MeV);", 100, 0, 310);
  auto hKEpip = new TH1D("hKEpip", "KE - pip; E_{CM} (MeV);", 100, 0, 310);
  auto hKEpipBg = new TH1D("hKEpipBg", "KE - pip; E_{CM} (MeV);", 100, 0, 310);
  hKEpip -> SetLineColor(2);
  hKEpipBg -> SetLineColor(8);
  */
  auto hproton2 = new TH2D("hproton2", "pT vs y; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);

  auto hproton3 = new TH2D("hproton3", "pT vs y; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);
  auto hproton4 = new TH2D("hproton4", "pT vs y; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);

/*
  auto hpimAngle = new TH1D("hpimAngle", "p - pim - #phi cut; p (MeV/c)", 100, 0, 500);
  auto hpipAngle = new TH1D("hpipAngle", "p - pip - #phi cut; p (MeV/c)", 100, 0, 500);
  auto hpipAngleBg = new TH1D("hpipAngleBg", "p - pip - #phi cut; p (MeV/c)", 100, 0, 500);
  hpipAngle -> SetLineColor(2);
  hpipAngleBg -> SetLineColor(8);
  auto hKEpimAngle = new TH1D("hKEpimAngle", "KE - pim - #phi cut; E_{CM} (MeV);", 100, 0, 310);
  auto hKEpipAngle = new TH1D("hKEpipAngle", "KE - pip - #phi cut; E_{CM} (MeV);", 100, 0, 310);
  auto hKEpipAngleBg = new TH1D("hKEpipAngleBg", "KE - pip - #phi cut; E_{CM} (MeV);", 100, 0, 310);
  hKEpipAngle -> SetLineColor(2);
  hKEpipAngleBg -> SetLineColor(8);
  */
  auto hproton2Angle = new TH2D("hproton2Angle", "pT vs y - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);

  auto hproton3Angle = new TH2D("hproton3Angle", "pT vs y - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);
  auto hproton4Angle = new TH2D("hproton4Angle", "pT vs y - #phi cut; y_{CM}/y_{Beam,CM};p_{T} (MeV/c)", 70, -1.5, 3, 70, 0, 1000);

/*
  auto legend = new TLegend(0.6, 0.6, 0.8, 0.8);
  legend -> SetLineColor(0);
  legend -> AddEntry(hpim, "#pi^{-}", "L");
  legend -> AddEntry(hpip, "#pi^{+}", "L");
//  legend -> AddEntry(hpipBg, "#pi^{+} BG", "L");

  auto legend2 = new TLegend(0.6, 0.6, 0.8, 0.8);
  egend2 -> SetLineColor(0);
  legend2 -> AddEntry(hpim, "#pi^{-}", "L");
  legend2 -> AddEntry(hpip, "#pi^{+}-BG", "L");
  */

  auto lowCut = 35;
  auto highCut = 55;
  auto numMult = multFile -> GetEntries();
  auto numEntries = tree -> GetEntries();
  auto countl = 0;
  auto countm = 0;
  auto counth = 0;
  auto iEntry = 0;
  for (Int_t iMult = 0; iMult < numMult; iMult++) {
    multFile -> GetEntry(iMult);

    if (mult4 == 0)
      continue;

    if (mult4 <= lowCut)
      countl++;
    else if (mult4 > lowCut && mult4 < highCut)
      countm++;
    else if (mult4 >= highCut)
      counth++;

    for (; iEntry < numEntries; iEntry++) {
      if (iEntry%1000000 == 0)
        cout << "Entry: " << iEntry << "/" << numEntries << " (" << (Double_t)iEntry/numEntries*100. << "%)" << endl;

      tree -> GetEntry(iEntry);

      if (run < multrun)
        continue;
      else if (run > multrun)
        break;

      if (eventid < multeventid)
        continue;
      else if (eventid > multeventid)
        break;

      if (run != multrun || eventid != multeventid)
        cout << "run:" << run << " multrun:" << multrun << " eventid:" << eventid << " multeventid:" << multeventid << endl;

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

      if (pFuncCut -> IsInside(mom/charge, dedx)) {
        if (mult4 <= lowCut) {
          hproton2 -> Fill(rapCM, pt);
        } else if (mult4 > lowCut && mult4 < highCut) {
          hproton3 -> Fill(rapCM, pt);
        } else if (mult4 >= highCut) {
          hproton4 -> Fill(rapCM, pt);
        }

        if (!((phi() > 0 && phi() < 40) || (phi() > 150 && phi() < 230) || (phi() > 310 && phi() < 360)))
          continue;

        if (mult4 < lowCut) {
          hproton2Angle -> Fill(rapCM, pt);
        } else if (mult4 > lowCut && mult4 < highCut) {
          hproton3Angle -> Fill(rapCM, pt);
        } else if (mult4 >= highCut) {
          hproton4Angle -> Fill(rapCM, pt);
        }
      } 
    }
  }

/*
  auto c1 = new TCanvas();
  hpim -> Draw();
  hpip -> Draw("same");
  hpipBg -> Draw("same");

  legend -> Draw();

  auto c1_1 = new TCanvas();
  hpim -> Draw();
  auto hpipCor = (TH1D *) hpip -> Clone("hpipCor");
  hpipCor -> Add(hpipBg, -1);
  hpipCor -> Draw("same");

  legend2 -> Draw();

  auto c2 = new TCanvas();
  hpimAngle -> Draw();
  hpipAngle -> Draw("same");
  hpipAngleBg -> Draw("same");

  legend -> Draw();

  auto c2_1 = new TCanvas();
  hpimAngle -> Draw();
  auto hpipAngleCor = (TH1D *) hpipAngle -> Clone("hpipAngleCor");
  hpipAngleCor -> Add(hpipAngleBg, -1);
  hpipAngleCor -> Draw("same");
  
  legend2 -> Draw();

  auto c3 = new TCanvas();
  hKEpim -> Draw();
  hKEpip -> Draw("same");
  hKEpipBg -> Draw("same");

  legend -> Draw();

  auto c3_1 = new TCanvas();
  hKEpim -> Draw();
  auto hKEpipCor = (TH1D *) hKEpip -> Clone("hKEpipCor");
  hKEpipCor -> Add(hKEpipBg, -1);
  hKEpipCor -> Draw("same");

  legend2 -> Draw();

  auto c4 = new TCanvas();
  hKEpimAngle -> Draw();
  hKEpipAngle -> Draw("same");
//  hKEpipAngleBg -> Draw("same");

  legend -> Draw();

  auto c4_1 = new TCanvas();
  hKEpimAngle -> Draw();
  auto hKEpipAngleCor = (TH1D *) hKEpipAngle -> Clone("hKEpipAngleCor");
  hKEpipAngleCor -> Add(hKEpipAngleBg, -1);
  hKEpipAngleCor -> Draw("same");

  legend2 -> Draw();
  */

/*
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
  */

  auto c13 = new TCanvas();
  hproton2 -> Draw("colz");

  auto c14 = new TCanvas();
  hproton2Angle -> Draw("colz");

  auto c15 = new TCanvas();
  hproton3 -> Draw("colz");

  auto c16 = new TCanvas();
  hproton3Angle -> Draw("colz");

  auto c17 = new TCanvas();
  hproton4 -> Draw("colz");

  auto c18 = new TCanvas();
  hproton4Angle -> Draw("colz");

  cout << "<=" << lowCut << ":" << countl++ << endl;
  cout << "<>:" << countm++ << endl;
  cout << ">=" << highCut << ":" << counth++ << endl;

/*
  c1 -> SaveAs("pDist.png");
  c1_1 -> SaveAs("pDistCor.png");
  c2 -> SaveAs("pDistAngle.png");
  c2_1 -> SaveAs("pDistAngleCor.png");
  c3 -> SaveAs("KEDist.png");
  c3_1 -> SaveAs("KEDistCor.png");
  c4 -> SaveAs("KEDistAngle.png");
  c4_1 -> SaveAs("KEDistAngleCor.png");
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
