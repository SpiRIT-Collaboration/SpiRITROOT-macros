auto mpi = 139.57018;
auto mp = 938.2720813;
auto mn = 939.565346;
auto md = 1875.612762;
auto mt = 2808.921112;
auto mhe3 = 2808.39132;
auto mal = 3727.379378;
auto me = 0.5109989461;

void drawpTVSy() {
  gStyle -> SetOptStat(0);
  gStyle -> SetPadRightMargin(gStyle -> GetPadRightMargin() + 0.05);

  auto cutFile = new TFile("pdFuncCut.root");
  auto pFuncCut = (TCutG *) cutFile -> FindObjectAny("pFuncCut");

  auto tree = new TChain("dedx");
  Bool_t sigma10, sigma15, sigma20;
  Double_t vx, vy, vz;
  Double_t dedx, mom, ndf, dx, dy, dz, dist, pt, rapCM, rapL;
  Int_t charge, pid, run, eventid, vid, parentvid;
  tree -> AddFile("../dedxSn132-LayerCut90.root");

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
  tree -> SetBranchAddress("rapL", &rapL);
  tree -> SetBranchAddress("rapCM", &rapCM);

  auto theta = [&](Bool_t deg = kTRUE) -> Double_t {
    return TMath::ATan2(TMath::Sqrt(dx*dx + dy*dy), dz)*(deg ? 180./TMath::Pi() : 1.);
  };

  auto phi = [&](Bool_t deg = kTRUE) -> Double_t {
    Double_t value = TMath::ATan2(dy, dx)*(deg ? 180./TMath::Pi() : 1.);

    return value + (value < 0 ? 360. : 0);
  };

  auto hist = new TH2D("hist", "", 500, -1.5, 3, 500, 0, 1000);
  auto hist2 = new TH2D("hist2", "", 500, -1.5, 3, 500, 0, 1000);

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

    if (!(pFuncCut -> IsInside(mom/charge, dedx)))
      continue;

    hist -> Fill(rapCM, pt);

    if (!((phi() > 0 && phi() < 40) || (phi() > 150 && phi() < 230) || (phi() > 310 && phi() < 360)))
      continue;

    hist2 -> Fill(rapCM, pt);
  }

  auto c1 = new TCanvas("cvs", "");
  hist -> Draw("colz");

  auto c2 = new TCanvas("cvs2", "");
  hist2 -> Draw("colz");
}
