void makeMultiplicity() {
  auto dedx = new TChain("dedx");
  dedx -> AddFile("../dedxSn132-LayerCut90.root");

  Int_t eventid = 0;
  Int_t run = 0;
  Int_t vid, parentvid;
  Bool_t sigma10 = 0, sigma15 = 0; sigma20 = 0;
  Double_t vx, vy, vz, dist, ndf, projx, projy, projz;
  dedx -> SetBranchAddress("run", &run);
  dedx -> SetBranchAddress("eventid", &eventid);
  dedx -> SetBranchAddress("vid", &vid);
  dedx -> SetBranchAddress("parentvid", &parentvid);
  dedx -> SetBranchAddress("vx", &vx);
  dedx -> SetBranchAddress("vy", &vy);
  dedx -> SetBranchAddress("vz", &vz);
  dedx -> SetBranchAddress("dist", &dist);
  dedx -> SetBranchAddress("ndf", &ndf);
  dedx -> SetBranchAddress("sigma10", &sigma10);
  dedx -> SetBranchAddress("sigma15", &sigma15);
  dedx -> SetBranchAddress("sigma20", &sigma20);
  dedx -> SetBranchAddress("projx", &projx);
  dedx -> SetBranchAddress("projy", &projy);
  dedx -> SetBranchAddress("projz", &projz);

  Bool_t oldsigma10, oldsigma15, oldsigma20;
  Double_t oldvx, oldvy, oldvz;
  Int_t oldRun = 0;
  Int_t oldEventid = 0;
  Int_t multiplicity[8] = {0};
  auto writeFile = new TFile("multSn132_5.root", "recreate");
  auto writeTree = new TTree("mult", "");
  writeTree -> Branch("run", &oldRun);
  writeTree -> Branch("eventid", &oldEventid);
  writeTree -> Branch("vx", &oldvx);
  writeTree -> Branch("vy", &oldvy);
  writeTree -> Branch("vz", &oldvz);
  writeTree -> Branch("sigma10", &oldsigma10);
  writeTree -> Branch("sigma15", &oldsigma15);
  writeTree -> Branch("sigma20", &oldsigma20);
  writeTree -> Branch("projx", &projx);
  writeTree -> Branch("projy", &projy);
  writeTree -> Branch("projz", &projz);
  writeTree -> Branch("mult0", multiplicity); // no cut
  writeTree -> Branch("mult1", multiplicity + 1); // sigma10
  writeTree -> Branch("mult2", multiplicity + 2); // sigma15
  writeTree -> Branch("mult3", multiplicity + 3); // sigma20
  writeTree -> Branch("mult4", multiplicity + 4); // sigma20 + vertex
  writeTree -> Branch("mult5", multiplicity + 5); // sigma20 + vertex + vid
  writeTree -> Branch("mult6", multiplicity + 6); // sigma20 + vertex + vid + dist
  writeTree -> Branch("mult7", multiplicity + 7); // sigma20 + vertex + vid + dist + ndf

  dedx -> GetEntry(0);
  oldsigma10 = sigma10;
  oldsigma15 = sigma15;
  oldsigma20 = sigma20;
  oldvx = vx;
  oldvy = vy;
  oldvz = vz;
  oldRun = run;
  oldEventid = eventid;

  auto numEntries = dedx -> GetEntries();
  for (auto iEntry = 0; iEntry < numEntries; iEntry++) {
    if (iEntry%100000 == 0)
      cout << "Entry: " << iEntry << "/" << numEntries << " (" << (Double_t)iEntry/numEntries*100. << "%)" << endl;
    dedx -> GetEntry(iEntry);

    if (oldEventid == eventid) {
      multiplicity[0]++;
      if (sigma10)
        multiplicity[1]++;
      if (sigma15)
        multiplicity[2]++;
      if (sigma20) {
        multiplicity[3]++;

        if (!(vz < -9.49569 && vz > -12.80121))
          continue;

        if (!(vx > -15 && vx < 15 && vy < -206.06 && vy > -246.06))
          continue;

        multiplicity[4]++;

        if (!(vid == parentvid))
          continue;

        multiplicity[5]++;

        if (dist > 5.)
          continue;

        multiplicity[6]++;

        if (ndf < 30)
          continue;

        multiplicity[7]++;
      }
    }
    else {
//      cout << "Run: " << oldRun << " Event: " << oldEventid << " Multiplicity: " << multiplicity << endl;

      writeTree -> Fill();

      oldsigma10 = sigma10;
      oldsigma15 = sigma15;
      oldsigma20 = sigma20;
      oldvx = vx;
      oldvy = vy;
      oldvz = vz;
      oldRun = run;
      oldEventid = eventid;
      memset(multiplicity, 0, sizeof(Int_t)*8);
      iEntry--;
    }
  }

  writeTree -> Write();
}
