void refine_beam(){
  TCanvas *cvs2 = new TCanvas("cvspid2", "", 0, 0, 800, 580);
  TF2 *sigmaLine = new TF2("sigmaLine", "bigaus");
  sigmaLine -> SetFillStyle(1000);
  sigmaLine -> SetLineWidth(2);
  sigmaLine -> SetRange(2.635, 49, 2.647, 51);
  sigmaLine -> SetParameters(0.253831, 2.64, 0.00127952, 49.9701, 0.204459, -0.0465522);
  sigmaLine -> SetNpx(1000);
  sigmaLine -> SetNpy(1000);

  double contours[3];
  double corr = sigmaLine->Eval(2.64, 49.9701)/ROOT::Math::bigaussian_pdf(0,0,1,1,0);
  contours[0] = corr*ROOT::Math::bigaussian_pdf(2,2,1,1,0); //0.05;
  contours[1] = corr*ROOT::Math::bigaussian_pdf(1.5,1.5,1,1,0);
  contours[2] = corr*ROOT::Math::bigaussian_pdf(1,1,1,1,0);
  sigmaLine -> SetContour(3,contours);

  sigmaLine -> Draw("cont list");
  cvs2 -> Update();

  TObjArray *myContours = (TObjArray *) gROOT -> GetListOfSpecials() -> FindObject("contours");
  auto graphSigma20 = (TGraph *) ((TList *) myContours -> At(0)) -> First();
  auto graphSigma15 = (TGraph *) ((TList *) myContours -> At(1)) -> First();
  auto graphSigma10 = (TGraph *) ((TList *) myContours -> At(2)) -> First();

  auto cutSigma20 = new TCutG("sigma20");
  auto cutSigma15 = new TCutG("sigma15");
  auto cutSigma10 = new TCutG("sigma10");
  Double_t x = 0, y = 0;
  for (Int_t iPoint = 0; iPoint < graphSigma20 -> GetN(); iPoint++) {
    graphSigma20 -> GetPoint(iPoint, x, y);
    cutSigma20 -> SetPoint(iPoint, x, y);
  }
  for (Int_t iPoint = 0; iPoint < graphSigma15 -> GetN(); iPoint++) {
    graphSigma15 -> GetPoint(iPoint, x, y);
    cutSigma15 -> SetPoint(iPoint, x, y);
  }
  for (Int_t iPoint = 0; iPoint < graphSigma10 -> GetN(); iPoint++) {
    graphSigma10 -> GetPoint(iPoint, x, y);
    cutSigma10 -> SetPoint(iPoint, x, y);
  }
  cutSigma20 -> SetVarX("aoq");
  cutSigma20 -> SetVarY("z");
  cutSigma15 -> SetVarX("aoq");
  cutSigma15 -> SetVarY("z");
  cutSigma10 -> SetVarX("aoq");
  cutSigma10 -> SetVarY("z");


  Int_t firstRun = 2843;
  Int_t lastRun = 2894;
  //add runs to the chain
  for (Int_t iRun = firstRun; iRun <= lastRun; iRun++){
    if (!(gSystem -> IsFileInIncludePath(Form("output/beam/beam_run%i.ridf.root", iRun))))
      continue;

    auto beam = new TChain("TBeam");
    auto bdc = new TChain("TBDC");
    Double_t aoq = 0, z = 0, projmevu, proja, projb, projpx, projpy, projpz, projx, projy, projz;
    beam -> AddFile(Form("output/beam/beam_run%i.ridf.root", iRun));
    beam -> SetBranchAddress("aoq", &aoq);
    beam -> SetBranchAddress("z", &z);
    bdc -> AddFile(Form("output/beam/beam_run%i.ridf.root", iRun));
    bdc -> SetBranchAddress("ProjA", &proja);
    bdc -> SetBranchAddress("ProjB", &projb);
    bdc -> SetBranchAddress("ProjMeVu", &projmevu);
    bdc -> SetBranchAddress("ProjPX", &projpx);
    bdc -> SetBranchAddress("ProjPY", &projpy);
    bdc -> SetBranchAddress("ProjPZ", &projpz);
    bdc -> SetBranchAddress("ProjX", &projx);
    bdc -> SetBranchAddress("ProjY", &projy);
    bdc -> SetBranchAddress("ProjZ", &projz);

    auto numEntries = beam -> GetEntries();
    cout << " Entries = " << numEntries << endl;

    auto file = new TFile(Form("refined/run%i.refined.root", iRun), "recreate");
    Bool_t sigma10 = kFALSE, sigma15 = kFALSE, sigma20 = kFALSE;
    auto tree = new TTree("beam", "");
    tree -> Branch("aoq", &aoq);
    tree -> Branch("z", &z);
    tree -> Branch("sigma10", &sigma10);
    tree -> Branch("sigma15", &sigma15);
    tree -> Branch("sigma20", &sigma20);
    tree -> Branch("projmevu", &projmevu);
    tree -> Branch("projpx", &projpx);
    tree -> Branch("projpy", &projpy);
    tree -> Branch("projpz", &projpz);
    tree -> Branch("proja", &proja);
    tree -> Branch("projb", &projb);
    tree -> Branch("projx", &projx);
    tree -> Branch("projy", &projy);
    tree -> Branch("projz", &projz);

    for (Int_t iEntry = 0; iEntry < numEntries; iEntry++) {
      if(iEntry%10000 == 0)
        cout << "File read: " << iEntry*100./numEntries << "%" << endl;

      beam -> GetEvent(iEntry);
      bdc -> GetEntry(iEntry);

      sigma20 = graphSigma20 -> IsInside(aoq, z);
      sigma15 = graphSigma15 -> IsInside(aoq, z);
      sigma10 = graphSigma10 -> IsInside(aoq, z);

      tree -> Fill();
    }

    delete beam;
    delete bdc;

    tree -> Write();
  }
}
