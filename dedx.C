void dedx()
{
  auto tree = new TChain("cbmsim");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s0.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s1.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s2.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s3.reco.root");

  TClonesArray *trackArray = nullptr;
  tree -> SetBranchAddress("STRecoTrack", &trackArray);

  auto hist = new TH1D("hist","",100,0,100);
  auto hist2 = new TH1D("hist2","",200,0,2000);

  Int_t nEvents = tree -> GetEntries();
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (iEvent % 500 == 0)
      cout << "Event " << iEvent << endl;

    tree -> GetEntry(iEvent);

    Int_t nTracks = trackArray -> GetEntries();
    for (auto iTrack = 0; iTrack < nTracks; ++iTrack) {
      auto track = (STRecoTrack *) trackArray -> At(iTrack);
      auto array = track -> GetdEdxPointArray();
      hist -> Fill(array -> size());
      auto l = 0.;
      for (auto dedx : *array)
        l += dedx.fdx;
      hist2 -> Fill(l);
    }
  }

  new TCanvas();
  hist -> Draw();
  new TCanvas();
  hist2 -> Draw();
}
