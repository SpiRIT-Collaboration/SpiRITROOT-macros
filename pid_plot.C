void pid_plot
(
  bool drawPIDLine = true,
  STPID::PID selectPID = STPID::kNon, // kPion, kProton, kDeuteron, kTriton, k3He, k4He
  Int_t nbins = 200,
  Int_t p1 = -500,
  Int_t p2 = 2500,
  Int_t dedx1 = 0,
  Int_t dedx2 = 800
)
{
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.10);
  gStyle -> SetPadLeftMargin(0.11);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetTitleFontSize(0.06);

  auto tree = new TChain("cbmsim");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s0.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s1.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s2.reco.root");
  tree -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s3.reco.root");

  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;

  tree -> SetBranchAddress("STRecoTrack", &trackArray);
  tree -> SetBranchAddress("STVertex", &vertexArray);

  auto hVertexXY = new TH2D("hVertexXY",";x (mm);y (mm)",100,-100,100,100,-300,-150);
  auto hVertexZY = new TH1D("hVertexZY",";z (mm)",200,-100,1344);

  auto histPID = new TH2D("histPID",";p (MeV); dEdx (ADC/mm)",nbins,p1,p2,nbins,dedx1,dedx2);
  histPID -> SetTitleSize(0.04,"xy");
  histPID -> SetTitleOffset(1.4,"y");
  histPID -> SetTitleOffset(1.1,"x");
  histPID -> GetXaxis() -> CenterTitle();
  histPID -> GetYaxis() -> CenterTitle();

  Int_t feventID = -1;
  Int_t ftrackID = -1;
  Int_t fpid = -1;
  Double_t fp = -999;
  Double_t fdedx = -999;

  auto pidTest = new STPIDTest();

  Int_t nEvents = tree -> GetEntries();
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (iEvent % 500 == 0)
      cout << "Event " << iEvent << endl;

    tree -> GetEntry(iEvent);

    if (vertexArray -> GetEntries() != 1) continue;
    auto vertex = (STVertex *) vertexArray -> At(0);
    auto posVertex = vertex -> GetPos();

    hVertexXY -> Fill(posVertex.X(),posVertex.Y());
    hVertexZY -> Fill(posVertex.Z());

    Int_t nTracks = trackArray -> GetEntries();
    for (auto iTrack = 0; iTrack < nTracks; ++iTrack) {
      auto track = (STRecoTrack *) trackArray -> At(iTrack);
      auto pid = track -> GetPID();

      if (selectPID != STPID::kNon && pid != selectPID) continue;
      if (track -> GetVertexID() < 0) continue;

      auto dedxArray = track -> GetdEdxPointArray();
      if (dedxArray -> size() < 20) continue;

      auto p = track -> GetCharge() * track -> GetMomentum().Mag();
      auto dedx = track -> GetdEdxWithCut(0,0.7);
      histPID -> Fill(p, dedx);
    }
  }

  auto cvsPID = new TCanvas("cvsPID","",20,20,1200,700);
  histPID -> Draw("colz"); 

  for (auto ipid = 0; ipid < STPID::GetNUMSTPID(); ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    pidTest -> GetdEdxFunction(pid) -> Draw("same");
  }

  auto cvsVXY = new TCanvas("cvsVXY","",20,20,400,300);
  hVertexXY -> Draw("colz");
  auto cvsVZY = new TCanvas("cvsVZY","",420,20,400,300);
  cvsVZY -> SetLogy();
  hVertexZY -> Draw("colz");
}
