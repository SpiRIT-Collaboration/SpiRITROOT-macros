auto mpi = 139.57018;
auto mp = 938.2720813;
auto mn = 939.565346;
auto md = 1875.612762;
auto mt = 2808.921112;
auto mhe3 = 2808.39132;
auto mal = 3727.379378;

void dedxSn132() {
  Int_t runs[] = {2843, 2844, 2845, 2846, 2848, 2849, 2850, 2851, 2852, // 9
                  2855, 2856, 2857, 2858, 2859, 2860, 2861,             // 7
                  2877, 2878, 2879, 2880, 2881, 2882, 2883, 2884,       // 8
                  2887, 2888, 2889, 2890, 2891, 2892, 2893, 2894};      // 8
  Int_t numFiles[] = {42, 10, 38, 38, 5, 37, 36, 38, 42,
                      38, 6, 43, 45, 38, 43, 43,
                      46, 44, 37, 41, 38, 44, 46, 31,
                      46, 46, 45, 44, 38, 7, 37, 48}; 
  Int_t numData = sizeof(runs)/sizeof(Int_t);
  TChain **chain = new TChain*[numData];
  TChain **beam = new TChain*[numData];
  for (Int_t i = 0; i < numData; i++) {
    chain[i] = new TChain("cbmsim");
    beam[i] = new TChain("beam");
  }

  auto pCutFile = new TFile("pdFuncCut.root");
  auto pFuncCut = (TCutG *) pCutFile -> FindObjectAny("pFuncCut");

  auto piCutFile = new TFile("pionCut.root");
  auto pimCut = (TCutG *) piCutFile -> FindObjectAny("pim");
  auto pipCut = (TCutG *) piCutFile -> FindObjectAny("pip");

  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;
  STEventHeader *header = nullptr;
  Bool_t sigma10, sigma15, sigma20, proton, pip, pim, pipbg;
  Double_t vx, vy, vz, projmevu, projpx, projpy, projpz, proja, projb, projx, projy, projz;

  for (Int_t iRun = 0; iRun < numData; iRun++) {
    for (Int_t iFile = 0; iFile < numFiles[iRun]; iFile++) {
      chain[iRun] -> AddFile(Form("Sn132-LayerCut90/run%d_s%d.reco.v1.04.root", runs[iRun], iFile));
    }
    beam[iRun] -> AddFile(Form("/mnt/spirit/analysis/changj/BeamAnalysis/macros/refined/Sn132/run%d.refined.root", runs[iRun]));
  }

  for (Int_t i = 0; i < numData; i++) {
    chain[i] -> SetBranchAddress("STEventHeader", &header);
    chain[i] -> SetBranchAddress("STRecoTrack", &trackArray);
    chain[i] -> SetBranchAddress("STVertex", &vertexArray);
    beam[i] -> SetBranchAddress("sigma10", &sigma10);
    beam[i] -> SetBranchAddress("sigma15", &sigma15);
    beam[i] -> SetBranchAddress("sigma20", &sigma20);
    beam[i] -> SetBranchAddress("projmevu", &projmevu);
    beam[i] -> SetBranchAddress("projpx", &projpx);
    beam[i] -> SetBranchAddress("projpy", &projpy);
    beam[i] -> SetBranchAddress("projpz", &projpz);
    beam[i] -> SetBranchAddress("proja", &proja);
    beam[i] -> SetBranchAddress("projb", &projb);
    beam[i] -> SetBranchAddress("projx", &projx);
    beam[i] -> SetBranchAddress("projy", &projy);
    beam[i] -> SetBranchAddress("projz", &projz);
  }

  auto file = new TFile("dedxSn132-LayerCut90.root", "recreate");
  Double_t dedx, mom, ndf, dx, dy, dz, dist, pt, pzCM, rapL, rapCM, KECM, phiL, thetaL, phiCM, thetaCM;
  Int_t charge, pid, run, eventid, vid, parentvid, mult;
  auto tree = new TTree("dedx", "");
  tree -> Branch("run", &run);
  tree -> Branch("eventid", &eventid);
  tree -> Branch("dedx", &dedx);
  tree -> Branch("mom", &mom);
  tree -> Branch("charge", &charge);
  tree -> Branch("pid", &pid);
  tree -> Branch("ndf", &ndf);
  tree -> Branch("dx", &dx); // mom direction
  tree -> Branch("dy", &dy);
  tree -> Branch("dz", &dz);
  tree -> Branch("vid", &vid); // vertex id 
  tree -> Branch("parentvid", &parentvid); // track's parent vertex id 
  tree -> Branch("vx", &vx); // vertex position
  tree -> Branch("vy", &vy);
  tree -> Branch("vz", &vz);
  tree -> Branch("dist", &dist); // poca and vertex distance
  tree -> Branch("sigma10", &sigma10);
  tree -> Branch("sigma15", &sigma15);
  tree -> Branch("sigma20", &sigma20);
  tree -> Branch("mult", &mult);
  tree -> Branch("pt", &pt);
  tree -> Branch("pzCM", &pzCM);
  tree -> Branch("KECM", &KECM);
  tree -> Branch("rapL", &rapL);
  tree -> Branch("rapCM", &rapCM);
  tree -> Branch("projx", &projx);
  tree -> Branch("projy", &projy);
  tree -> Branch("projz", &projz);
  tree -> Branch("phiL", &phiL);
  tree -> Branch("thetaL", &thetaL);
  tree -> Branch("phiCM", &phiCM);
  tree -> Branch("thetaCM", &thetaCM);
  tree -> Branch("proton", &proton);
  tree -> Branch("pip", &pip);
  tree -> Branch("pim", &pim);

  Double_t amu2mev = 931.494028;
  Double_t beamZ = 50;
  Double_t beamN = 82;
  Double_t beamA = beamZ + beamN;
  Double_t targetZ = 50;
  Double_t targetN = 74;
  Double_t targetA = targetZ + targetN;
  Double_t aMass = 0.5*(mp + mn) - 8; // MeV

  for (auto iRun = 0; iRun < numData; iRun++) {
    run = runs[iRun];
    auto numEntries = chain[iRun] -> GetEntries();
    cout << "Run: " << runs[iRun] << " Entries: " << numEntries << endl;
    for (auto iEntry = 0; iEntry < numEntries; iEntry++) {
      proton = 0;
      pim = 0;
      pip = 0;
      pipbg = 0;

      beam[iRun] -> GetEntry(iEntry);
      chain[iRun] -> GetEntry(iEntry);
      eventid = header -> GetEventID();

      Double_t beamKE = projmevu*0.99917; // MeV
      Double_t beamM = beamA*aMass;
      Double_t beamE = beamKE*beamA + beamM;
      Double_t beamP = TMath::Sqrt(beamE*beamE - beamM*beamM); // MeV/c
      Double_t beamRapidity = 0.5*log((beamE + beamP)/(beamE - beamP));
      Double_t E_CM = beamE + targetA*aMass; 
      Double_t P_CM = beamP;
      Double_t rapidity_CM = 0.5*log((E_CM + P_CM)/(E_CM - P_CM));
      Double_t beta = -(P_CM/E_CM);
      TVector3 betaVec(0., 0., beta);

      if (iEntry%10000 == 0)
        cout << "Current: " << iEntry << "/" << numEntries << endl;

      auto numTracks = trackArray -> GetEntries();
      mult = numTracks;
      for (auto iTrack = 0; iTrack < numTracks; iTrack++) {
        auto track = (STRecoTrack *) trackArray -> At(iTrack);

        parentvid = track -> GetVertexID();
        if (parentvid < 0) {
          vid = -99998;
          vx = -99998;
          vy = -99998;
          vz = -99998;
          dist = -99998;
        } else {
          vid = parentvid;
          auto vertex = (STVertex *) vertexArray -> At(parentvid);
          auto posVertex = vertex -> GetPos();

          vx = posVertex.X();
          vy = posVertex.Y();
          vz = posVertex.Z();

          auto pocaVertex = track -> GetPOCAVertex();
          dist = (pocaVertex - posVertex).Mag();
        }

        dedx = track -> GetdEdxWithCut(0, 0.7);
        charge = track -> GetCharge();
        mom = track -> GetMomentum().Mag();
        pid = STPID::GetPDG(track -> GetPID());

        Double_t mass = 0;
        switch (pid) {
          case 2212:
            mass = mp;
            break;

          case 211:
          case -211:
            mass = mpi;
            break;

          case 1000010020:
            mass = md;
            break;

          case 1000010030:
            mass = mt;
            break;

          case 1000020030:
            charge = 2;
            mass = mhe3;
            break;

          case 1000020040:
            charge = 2;
            mass = mal;
            break;

          default:
            mass = 0;
            break;
        }

        auto momVec = track -> GetMomentum();

        if (momVec.Z() < 0)
          momVec = -momVec;

        momVec.RotateY(-proja/1000.);
        momVec.RotateX(projb/1000.);
        pt = momVec.Pt();
        Double_t e = TMath::Sqrt(mass*mass + mom*mom);
        TLorentzVector vec(momVec, e);
        rapL = vec.Rapidity();
        vec.Boost(betaVec);
        pzCM = vec.Pz();
        KECM = vec.E() - vec.M();
        rapCM = vec.Rapidity()/(beamRapidity - rapidity_CM);

        momVec = momVec.Unit();
        ndf = track -> GetClusterIDArray() -> size();
        dx = momVec.X();
        dy = momVec.Y();
        dz = momVec.Z();
        phiL = momVec.Phi()*180./TMath::Pi();
        phiL = (phiL < 0 ? phiL + 360 : phiL);
        thetaL = momVec.Theta()*180./TMath::Pi();
        phiCM = vec.Phi()*180./TMath::Pi();
        phiCM = (phiCM < 0 ? phiCM + 360 : phiCM);
        thetaCM = momVec.Theta()*180./TMath::Pi();

        proton = pFuncCut -> IsInside(mom/charge, dedx);
        pim = pimCut -> IsInside(mom/charge, dedx);
        pip = pipCut -> IsInside(mom/charge, dedx);

        tree -> Fill();
      }
    }
  }

  file -> cd();
  tree -> Write();
}
