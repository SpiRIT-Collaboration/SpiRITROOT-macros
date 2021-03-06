void run_reco_experiment
(
  Int_t fRunNo = 2894,
  Int_t fSplitNo = 0,
  Int_t fNumEventsInRun = 94285,
  Int_t fNumEventsInSplit = 2,
  TString fTag = "",

  TString fParameterFile = "ST.parameters.Commissioning_201604.par",
  Double_t fPSAThreshold = 30,
  TString fGCData = "",
  TString fGGData = "",

  TString fPathToData = "/data/G17033/recoData/jungwoo/data",
  Bool_t fUseMeta = true,
  TString fSupplePath = "/data/G17033/rawdataSupplement"
)
{
  Int_t start = fSplitNo * fNumEventsInSplit;
  if (start >= fNumEventsInRun) return;
  if (start + fNumEventsInSplit > fNumEventsInRun)
    fNumEventsInSplit = fNumEventsInRun - start;

  TString sRunNo   = TString::Itoa(fRunNo, 10);
  TString sSplitNo = TString::Itoa(fSplitNo, 10);

  if (!fTag.IsNull())
    fTag = TString(".") + fTag;

  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  if (fPathToData.IsNull())
    fPathToData = spiritroot+"macros/data/";

  TString version; {
    TString name = spiritroot + "VERSION.compiled";
    std::ifstream vfile(name);
    vfile >> version;
    vfile.close();
  }
  TString par = spiritroot+"parameters/"+fParameterFile;
  TString geo = spiritroot+"geometry/geomSpiRIT.man.root";
  TString raw = TString(gSystem -> Getenv("PWD"))+"/list_run"+sRunNo+".txt";
  TString out = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco"+fTag+"."+version+".root";
  TString log = fPathToData+"run"+sRunNo+"_s"+sSplitNo+fTag+"."+version+".log";

  if (TString(gSystem -> Which(".", raw)).IsNull() && !fUseMeta)
    gSystem -> Exec("./createList.sh "+sRunNo);

  TString metaFile;
  if (fUseMeta) {
    raw = Form("%s/run_%04d/dataList.txt", fSupplePath.Data(), fRunNo);
    metaFile = Form("%s/run_%04d/metadataList.txt", fSupplePath.Data(), fRunNo);
  }

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(out);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  STDecoderTask *decoder = new STDecoderTask();
  decoder -> SetUseSeparatedData(true);
  decoder -> SetPersistence(false);
  if (fGCData.IsNull())
    decoder -> SetUseGainCalibration(false);
  else {
    decoder -> SetUseGainCalibration(true);
    decoder -> SetGainCalibrationData(fGCData);
  }
  decoder -> SetGGNoiseData(fGGData);
  decoder -> SetDataList(raw);
  decoder -> SetEventID(start);

  if (fUseMeta) {
    std::ifstream metalistFile(metaFile.Data());
    TString dataFileWithPath;
    for (Int_t iCobo = 0; iCobo < 12; iCobo++) {
      dataFileWithPath.ReadLine(metalistFile);
      dataFileWithPath = Form("%s/run_%04d/%s", fSupplePath.Data(), fRunNo, dataFileWithPath.Data());
      decoder -> SetMetaData(dataFileWithPath, iCobo);
    }
  }

  auto preview = new STEventPreviewTask();
  preview -> SetPersistence(true);

  auto psa = new STPSAETask();
  psa -> SetPersistence(false);
  psa -> SetThreshold(fPSAThreshold);
  psa -> SetLayerCut(-1, 90);
  psa -> SetPulserData("pulser_117ns.dat");

  auto helix = new STHelixTrackingTask();
  helix -> SetPersistence(false);
  helix -> SetClusterPersistence(false);
  helix -> SetClusteringOption(2);

  auto st_genfit = new STGenfitPIDTask();
  st_genfit -> SetPersistence(true);
  st_genfit -> SetConstantField();
  
  auto pidmatching = new STPIDCorrelatorTask(false);

  auto st_genfit2 = new STGenfitSinglePIDTask();
  st_genfit2 -> SetPersistence(true);

  run -> AddTask(decoder);
  run -> AddTask(preview);
  run -> AddTask(psa);
  run -> AddTask(helix);
  run -> AddTask(st_genfit);

  auto outFile = FairRootManager::Instance() -> GetOutFile();
  auto recoHeader = new STRecoHeader("RecoHeader","");
  recoHeader -> SetPar("version", version);
  recoHeader -> SetPar("eventStart", start);
  recoHeader -> SetPar("numEvents", fNumEventsInSplit);
  recoHeader -> SetPar("parameter", fParameterFile);
  recoHeader -> SetPar("GCData", fGCData);
  recoHeader -> SetPar("GGData", fGGData);
  recoHeader -> Write("RecoHeader");

  run -> Init();

  run -> Run(0, fNumEventsInSplit);

  cout << "Log    : " << log << endl;
  cout << "Input  : " << raw << endl;
  cout << "Output : " << out << endl;

  gApplication -> Terminate();
}
