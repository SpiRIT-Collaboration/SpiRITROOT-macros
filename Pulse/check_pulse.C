/*
 * This macro draw ADC distribution and fitted pulse of pad.
 * For input, raw-data (or meta-data) is needed.
 * One can set reco-file to fChain. This will activate to draw fitted pulses.
*/



// Put the parameter file name. Path is automatically concatenated.
TString fParameterFile = "ST.parameters.Commissioning_201604.par";

// Set use the meta data files
Bool_t fUseMetadata = kTRUE;

// Set the supplement path which contains data list and meta data
// Only works when fUseMetadata is kTRUE
TString fSupplementPath = "/mnt/spirit/rawdata/misc/rawdataSupplement"; // for nscl flagtial
//TString fSupplementPath = "/data/Q16264/rawdataSupplement";

// Set use the gain calibration data file.
Bool_t fUseGainCalibration = kFALSE;

// FPN pedestal range selection threshold
Int_t fFPNThreshold = 5;

// Set the gating grid noise data. If left blank, it will use FPN pedestal.
TString fGGNoiseData = "";



STCore *fCore = nullptr;
STRawEvent *fRawEvent = nullptr;

TChain *fChain = nullptr;
TClonesArray *fHitArray = nullptr;
STPadPlaneMap *fPadPlaneMap = nullptr;
STPulse *fPulse = nullptr;
TCanvas *fCvsPad = nullptr;
TH1D *fHistPad = nullptr;



void SetCore(Int_t runNo);
void event(Int_t event); // set event
void pad(Int_t row, Int_t layer); // draw pad

void check_pulse (
  Int_t runNo = 2894, // TODO 
  Int_t eventNo = 0   // TODO, one can set event later with function event(eventNo)
)
{
  fChain = new TChain("cbmsim");

  //TODO
  fChain -> Add("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s0.reco.root");
  //for (auto iEvent = 0; iEvent < 6; ++iEvent)
    //fChain -> Add(Form("/mnt/spirit/analysis/leej/SpiRITROOT/macros/data/run2894_s%d.reco.v1.04.root",iEvent));

  if (fChain) {
    fChain -> SetBranchAddress("STHit", &fHitArray);
    fPadPlaneMap = new STPadPlaneMap();
    fPulse = new STPulse();
  }

  fCvsPad = new TCanvas("cvsPad","",800,600);
  fHistPad = new TH1D("histPad",";time bucket;ADC",512,0,512);
  fHistPad -> SetStats(0);
  fHistPad -> SetTitleOffset(1.1,"x");
  fHistPad -> SetTitleOffset(1.3,"y");

  SetCore(runNo);

  cout << "===========================================================" << endl;
  cout << " To set event: event(eventNo)" << endl;
  cout << " To draw pad ADC-dist. and pulse-fit: pad(row, layer) " << endl;
  cout << "===========================================================" << endl;

  event(eventNo);
}

void pad(Int_t row, Int_t layer)
{
  auto pad = fRawEvent -> GetPad(row, layer);
  auto adc = pad -> GetADC();

  fHistPad -> SetTitle(Form("pad (%d,%d)",row, layer));

  for (auto tb = 0; tb < 512; ++tb)
    fHistPad -> SetBinContent(tb+1, adc[tb]);
  fHistPad -> GetXaxis() -> SetRangeUser(0,280);
  fHistPad -> Draw();

  if (fChain) {
    vector<STHit *> hits;
    fPadPlaneMap -> GetPad(row, layer) -> GetHits(&hits);
    for (auto hit : hits) {
      auto pulse = fPulse -> GetPulseFunction(hit); // TF1
      pulse -> SetNpx(400);
      pulse -> Draw("lsame");

      cout << "STHit: tb(" << hit -> GetTb() << ") adc(" << hit -> GetCharge() << ")" << endl;
    }
  }
}

void event(Int_t eventNo)
{
  fRawEvent = fCore -> GetRawEvent(eventNo);
  if (fRawEvent == nullptr) {
    cout << "End of Run!" << endl;
    return;
  }

  if (fChain == nullptr)
    return;
  else if (fChain -> GetEntriesFast() <= eventNo) {
    cout << "TChain has no event-" << eventNo << endl;
    return;
  }

  fChain -> GetEntry(eventNo);

  fPadPlaneMap -> Clear();
  auto nHits = fHitArray -> GetEntriesFast();
  for (auto iHit = 0; iHit < nHits; ++iHit)
    fPadPlaneMap -> AddHit((STHit *)fHitArray->At(iHit));

  fCvsPad -> SetTitle(Form("event %d",eventNo));
  fCvsPad -> Clear();

  pad(51,23);
}

void SetCore(Int_t runNo)
{
  if (fCore != nullptr)
    return;

  TString dataFile = "";
  TString metaFile = "";
  if (fUseMetadata) {
    dataFile = Form("%s/run_%04d/dataList.txt", fSupplementPath.Data(), runNo);
    metaFile = Form("%s/run_%04d/metadataList.txt", fSupplementPath.Data(), runNo);
  } else {
    gSystem -> Exec(Form("./createList.sh %d", runNo));

    dataFile = Form("list_run%04d.txt", runNo);
  }

  TString workDir = gSystem -> Getenv("VMCWORKDIR");
  TString parameterDir = workDir + "/parameters/";

  STParReader *fPar = new STParReader(parameterDir + fParameterFile);

  Bool_t fUseSeparatedData = kFALSE;
  if (dataFile.EndsWith(".txt"))
    fUseSeparatedData = kTRUE;

  if (!fUseSeparatedData) {
    fCore = new STCore(dataFile);
  } else {
    fCore = new STCore();
    fCore -> SetUseSeparatedData(fUseSeparatedData);

    TString dataFileWithPath = dataFile;
    std::ifstream listFile(dataFileWithPath.Data());
    TString buffer;
    Int_t iCobo = -1;
    while (dataFileWithPath.ReadLine(listFile)) {
      if (dataFileWithPath.Contains("s."))
        fCore -> AddData(dataFileWithPath, iCobo);
      else {
        iCobo++;
        fCore -> AddData(dataFileWithPath, iCobo);
      }
    }
  }

  if (fUseGainCalibration) {
    fCore -> SetGainCalibrationData(fPar -> GetFilePar(fPar -> GetIntPar("GainCalibrationDataFile")));
    fCore -> SetGainReference(fPar -> GetDoublePar("GCConstant"), fPar -> GetDoublePar("GCLinear"), fPar -> GetDoublePar("GCQuadratic"));
  }

  fCore -> SetUAMap(fPar -> GetFilePar(fPar -> GetIntPar("UAMapFile")));
  fCore -> SetAGETMap(fPar -> GetFilePar(fPar -> GetIntPar("AGETMapFile")));
  fCore -> SetFPNPedestal(fFPNThreshold);
  fCore -> SetData(0);

  fCore -> SetNumTbs(fPar -> GetIntPar("NumTbs"));

  if (!fGGNoiseData.IsNull()) {
    fCore -> SetGGNoiseData(fGGNoiseData);
    fCore -> InitGGNoiseSubtractor();
  }

  if (fUseMetadata) {
    std::ifstream metalistFile(metaFile.Data());
    TString dataFileWithPath;
    for (Int_t iCobo = 0; iCobo < 12; iCobo++) {
      dataFileWithPath.ReadLine(metalistFile);
      dataFileWithPath = Form("%s/run_%04d/%s", fSupplementPath.Data(), runNo, dataFileWithPath.Data());
      fCore -> LoadMetaData(dataFileWithPath, iCobo);
    }
  }
}
