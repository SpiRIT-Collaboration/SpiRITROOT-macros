/*
 * Macro to check pid probability density function in p vs dEdx plot.
 * Data are from simulation result.
 * List of supported partices are pion, proton, deuteron, triton, 3He, 4He.
 *
 * Jung Woo
 */

void pid_pdf
(
  Int_t nbins = 200,
  Double_t x1 = 0,
  Double_t x2 = 2500,
  Double_t y1 = 0,
  Double_t y2 = 800,
  Double_t probabilityAcceptance = 0.3,
  bool saveFigures = false,
  TString figurePath = "figures" // the path should exist
)
{
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.10);
  gStyle -> SetPadLeftMargin(0.11);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetTitleFontSize(0.06);

  Double_t dx = (x2-x1)/nbins;
  Double_t dy = (y2-y1)/nbins;
  figurePath = figurePath+"/";

  auto numPID = STPID::GetNUMSTPID();
  auto pidTest = new STPIDTest();

  TString pdfMapName = Form("PID_map_acceptance(%.1f)",probabilityAcceptance);
  auto pidMap = new TH2D("pid_map",pdfMapName+";p (MeV/c); dEdx (ADC/mm)",nbins,x1,x2,nbins,y1,y2);

  auto FillHist = [pidTest, nbins, dx, dy](TH2D *hist, STPID::PID stpid) {
    for (auto imom = 0; imom < nbins; ++imom) {
      auto mom = imom * dx;
      for (auto idedx = 0; idedx < nbins; ++idedx) {
        auto dedx = idedx * dy;
        hist -> Fill(mom,dedx,pidTest->GetProbability(stpid,mom,dedx));
      }
    }
  };

  TH2D *histPDF[numPID];
  for (auto ipid = 0; ipid < numPID; ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    TString pdfName = TString("PDF_")+STPID::GetName(pid);

    histPDF[ipid] = new TH2D(Form("histPDF%d",ipid),pdfName+";p (MeV/c); dEdx (ADC/mm)",nbins,x1,x2,nbins,y1,y2);

    auto cvsPID = new TCanvas(Form("cvs%d",ipid),pdfName,(ipid+1)*10,(ipid+1)*10,800,600);
    FillHist(histPDF[ipid], pid);
    histPDF[ipid] -> Draw("colz");

    if (saveFigures)
      cvsPID -> SaveAs(figurePath+pdfName+".png");
  }

  for (auto x = x1; x < x2; x+=dx) {
    for (auto y = y1; y < y2; y+=dy) {
      auto bin = pidMap -> FindBin(x, y);

      Int_t index = -1;
      Double_t maxProbability = 0.; 

      for (auto ipid = 0; ipid < numPID; ++ipid) {
        auto probability = histPDF[ipid] -> GetBinContent(bin);
        if (maxProbability < probability) {
          maxProbability = probability;
          index = ipid;
        }
      }
      if (index == -1)
        continue;

      if (maxProbability > probabilityAcceptance)
        pidMap -> SetBinContent(bin, index+1);
    }
  }

  auto cvs = new TCanvas("cvs",pdfMapName,800,600);
  pidMap -> Draw("col");

  for (auto ipid = 0; ipid < numPID; ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    pidTest -> GetdEdxFunction(pid) -> Draw("same");
  }

  if (saveFigures)
    cvs -> SaveAs(figurePath+pdfMapName+".png");
}
