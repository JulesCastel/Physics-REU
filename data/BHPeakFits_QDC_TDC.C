// BHPeakFits_QDC_TDC.C
//
// ROOT macro for fitting Gaussian peaks in BH QDC/TDC histograms.
//
// It iterates over run35566.root ... run35572.root by default.
// For every run, it loops over:
//   BHC paddles 00--15
//   BHD paddles 00--12
// and fits:
//   QDC/qdc_hit_[BHC/BHD][paddle#]left
//   TDC/tdc_coin_left_[BHC/BHD][paddle#]-RF
//
// Output:
//   1. CSV table with peak positions and uncertainties.
//   2. Individual PNG plots for every successful/attempted histogram.
//   3. Multipage PDFs grouped by run, plane, and QDC/TDC type.
//
// Example usage from a terminal:
//   root -l -b -q 'BHPeakFits_QDC_TDC.C(".", "BH_peak_output")'
//
// Example usage inside ROOT:
//   .L BHPeakFits_QDC_TDC.C
//   BHPeakFits_QDC_TDC(".", "BH_peak_output");
//
// If your files live somewhere else, replace "." with that directory,
// for example:
//   BHPeakFits_QDC_TDC("/data2/processed_fast_reduced/BH_detail", "BH_peak_output");

#include "TCanvas.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

struct PeakFitResult {
    int run = -1;
    TString plane = "";
    int paddle = -1;
    TString type = "";       // "QDC" or "TDC"
    TString histName = "";
    TString histPath = "";
    bool found = false;
    bool ok = false;
    int fitStatus = -999;

    double entries = 0.0;
    double amplitude = std::numeric_limits<double>::quiet_NaN();
    double amplitudeErr = std::numeric_limits<double>::quiet_NaN();
    double peak = std::numeric_limits<double>::quiet_NaN();       // Gaussian mean
    double peakErr = std::numeric_limits<double>::quiet_NaN();    // uncertainty on Gaussian mean
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigmaErr = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double ndf = std::numeric_limits<double>::quiet_NaN();
    double fitLow = std::numeric_limits<double>::quiet_NaN();
    double fitHigh = std::numeric_limits<double>::quiet_NaN();

    TString message = "";
};

void MakeDirectoryIfNeeded(const TString& dir)
{
    if (dir.Length() > 0) {
        gSystem->mkdir(dir.Data(), kTRUE);
    }
}

TString BuildHistName(const TString& plane, int paddle, const TString& type)
{
    TString paddleTag = Form("%s%02d", plane.Data(), paddle);

    if (type == "QDC") {
        return Form("qdc_hit_%sleft", paddleTag.Data());
    }

    if (type == "TDC") {
        return Form("tdc_coin_left_%s-RF", paddleTag.Data());
    }

    return "";
}

TH1* AsHistogram(TObject* obj)
{
    if (obj == nullptr) {
        return nullptr;
    }

    if (obj->InheritsFrom(TH1::Class())) {
        return dynamic_cast<TH1*>(obj);
    }

    return nullptr;
}

// Recursive fallback search. This is useful if the ROOT file directory layout is
// slightly different from the assumed path, but the histogram name is correct.
TH1* FindHistRecursive(TDirectory* dir, const TString& wantedName, TString& foundPath)
{
    if (dir == nullptr) {
        return nullptr;
    }

    TIter nextKey(dir->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = dynamic_cast<TKey*>(nextKey()))) {
        TString keyName = key->GetName();

        // Read the object represented by this key.
        TObject* obj = key->ReadObj();
        if (obj == nullptr) {
            continue;
        }

        TString currentPath = Form("%s/%s", dir->GetPath(), keyName.Data());

        if (keyName == wantedName && obj->InheritsFrom(TH1::Class())) {
            foundPath = currentPath;
            return dynamic_cast<TH1*>(obj);
        }

        if (obj->InheritsFrom(TDirectory::Class())) {
            TH1* h = FindHistRecursive(dynamic_cast<TDirectory*>(obj), wantedName, foundPath);
            if (h != nullptr) {
                return h;
            }
        }
    }

    return nullptr;
}

TH1* GetPaddleHistogram(TFile* file,
                        const TString& plane,
                        int paddle,
                        const TString& type,
                        TString& usedPath)
{
    if (file == nullptr || file->IsZombie()) {
        return nullptr;
    }

    TString histName = BuildHistName(plane, paddle, type);
    TString paddleTag = Form("%s%02d", plane.Data(), paddle);

    std::vector<TString> candidatePaths;
    candidatePaths.push_back(Form("%s/%s/%s/%s", plane.Data(), paddleTag.Data(), type.Data(), histName.Data()));
    candidatePaths.push_back(Form("%s/%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));
    candidatePaths.push_back(Form("%s/Paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));
    candidatePaths.push_back(Form("%s/paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));
    candidatePaths.push_back(Form("%s/%s", plane.Data(), histName.Data()));

    // First try the most likely direct paths.
    for (const auto& path : candidatePaths) {
        TObject* obj = file->Get(path.Data());
        TH1* h = AsHistogram(obj);
        if (h != nullptr) {
            usedPath = path;
            return h;
        }
    }

    // Then try a recursive search under the plane directory.
    TDirectory* planeDir = dynamic_cast<TDirectory*>(file->Get(plane.Data()));
    if (planeDir != nullptr) {
        TH1* h = FindHistRecursive(planeDir, histName, usedPath);
        if (h != nullptr) {
            return h;
        }
    }

    // Last resort: recursive search in the whole file.
    TH1* h = FindHistRecursive(file, histName, usedPath);
    if (h != nullptr) {
        return h;
    }

    usedPath = "";
    return nullptr;
}

bool IsFinite(double x)
{
    return std::isfinite(x);
}

// Estimate a fit window centered on the strongest peak.
// This avoids hard-coding QDC/TDC ranges for every run and paddle.
void EstimateGaussianWindow(TH1* h,
                            double& peakX,
                            double& sigmaGuess,
                            double& fitLow,
                            double& fitHigh)
{
    const int nbins = h->GetNbinsX();
    const int peakBin = h->GetMaximumBin();

    const double xmin = h->GetXaxis()->GetXmin();
    const double xmax = h->GetXaxis()->GetXmax();
    const double fullRange = xmax - xmin;
    const double binWidth = h->GetXaxis()->GetBinWidth(peakBin);

    peakX = h->GetBinCenter(peakBin);
    const double peakY = h->GetBinContent(peakBin);
    const double halfMax = 0.5 * peakY;

    int leftBin = peakBin;
    int rightBin = peakBin;

    while (leftBin > 1 && h->GetBinContent(leftBin) > halfMax) {
        leftBin--;
    }
    while (rightBin < nbins && h->GetBinContent(rightBin) > halfMax) {
        rightBin++;
    }

    const double fwhmLow = h->GetBinLowEdge(leftBin);
    const double fwhmHigh = h->GetBinLowEdge(rightBin + 1);
    double fwhm = fwhmHigh - fwhmLow;

    // If the FWHM estimate is bad, fall back to a conservative fraction of the axis.
    if (!IsFinite(fwhm) || fwhm <= 2.0 * binWidth || fwhm > 0.8 * fullRange) {
        fwhm = 0.10 * fullRange;
    }

    sigmaGuess = fwhm / 2.355;
    if (!IsFinite(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = std::max(h->GetRMS(), 5.0 * binWidth);
    }

    if (!IsFinite(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = 0.05 * fullRange;
    }

    double window = std::max(3.0 * sigmaGuess, 6.0 * binWidth);

    // Keep the automatic fit from accidentally covering almost the whole histogram.
    // If the peak is broad, this still leaves a broad enough window for the fit.
    window = std::min(window, 0.30 * fullRange);

    fitLow = std::max(xmin, peakX - window);
    fitHigh = std::min(xmax, peakX + window);

    if (fitHigh <= fitLow) {
        fitLow = xmin;
        fitHigh = xmax;
    }
}

void DrawFitText(PeakFitResult& r)
{
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextSize(0.032);

    double y = 0.86;
    latex.DrawLatex(0.58, y, Form("Peak = %.6g #pm %.3g", r.peak, r.peakErr));
    y -= 0.045;
    latex.DrawLatex(0.58, y, Form("#sigma = %.6g #pm %.3g", r.sigma, r.sigmaErr));
    y -= 0.045;

    if (r.ndf > 0.0) {
        latex.DrawLatex(0.58, y, Form("#chi^{2}/ndf = %.3g / %.0f", r.chi2, r.ndf));
    } else {
        latex.DrawLatex(0.58, y, "#chi^{2}/ndf unavailable");
    }
}

PeakFitResult FitAndPlotHistogram(TH1* h,
                                  int run,
                                  const TString& plane,
                                  int paddle,
                                  const TString& type,
                                  const TString& histPath,
                                  const TString& outputDir,
                                  const TString& pdfFile)
{
    PeakFitResult r;
    r.run = run;
    r.plane = plane;
    r.paddle = paddle;
    r.type = type;
    r.histName = BuildHistName(plane, paddle, type);
    r.histPath = histPath;
    r.found = (h != nullptr);

    TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    TString plotName = Form("run%d_%s_%s", run, paddleTag.Data(), type.Data());

    TCanvas* c = new TCanvas(Form("c_%s", plotName.Data()),
                             Form("%s Gaussian fit", plotName.Data()),
                             950, 700);
    c->SetGrid();

    if (h == nullptr) {
        r.message = "histogram not found";
        c->cd();
        TLatex latex;
        latex.SetNDC(kTRUE);
        latex.SetTextSize(0.040);
        latex.DrawLatex(0.12, 0.75, Form("Missing histogram: %s", r.histName.Data()));
        latex.DrawLatex(0.12, 0.68, Form("run %d, %s paddle %02d, %s", run, plane.Data(), paddle, type.Data()));

        TString png = Form("%s/%s_MISSING.png", outputDir.Data(), plotName.Data());
        c->SaveAs(png.Data());
        if (pdfFile.Length() > 0) {
            c->Print(pdfFile.Data());
        }

        delete c;
        return r;
    }

    r.entries = h->GetEntries();

    if (r.entries < 10) {
        r.message = "too few entries";
        h->SetTitle(Form("run %d %s%02d %s: %s;ADC/TDC channel;Counts",
                         run, plane.Data(), paddle, type.Data(), r.histName.Data()));
        h->Draw("hist");

        TLatex latex;
        latex.SetNDC(kTRUE);
        latex.SetTextSize(0.040);
        latex.DrawLatex(0.15, 0.82, "Too few entries for a reliable Gaussian fit.");

        TString png = Form("%s/%s_TOO_FEW_ENTRIES.png", outputDir.Data(), plotName.Data());
        c->SaveAs(png.Data());
        if (pdfFile.Length() > 0) {
            c->Print(pdfFile.Data());
        }

        delete c;
        return r;
    }

    const int peakBin = h->GetMaximumBin();
    const double peakY = h->GetBinContent(peakBin);

    if (peakY <= 0.0) {
        r.message = "histogram maximum is non-positive";
        h->Draw("hist");

        TString png = Form("%s/%s_NO_POSITIVE_PEAK.png", outputDir.Data(), plotName.Data());
        c->SaveAs(png.Data());
        if (pdfFile.Length() > 0) {
            c->Print(pdfFile.Data());
        }

        delete c;
        return r;
    }

    double peakGuess = 0.0;
    double sigmaGuess = 0.0;
    EstimateGaussianWindow(h, peakGuess, sigmaGuess, r.fitLow, r.fitHigh);

    TString fitName = Form("gaus_%d_%s_%02d_%s", run, plane.Data(), paddle, type.Data());
    TF1* fGaus = new TF1(fitName.Data(), "gaus", r.fitLow, r.fitHigh);
    fGaus->SetParameters(peakY, peakGuess, sigmaGuess);
    fGaus->SetParNames("Amplitude", "Mean", "Sigma");

    const double binWidth = h->GetXaxis()->GetBinWidth(peakBin);
    const double fullRange = h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin();
    fGaus->SetParLimits(2, std::max(0.05 * binWidth, 1.0e-9), fullRange);

    // Q = quiet, R = use fit range, 0 = do not draw during fit, S = return fit result.
    TFitResultPtr fitPtr = h->Fit(fGaus, "QR0S");
    r.fitStatus = static_cast<int>(fitPtr);

    r.amplitude = fGaus->GetParameter(0);
    r.amplitudeErr = fGaus->GetParError(0);
    r.peak = fGaus->GetParameter(1);
    r.peakErr = fGaus->GetParError(1);
    r.sigma = std::fabs(fGaus->GetParameter(2));
    r.sigmaErr = fGaus->GetParError(2);
    r.chi2 = fGaus->GetChisquare();
    r.ndf = fGaus->GetNDF();

    r.ok = (r.fitStatus == 0 &&
            IsFinite(r.peak) &&
            IsFinite(r.peakErr) &&
            r.peakErr > 0.0 &&
            IsFinite(r.sigma) &&
            r.sigma > 0.0);

    if (!r.ok) {
        r.message = Form("fit completed with status %d; inspect plot", r.fitStatus);
    } else {
        r.message = "ok";
    }

    c->cd();
    h->SetLineWidth(2);
    h->SetTitle(Form("run %d %s%02d %s: %s;%s channel;Counts",
                     run,
                     plane.Data(),
                     paddle,
                     type.Data(),
                     r.histName.Data(),
                     type.Data()));
    h->Draw("hist");

    fGaus->SetLineColor(kRed + 1);
    fGaus->SetLineWidth(3);
    fGaus->Draw("same");

    TLine* lowLine = new TLine(r.fitLow, 0.0, r.fitLow, peakY);
    TLine* highLine = new TLine(r.fitHigh, 0.0, r.fitHigh, peakY);
    lowLine->SetLineStyle(2);
    highLine->SetLineStyle(2);
    lowLine->SetLineColor(kGray + 2);
    highLine->SetLineColor(kGray + 2);
    lowLine->Draw("same");
    highLine->Draw("same");

    DrawFitText(r);

    TString statusLabel = r.ok ? "OK" : "CHECK";
    TString png = Form("%s/%s_%s.png", outputDir.Data(), plotName.Data(), statusLabel.Data());
    c->SaveAs(png.Data());

    if (pdfFile.Length() > 0) {
        c->Print(pdfFile.Data());
    }

    delete c;

    return r;
}

void WriteCsvRow(std::ofstream& csv, const PeakFitResult& r)
{
    csv << r.run << ","
        << r.plane << ","
        << std::setw(2) << std::setfill('0') << r.paddle << std::setfill(' ') << ","
        << r.type << ","
        << r.histName << ","
        << r.histPath << ","
        << (r.found ? 1 : 0) << ","
        << (r.ok ? 1 : 0) << ","
        << r.fitStatus << ","
        << std::setprecision(12)
        << r.entries << ","
        << r.peak << ","
        << r.peakErr << ","
        << r.sigma << ","
        << r.sigmaErr << ","
        << r.amplitude << ","
        << r.amplitudeErr << ","
        << r.chi2 << ","
        << r.ndf << ","
        << r.fitLow << ","
        << r.fitHigh << ","
        << r.message
        << "\n";
}

void BHPeakFits_QDC_TDC(const char* inputDir = ".",
                        const char* outputDir = "BH_peak_output",
                        int firstRun = 35566,
                        int lastRun = 35572)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);

    TString outDir(outputDir);
    MakeDirectoryIfNeeded(outDir);

    TString csvName = Form("%s/BH_QDC_TDC_peak_summary_runs%d_%d.csv",
                           outDir.Data(), firstRun, lastRun);

    std::ofstream csv(csvName.Data());
    if (!csv.is_open()) {
        std::cerr << "ERROR: could not open CSV output file: " << csvName << std::endl;
        return;
    }

    csv << "run,plane,paddle,type,histogram,hist_path,found,fit_ok,fit_status,"
        << "entries,peak,peak_uncertainty,sigma,sigma_uncertainty,"
        << "amplitude,amplitude_uncertainty,chi2,ndf,fit_low,fit_high,message\n";

    const std::vector<TString> planes = {"BHC", "BHD"};
    const std::vector<TString> types = {"QDC", "TDC"};

    int totalRequested = 0;
    int totalFound = 0;
    int totalGoodFits = 0;

    for (int run = firstRun; run <= lastRun; run++) {
        TString fileName = Form("%s/run%d.root", inputDir, run);
        TFile* file = TFile::Open(fileName.Data(), "READ");

        if (file == nullptr || file->IsZombie()) {
            std::cerr << "WARNING: could not open " << fileName << std::endl;
            if (file != nullptr) {
                file->Close();
                delete file;
            }
            continue;
        }

        std::cout << "\nProcessing " << fileName << std::endl;

        TString runDir = Form("%s/run%d", outDir.Data(), run);
        MakeDirectoryIfNeeded(runDir);

        for (const auto& plane : planes) {
            int nPaddles = (plane == "BHC") ? 16 : 13;

            TString planeDir = Form("%s/%s", runDir.Data(), plane.Data());
            MakeDirectoryIfNeeded(planeDir);

            for (const auto& type : types) {
                TString typeDir = Form("%s/%s", planeDir.Data(), type.Data());
                MakeDirectoryIfNeeded(typeDir);

                TString pdfFile = Form("%s/run%d_%s_%s_gaussian_fits.pdf",
                                       typeDir.Data(),
                                       run,
                                       plane.Data(),
                                       type.Data());

                TCanvas pdfBookStarter("pdfBookStarter", "pdfBookStarter", 10, 10);
                pdfBookStarter.Print(Form("%s[", pdfFile.Data()));

                for (int paddle = 0; paddle < nPaddles; paddle++) {
                    totalRequested++;

                    TString usedPath = "";
                    TH1* h = GetPaddleHistogram(file, plane, paddle, type, usedPath);

                    PeakFitResult r = FitAndPlotHistogram(h,
                                                          run,
                                                          plane,
                                                          paddle,
                                                          type,
                                                          usedPath,
                                                          typeDir,
                                                          pdfFile);

                    if (r.found) {
                        totalFound++;
                    }
                    if (r.ok) {
                        totalGoodFits++;
                    }

                    WriteCsvRow(csv, r);

                    std::cout << "  run " << run
                              << " " << plane << Form("%02d", paddle)
                              << " " << type
                              << " : peak = " << r.peak
                              << " +/- " << r.peakErr
                              << "  [" << r.message << "]"
                              << std::endl;
                }

                pdfBookStarter.Print(Form("%s]", pdfFile.Data()));
            }
        }

        file->Close();
        delete file;
    }

    csv.close();

    std::cout << "\nDone.\n"
              << "Requested histograms: " << totalRequested << "\n"
              << "Found histograms:     " << totalFound << "\n"
              << "Good fits:            " << totalGoodFits << "\n"
              << "CSV summary:          " << csvName << "\n"
              << "Plots directory:      " << outDir << "\n"
              << std::endl;
}
