// BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded.C
//
// Memory-safe ROOT macro for Gaussian peak fits in BH QDC/TDC histograms, adjusted for the run35567.root directory layout.
// This version formats all finite non-integer numeric output to the hundredths place.
//
// Main memory fixes compared with the first version:
//   1. Does NOT recursively read every object in the ROOT file by default.
//   2. If recursive search is enabled, it searches TKeys first and only reads
//      the matching histogram instead of loading every histogram into memory.
//   3. Reuses one TCanvas per plot book instead of creating hundreds of canvases.
//   4. Uses fit option "N" so ROOT does not attach every fit function to the histogram.
//   5. Deletes temporary TF1 objects and clears the canvas after every plot.
//
// Usage, one run with plots:
//   root -l -b -q 'BHPeakFits_QDC_TDC_memsafe_fixedpaths.C(".", "BH_peak_output_35566", 35566, 35566)'
//
// For run35567.root, the relevant structure is:
//   BHC/Paddle00/QDC/qdc_hit_BHC00left
//   BHC/Paddle00/TDC/tdc_coinc_left_BHC00-RF
//   BHD/Paddle00/QDC/qdc_hit_BHD00down
//   BHD/Paddle00/TDC/tdc_coinc_down_BHD00-RF
// Note: the TDC names use "tdc_coinc", not "tdc_coin".
//       BHD uses down/up readout names, not left/right.
//
// Usage, one run CSV only, no plots:
//   root -l -b -q 'BHPeakFits_QDC_TDC_memsafe_fixedpaths.C(".", "BH_peak_csv_35566", 35566, 35566, false)'
//
// If histograms are not found because the directory structure is different,
// enable the memory-safe recursive search as the 6th argument:
//   root -l -b -q 'BHPeakFits_QDC_TDC_memsafe_fixedpaths.C(".", "BH_peak_recursive_35566", 35566, 35566, false, true)'

#include "TCanvas.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResultPtr.h"
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

struct PeakFitResultMemSafe {
    int run = -1;
    TString plane = "";
    int paddle = -1;
    TString type = "";
    TString histName = "";
    TString histPath = "";
    bool found = false;
    bool ok = false;
    int fitStatus = -999;

    double entries = 0.0;
    double amplitude = std::numeric_limits<double>::quiet_NaN();
    double amplitudeErr = std::numeric_limits<double>::quiet_NaN();
    double peak = std::numeric_limits<double>::quiet_NaN();
    double peakErr = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigmaErr = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double ndf = std::numeric_limits<double>::quiet_NaN();
    double fitLow = std::numeric_limits<double>::quiet_NaN();
    double fitHigh = std::numeric_limits<double>::quiet_NaN();

    TString message = "";
};

void MakeDirectoryIfNeededMemSafe(const TString& dir)
{
    if (dir.Length() > 0) {
        gSystem->mkdir(dir.Data(), kTRUE);
    }
}

bool IsFiniteMemSafe(double x)
{
    return std::isfinite(x);
}


TString FormatNumberHundredthMemSafe(double x)
{
    // Keep missing/failed fits obvious in the CSV and terminal output.
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return (x > 0.0) ? "inf" : "-inf";

    // Round finite values to the hundredths place for display/output.
    double rounded = std::round(x * 100.0) / 100.0;

    // Avoid printing tiny negative zeros such as -0.00.
    if (std::fabs(rounded) < 0.005) rounded = 0.0;

    // Integer-valued quantities such as entries and ndf should remain integers.
    double nearestInteger = std::round(rounded);
    if (std::fabs(rounded - nearestInteger) < 1.0e-9) {
        return Form("%.0f", rounded);
    }

    return Form("%.2f", rounded);
}

TString DefaultSideForPlaneMemSafe(const TString& plane)
{
    // In the uploaded run35567.root file:
    //   BHC paddles use left/right readout names.
    //   BHD paddles use down/up readout names.
    // This macro follows the user's requested "left" side for BHC, and uses
    // the corresponding first BHD side, "down", because BHD has no "left" histograms.
    if (plane == "BHD") return "down";
    return "left";
}

std::vector<TString> BuildHistNameCandidatesMemSafe(const TString& plane, int paddle, const TString& type)
{
    std::vector<TString> names;

    TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    TString side = DefaultSideForPlaneMemSafe(plane);

    if (type == "QDC") {
        // Actual names seen in run35567.root:
        //   qdc_hit_BHC00left
        //   qdc_hit_BHD00down
        names.push_back(Form("qdc_hit_%s%s", paddleTag.Data(), side.Data()));

        // Fallbacks, useful if a later file uses the opposite SiPM side.
        if (plane == "BHC") names.push_back(Form("qdc_hit_%sright", paddleTag.Data()));
        if (plane == "BHD") names.push_back(Form("qdc_hit_%sup", paddleTag.Data()));
    }

    if (type == "TDC") {
        // Actual names seen in run35567.root use "tdc_coinc" with a c:
        //   tdc_coinc_left_BHC00-RF
        //   tdc_coinc_down_BHD00-RF
        names.push_back(Form("tdc_coinc_%s_%s-RF", side.Data(), paddleTag.Data()));

        // Fallback for the earlier typo/name assumption, "tdc_coin" without c.
        names.push_back(Form("tdc_coin_%s_%s-RF", side.Data(), paddleTag.Data()));

        // Fallbacks for the opposite SiPM side.
        if (plane == "BHC") {
            names.push_back(Form("tdc_coinc_right_%s-RF", paddleTag.Data()));
            names.push_back(Form("tdc_coin_right_%s-RF", paddleTag.Data()));
        }
        if (plane == "BHD") {
            names.push_back(Form("tdc_coinc_up_%s-RF", paddleTag.Data()));
            names.push_back(Form("tdc_coin_up_%s-RF", paddleTag.Data()));
        }
    }

    // Remove accidental duplicates while preserving order.
    std::vector<TString> uniqueNames;
    for (const auto& name : names) {
        bool seen = false;
        for (const auto& existing : uniqueNames) {
            if (existing == name) { seen = true; break; }
        }
        if (!seen && name.Length() > 0) uniqueNames.push_back(name);
    }

    return uniqueNames;
}

TString BuildHistNameMemSafe(const TString& plane, int paddle, const TString& type)
{
    std::vector<TString> names = BuildHistNameCandidatesMemSafe(plane, paddle, type);
    if (!names.empty()) return names.front();
    return "";
}

TH1* AsHistogramMemSafe(TObject* obj)
{
    if (obj == nullptr) {
        return nullptr;
    }
    if (obj->InheritsFrom(TH1::Class())) {
        return dynamic_cast<TH1*>(obj);
    }
    return nullptr;
}

bool KeyIsDirectoryMemSafe(TKey* key)
{
    if (key == nullptr) return false;
    TClass* cl = gROOT->GetClass(key->GetClassName());
    return (cl != nullptr && cl->InheritsFrom(TDirectory::Class()));
}

bool KeyIsHistogramMemSafe(TKey* key)
{
    if (key == nullptr) return false;
    TClass* cl = gROOT->GetClass(key->GetClassName());
    return (cl != nullptr && cl->InheritsFrom(TH1::Class()));
}

// This recursive search is much safer than reading every object.
// It reads directory objects to descend into them, but it only reads histogram
// objects when the key name already matches the requested histogram name.
TH1* FindHistRecursiveByKeyMemSafe(TDirectory* dir,
                                   const TString& wantedName,
                                   TString currentPath,
                                   TString& foundPath,
                                   int depth = 0,
                                   int maxDepth = 8)
{
    if (dir == nullptr || depth > maxDepth) {
        return nullptr;
    }

    TIter nextKey(dir->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = dynamic_cast<TKey*>(nextKey()))) {
        TString keyName = key->GetName();
        TString childPath = currentPath.Length() > 0
                            ? Form("%s/%s", currentPath.Data(), keyName.Data())
                            : keyName;

        if (keyName == wantedName && KeyIsHistogramMemSafe(key)) {
            TObject* obj = key->ReadObj();
            TH1* h = AsHistogramMemSafe(obj);
            if (h != nullptr) {
                foundPath = childPath;
                return h;
            }
            if (obj != nullptr) {
                delete obj;
            }
        }

        if (KeyIsDirectoryMemSafe(key)) {
            TObject* obj = key->ReadObj();
            TDirectory* subdir = dynamic_cast<TDirectory*>(obj);
            if (subdir == nullptr) {
                if (obj != nullptr) delete obj;
                continue;
            }

            TH1* h = FindHistRecursiveByKeyMemSafe(subdir,
                                                   wantedName,
                                                   childPath,
                                                   foundPath,
                                                   depth + 1,
                                                   maxDepth);
            if (h != nullptr) {
                return h;
            }

            // We only need this subdirectory temporarily if it did not contain the histogram.
            delete subdir;
        }
    }

    return nullptr;
}

TH1* GetPaddleHistogramMemSafe(TFile* file,
                               const TString& plane,
                               int paddle,
                               const TString& type,
                               TString& usedPath,
                               bool allowRecursiveSearch)
{
    if (file == nullptr || file->IsZombie()) {
        return nullptr;
    }

    TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    std::vector<TString> histNames = BuildHistNameCandidatesMemSafe(plane, paddle, type);

    for (const auto& histName : histNames) {
        std::vector<TString> candidatePaths;

        // Actual layout in run35567.root: BHC/Paddle00/QDC/<hist> and BHC/Paddle00/TDC/<hist>.
        candidatePaths.push_back(Form("%s/Paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));

        // Older/alternate guesses kept as fallbacks.
        candidatePaths.push_back(Form("%s/%s/%s/%s", plane.Data(), paddleTag.Data(), type.Data(), histName.Data()));
        candidatePaths.push_back(Form("%s/%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));
        candidatePaths.push_back(Form("%s/paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), histName.Data()));
        candidatePaths.push_back(Form("%s/%s", plane.Data(), histName.Data()));

        for (const auto& path : candidatePaths) {
            TObject* obj = file->Get(path.Data());
            TH1* h = AsHistogramMemSafe(obj);
            if (h != nullptr) {
                usedPath = path;
                return h;
            }
        }
    }

    if (!allowRecursiveSearch) {
        usedPath = "";
        return nullptr;
    }

    TDirectory* planeDir = dynamic_cast<TDirectory*>(file->Get(plane.Data()));
    if (planeDir != nullptr) {
        for (const auto& histName : histNames) {
            TH1* h = FindHistRecursiveByKeyMemSafe(planeDir, histName, plane, usedPath);
            if (h != nullptr) {
                return h;
            }
        }
    }

    for (const auto& histName : histNames) {
        TH1* h = FindHistRecursiveByKeyMemSafe(file, histName, "", usedPath);
        if (h != nullptr) {
            return h;
        }
    }

    usedPath = "";
    return nullptr;
}

void EstimateGaussianWindowMemSafe(TH1* h,
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

    if (!IsFiniteMemSafe(fwhm) || fwhm <= 2.0 * binWidth || fwhm > 0.8 * fullRange) {
        fwhm = 0.10 * fullRange;
    }

    sigmaGuess = fwhm / 2.355;
    if (!IsFiniteMemSafe(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = std::max(h->GetRMS(), 5.0 * binWidth);
    }
    if (!IsFiniteMemSafe(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = 0.05 * fullRange;
    }

    double window = std::max(3.0 * sigmaGuess, 6.0 * binWidth);
    window = std::min(window, 0.30 * fullRange);

    fitLow = std::max(xmin, peakX - window);
    fitHigh = std::min(xmax, peakX + window);

    if (fitHigh <= fitLow) {
        fitLow = xmin;
        fitHigh = xmax;
    }
}

void DrawFitTextMemSafe(const PeakFitResultMemSafe& r)
{
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextSize(0.032);

    TString peakText = FormatNumberHundredthMemSafe(r.peak);
    TString peakErrText = FormatNumberHundredthMemSafe(r.peakErr);
    TString sigmaText = FormatNumberHundredthMemSafe(r.sigma);
    TString sigmaErrText = FormatNumberHundredthMemSafe(r.sigmaErr);
    TString chi2Text = FormatNumberHundredthMemSafe(r.chi2);
    TString ndfText = FormatNumberHundredthMemSafe(r.ndf);

    double y = 0.86;
    latex.DrawLatex(0.58, y, Form("Peak = %s #pm %s", peakText.Data(), peakErrText.Data()));
    y -= 0.045;
    latex.DrawLatex(0.58, y, Form("#sigma = %s #pm %s", sigmaText.Data(), sigmaErrText.Data()));
    y -= 0.045;

    if (r.ndf > 0.0) {
        latex.DrawLatex(0.58, y, Form("#chi^{2}/ndf = %s / %s", chi2Text.Data(), ndfText.Data()));
    } else {
        latex.DrawLatex(0.58, y, "#chi^{2}/ndf unavailable");
    }
}

PeakFitResultMemSafe FitAndPlotHistogramMemSafe(TH1* h,
                                                int run,
                                                const TString& plane,
                                                int paddle,
                                                const TString& type,
                                                const TString& histPath,
                                                const TString& outputDir,
                                                const TString& pdfFile,
                                                TCanvas* c,
                                                bool savePlots)
{
    PeakFitResultMemSafe r;
    r.run = run;
    r.plane = plane;
    r.paddle = paddle;
    r.type = type;
    r.histName = BuildHistNameMemSafe(plane, paddle, type);
    r.histPath = histPath;
    r.found = (h != nullptr);

    TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    TString plotName = Form("run%d_%s_%s", run, paddleTag.Data(), type.Data());

    if (h == nullptr) {
        r.message = "histogram not found";

        if (savePlots && c != nullptr) {
            c->Clear();
            c->cd();
            TLatex latex;
            latex.SetNDC(kTRUE);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.12, 0.75, Form("Missing histogram: %s", r.histName.Data()));
            latex.DrawLatex(0.12, 0.68, Form("run %d, %s paddle %02d, %s", run, plane.Data(), paddle, type.Data()));

            TString png = Form("%s/%s_MISSING.png", outputDir.Data(), plotName.Data());
            c->SaveAs(png.Data());
            if (pdfFile.Length() > 0) c->Print(pdfFile.Data());
            c->Clear();
        }

        return r;
    }

    r.entries = h->GetEntries();

    if (r.entries < 10) {
        r.message = "too few entries";

        if (savePlots && c != nullptr) {
            c->Clear();
            c->cd();
            h->SetTitle(Form("run %d %s%02d %s: %s;%s channel;Counts",
                             run, plane.Data(), paddle, type.Data(), r.histName.Data(), type.Data()));
            h->Draw("hist");
            TLatex latex;
            latex.SetNDC(kTRUE);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.15, 0.82, "Too few entries for a reliable Gaussian fit.");
            TString png = Form("%s/%s_TOO_FEW_ENTRIES.png", outputDir.Data(), plotName.Data());
            c->SaveAs(png.Data());
            if (pdfFile.Length() > 0) c->Print(pdfFile.Data());
            c->Clear();
        }

        return r;
    }

    const int peakBin = h->GetMaximumBin();
    const double peakY = h->GetBinContent(peakBin);

    if (peakY <= 0.0) {
        r.message = "histogram maximum is non-positive";

        if (savePlots && c != nullptr) {
            c->Clear();
            c->cd();
            h->Draw("hist");
            TString png = Form("%s/%s_NO_POSITIVE_PEAK.png", outputDir.Data(), plotName.Data());
            c->SaveAs(png.Data());
            if (pdfFile.Length() > 0) c->Print(pdfFile.Data());
            c->Clear();
        }

        return r;
    }

    double peakGuess = 0.0;
    double sigmaGuess = 0.0;
    EstimateGaussianWindowMemSafe(h, peakGuess, sigmaGuess, r.fitLow, r.fitHigh);

    TString fitName = Form("gaus_%d_%s_%02d_%s", run, plane.Data(), paddle, type.Data());
    TF1* fGaus = new TF1(fitName.Data(), "gaus", r.fitLow, r.fitHigh);
    fGaus->SetParameters(peakY, peakGuess, sigmaGuess);
    fGaus->SetParNames("Amplitude", "Mean", "Sigma");

    const double binWidth = h->GetXaxis()->GetBinWidth(peakBin);
    const double fullRange = h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin();
    fGaus->SetParLimits(2, std::max(0.05 * binWidth, 1.0e-9), fullRange);

    // Q = quiet, R = fit range, 0 = do not draw during fit, N = do not store fit function.
    TFitResultPtr fitPtr = h->Fit(fGaus, "QR0N");
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
            IsFiniteMemSafe(r.peak) &&
            IsFiniteMemSafe(r.peakErr) &&
            r.peakErr > 0.0 &&
            IsFiniteMemSafe(r.sigma) &&
            r.sigma > 0.0);

    if (!r.ok) {
        r.message = Form("fit completed with status %d; inspect plot", r.fitStatus);
    } else {
        r.message = "ok";
    }

    if (savePlots && c != nullptr) {
        c->Clear();
        c->cd();
        c->SetGrid();

        h->SetLineWidth(2);
        h->SetTitle(Form("run %d %s%02d %s: %s;%s channel;Counts",
                         run, plane.Data(), paddle, type.Data(), r.histName.Data(), type.Data()));
        h->Draw("hist");

        fGaus->SetLineColor(kRed + 1);
        fGaus->SetLineWidth(3);
        fGaus->Draw("same");

        TLine lowLine(r.fitLow, 0.0, r.fitLow, peakY);
        TLine highLine(r.fitHigh, 0.0, r.fitHigh, peakY);
        lowLine.SetLineStyle(2);
        highLine.SetLineStyle(2);
        lowLine.SetLineColor(kGray + 2);
        highLine.SetLineColor(kGray + 2);
        lowLine.Draw("same");
        highLine.Draw("same");

        DrawFitTextMemSafe(r);

        TString statusLabel = r.ok ? "OK" : "CHECK";
        TString png = Form("%s/%s_%s.png", outputDir.Data(), plotName.Data(), statusLabel.Data());
        c->SaveAs(png.Data());
        if (pdfFile.Length() > 0) c->Print(pdfFile.Data());
        c->Clear();
    }

    delete fGaus;
    return r;
}

void WriteCsvRowMemSafe(std::ofstream& csv, const PeakFitResultMemSafe& r)
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
        << FormatNumberHundredthMemSafe(r.entries).Data() << ","
        << FormatNumberHundredthMemSafe(r.peak).Data() << ","
        << FormatNumberHundredthMemSafe(r.peakErr).Data() << ","
        << FormatNumberHundredthMemSafe(r.sigma).Data() << ","
        << FormatNumberHundredthMemSafe(r.sigmaErr).Data() << ","
        << FormatNumberHundredthMemSafe(r.amplitude).Data() << ","
        << FormatNumberHundredthMemSafe(r.amplitudeErr).Data() << ","
        << FormatNumberHundredthMemSafe(r.chi2).Data() << ","
        << FormatNumberHundredthMemSafe(r.ndf).Data() << ","
        << FormatNumberHundredthMemSafe(r.fitLow).Data() << ","
        << FormatNumberHundredthMemSafe(r.fitHigh).Data() << ","
        << r.message
        << "\n";
}

void BHPeakFits_QDC_TDC_memsafe_fixedpaths(const char* inputDir = ".",
                                const char* outputDir = "BH_peak_output",
                                int firstRun = 35566,
                                int lastRun = 35572,
                                bool savePlots = true,
                                bool allowRecursiveSearch = false)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatFormat("6.2f");
    gStyle->SetFitFormat("6.2f");

    // Avoid ROOT automatically attaching newly created histograms to directories.
    // This helps prevent hidden ownership/memory accumulation in long macros.
    TH1::AddDirectory(kFALSE);

    TString outDir(outputDir);
    MakeDirectoryIfNeededMemSafe(outDir);

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
        MakeDirectoryIfNeededMemSafe(runDir);

        for (const auto& plane : planes) {
            int nPaddles = (plane == "BHC") ? 16 : 13;

            TString planeDir = Form("%s/%s", runDir.Data(), plane.Data());
            MakeDirectoryIfNeededMemSafe(planeDir);

            for (const auto& type : types) {
                TString typeDir = Form("%s/%s", planeDir.Data(), type.Data());
                MakeDirectoryIfNeededMemSafe(typeDir);

                TString pdfFile = "";
                TCanvas* c = nullptr;

                if (savePlots) {
                    pdfFile = Form("%s/run%d_%s_%s_gaussian_fits.pdf",
                                   typeDir.Data(), run, plane.Data(), type.Data());
                    c = new TCanvas(Form("c_run%d_%s_%s", run, plane.Data(), type.Data()),
                                    Form("run %d %s %s Gaussian fits", run, plane.Data(), type.Data()),
                                    950, 700);
                    c->Print(Form("%s[", pdfFile.Data()));
                }

                for (int paddle = 0; paddle < nPaddles; paddle++) {
                    totalRequested++;

                    TString usedPath = "";
                    TH1* h = GetPaddleHistogramMemSafe(file, plane, paddle, type, usedPath, allowRecursiveSearch);

                    PeakFitResultMemSafe r = FitAndPlotHistogramMemSafe(h,
                                                                        run,
                                                                        plane,
                                                                        paddle,
                                                                        type,
                                                                        usedPath,
                                                                        typeDir,
                                                                        pdfFile,
                                                                        c,
                                                                        savePlots);

                    if (r.found) totalFound++;
                    if (r.ok) totalGoodFits++;

                    WriteCsvRowMemSafe(csv, r);

                    std::cout << "  run " << run
                              << " " << plane << Form("%02d", paddle)
                              << " " << type
                              << " : peak = " << FormatNumberHundredthMemSafe(r.peak).Data()
                              << " +/- " << FormatNumberHundredthMemSafe(r.peakErr).Data()
                              << "  [" << r.message << "]"
                              << std::endl;
                }

                if (savePlots && c != nullptr) {
                    c->Print(Form("%s]", pdfFile.Data()));
                    delete c;
                    c = nullptr;
                }
            }
        }

        file->Close();
        delete file;
        file = nullptr;

        // ROOT can hold onto deleted objects in internal cleanup lists.
        // This encourages cleanup between runs.
        gDirectory = nullptr;
    }

    csv.close();

    std::cout << "\nDone.\n"
              << "Requested histograms: " << totalRequested << "\n"
              << "Found histograms:     " << totalFound << "\n"
              << "Good fits:            " << totalGoodFits << "\n"
              << "CSV summary:          " << csvName << "\n"
              << "Plots directory:      " << outDir << "\n"
              << "savePlots:            " << (savePlots ? "true" : "false") << "\n"
              << "recursive search:     " << (allowRecursiveSearch ? "true" : "false") << "\n"
              << std::endl;
}

// Wrapper matching this file name, so ROOT can run it directly with:
//   root -l -b -q 'BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded.C(".", "BH_peak_output", 35566, 35572)'
void BHPeakFits(const char* inputDir = ".",
                                                   const char* outputDir = "BH_peak_output",
                                                   int firstRun = 35566,
                                                   int lastRun = 35572,
                                                   bool savePlots = true,
                                                   bool allowRecursiveSearch = false)
{
    BHPeakFits_QDC_TDC_memsafe_fixedpaths(inputDir, outputDir, firstRun, lastRun, savePlots, allowRecursiveSearch);
}

// Backward-compatible wrapper: this lets you call the fixed version using the old function name.
void BHPeakFits_QDC_TDC_memsafe(const char* inputDir = ".",
                                const char* outputDir = "BH_peak_output",
                                int firstRun = 35566,
                                int lastRun = 35572,
                                bool savePlots = true,
                                bool allowRecursiveSearch = false)
{
    BHPeakFits_QDC_TDC_memsafe_fixedpaths(inputDir, outputDir, firstRun, lastRun, savePlots, allowRecursiveSearch);
}
