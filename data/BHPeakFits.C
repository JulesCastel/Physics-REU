// BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded.C
//
// Memory-safe Gaussian fitting macro for BHC/BHD QDC and TDC histograms.
//
// This revision preserves the ORIGINAL CSV columns and appends a complete
// two-peak analysis for every QDC histogram:
//
//   pedestal peak = left QDC peak
//   signal peak   = right QDC / high-voltage-dependent peak
//   separation    = signal_peak - pedestal_peak
//
// The uncertainty on the separation is propagated in quadrature:
//
//   separation_uncertainty = sqrt(signal_peak_uncertainty^2
//                               + pedestal_peak_uncertainty^2)
//
// The legacy columns "peak", "peak_uncertainty", etc. are intentionally
// retained with their original meaning: a Gaussian fit around the histogram's
// largest bin. This keeps old scripts compatible. The appended, explicitly
// named pedestal/signal columns should be used for the QDC high-voltage study.
//
// TDC histograms still receive the original single-Gaussian analysis. Their
// appended two-peak fields are written as nan / not applicable.
//
// Default run range: 35566 through 35572.
//
// Run from a shell:
//   root -l -b -q 'BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded.C(".", "BH_peak_output", 35566, 35572)'
//
// CSV output:
//   BH_peak_output/BH_QDC_TDC_peak_summary_runs35566_35572.csv

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
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace BHPeakFitsTwoPeakDetail {

struct PeakCandidate {
    int bin = -1;
    double x = std::numeric_limits<double>::quiet_NaN();
    double height = 0.0;
};

struct GaussianFit {
    bool ok = false;
    int status = -999;

    double amplitude = std::numeric_limits<double>::quiet_NaN();
    double amplitudeError = std::numeric_limits<double>::quiet_NaN();
    double mean = std::numeric_limits<double>::quiet_NaN();
    double meanError = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigmaError = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double ndf = std::numeric_limits<double>::quiet_NaN();
    double fitLow = std::numeric_limits<double>::quiet_NaN();
    double fitHigh = std::numeric_limits<double>::quiet_NaN();

    TString message = "";
};

struct HistogramResult {
    int run = -1;
    TString plane = "";
    int paddle = -1;
    TString type = "";
    TString histogram = "";
    TString histogramPath = "";

    bool found = false;
    double entries = 0.0;

    // Legacy fit around the histogram's largest bin.
    GaussianFit legacy;

    // QDC-only two-peak analysis.
    bool twoPeaksFound = false;
    bool twoPeakFitOk = false;
    int valleyBin = -1;
    double valleyX = std::numeric_limits<double>::quiet_NaN();
    GaussianFit pedestal;
    GaussianFit signal;
    double separation = std::numeric_limits<double>::quiet_NaN();
    double separationError = std::numeric_limits<double>::quiet_NaN();
    TString twoPeakMessage = "";
};

bool IsFinite(double x)
{
    return std::isfinite(x);
}

TString FormatHundredth(double x)
{
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return (x > 0.0) ? "inf" : "-inf";

    double rounded = std::round(x * 100.0) / 100.0;
    if (std::fabs(rounded) < 0.005) rounded = 0.0;

    const double nearestInteger = std::round(rounded);
    if (std::fabs(rounded - nearestInteger) < 1.0e-9) {
        return Form("%.0f", rounded);
    }
    return Form("%.2f", rounded);
}

TString CsvSafe(TString text)
{
    text.ReplaceAll(",", ";");
    text.ReplaceAll("\n", " ");
    text.ReplaceAll("\r", " ");
    return text;
}

void MakeDirectory(const TString& path)
{
    if (path.Length() > 0) gSystem->mkdir(path.Data(), kTRUE);
}

TString DefaultSide(const TString& plane)
{
    return (plane == "BHD") ? "down" : "left";
}

std::vector<TString> HistogramNameCandidates(const TString& plane,
                                             int paddle,
                                             const TString& type)
{
    const TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    const TString side = DefaultSide(plane);
    std::vector<TString> names;

    if (type == "QDC") {
        names.push_back(Form("qdc_hit_%s%s", paddleTag.Data(), side.Data()));
        if (plane == "BHC") names.push_back(Form("qdc_hit_%sright", paddleTag.Data()));
        if (plane == "BHD") names.push_back(Form("qdc_hit_%sup", paddleTag.Data()));
    } else if (type == "TDC") {
        names.push_back(Form("tdc_coinc_%s_%s-RF", side.Data(), paddleTag.Data()));
        names.push_back(Form("tdc_coin_%s_%s-RF", side.Data(), paddleTag.Data()));

        if (plane == "BHC") {
            names.push_back(Form("tdc_coinc_right_%s-RF", paddleTag.Data()));
            names.push_back(Form("tdc_coin_right_%s-RF", paddleTag.Data()));
        } else if (plane == "BHD") {
            names.push_back(Form("tdc_coinc_up_%s-RF", paddleTag.Data()));
            names.push_back(Form("tdc_coin_up_%s-RF", paddleTag.Data()));
        }
    }

    std::vector<TString> unique;
    for (const TString& name : names) {
        bool seen = false;
        for (const TString& existing : unique) {
            if (existing == name) {
                seen = true;
                break;
            }
        }
        if (!seen && name.Length() > 0) unique.push_back(name);
    }
    return unique;
}

bool KeyIsDirectory(TKey* key)
{
    if (key == nullptr) return false;
    TClass* objectClass = gROOT->GetClass(key->GetClassName());
    return objectClass != nullptr && objectClass->InheritsFrom(TDirectory::Class());
}

bool KeyIsHistogram(TKey* key)
{
    if (key == nullptr) return false;
    TClass* objectClass = gROOT->GetClass(key->GetClassName());
    return objectClass != nullptr && objectClass->InheritsFrom(TH1::Class());
}

TH1* CloneHistogram(TH1* source, const TString& cloneName)
{
    if (source == nullptr) return nullptr;
    TH1* clone = dynamic_cast<TH1*>(source->Clone(cloneName.Data()));
    if (clone != nullptr) clone->SetDirectory(nullptr);
    return clone;
}

TH1* FindHistogramRecursiveClone(TDirectory* directory,
                                 const TString& wantedName,
                                 const TString& currentPath,
                                 TString& foundPath,
                                 const TString& cloneName,
                                 int depth = 0,
                                 int maximumDepth = 8)
{
    if (directory == nullptr || depth > maximumDepth) return nullptr;

    TIter nextKey(directory->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = dynamic_cast<TKey*>(nextKey()))) {
        const TString keyName = key->GetName();
        const TString childPath = currentPath.Length() > 0
                                ? Form("%s/%s", currentPath.Data(), keyName.Data())
                                : keyName;

        if (keyName == wantedName && KeyIsHistogram(key)) {
            TObject* object = key->ReadObj();
            TH1* histogram = dynamic_cast<TH1*>(object);
            TH1* clone = CloneHistogram(histogram, cloneName);
            delete object;

            if (clone != nullptr) {
                foundPath = childPath;
                return clone;
            }
        }

        if (KeyIsDirectory(key)) {
            TObject* object = key->ReadObj();
            TDirectory* subdirectory = dynamic_cast<TDirectory*>(object);
            if (subdirectory == nullptr) {
                delete object;
                continue;
            }

            TH1* clone = FindHistogramRecursiveClone(subdirectory,
                                                     wantedName,
                                                     childPath,
                                                     foundPath,
                                                     cloneName,
                                                     depth + 1,
                                                     maximumDepth);
            delete subdirectory;
            if (clone != nullptr) return clone;
        }
    }

    return nullptr;
}

std::unique_ptr<TH1> GetHistogramClone(TFile* file,
                                       int run,
                                       const TString& plane,
                                       int paddle,
                                       const TString& type,
                                       TString& usedName,
                                       TString& usedPath,
                                       bool allowRecursiveSearch)
{
    usedName = "";
    usedPath = "";
    if (file == nullptr || file->IsZombie()) return nullptr;

    const std::vector<TString> names = HistogramNameCandidates(plane, paddle, type);
    const TString cloneName = Form("work_run%d_%s%02d_%s",
                                   run, plane.Data(), paddle, type.Data());

    for (const TString& name : names) {
        const TString paddleTag = Form("%s%02d", plane.Data(), paddle);
        const std::vector<TString> paths = {
            Form("%s/Paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), name.Data()),
            Form("%s/%s/%s/%s", plane.Data(), paddleTag.Data(), type.Data(), name.Data()),
            Form("%s/%02d/%s/%s", plane.Data(), paddle, type.Data(), name.Data()),
            Form("%s/paddle%02d/%s/%s", plane.Data(), paddle, type.Data(), name.Data()),
            Form("%s/%s", plane.Data(), name.Data())
        };

        for (const TString& path : paths) {
            TObject* object = file->Get(path.Data());
            TH1* histogram = dynamic_cast<TH1*>(object);
            if (histogram == nullptr) continue;

            TH1* clone = CloneHistogram(histogram, cloneName);
            if (clone != nullptr) {
                usedName = name;
                usedPath = path;
                return std::unique_ptr<TH1>(clone);
            }
        }
    }

    if (!allowRecursiveSearch) return nullptr;

    TDirectory* planeDirectory = dynamic_cast<TDirectory*>(file->Get(plane.Data()));
    if (planeDirectory != nullptr) {
        for (const TString& name : names) {
            TH1* clone = FindHistogramRecursiveClone(planeDirectory,
                                                     name,
                                                     plane,
                                                     usedPath,
                                                     cloneName);
            if (clone != nullptr) {
                usedName = name;
                return std::unique_ptr<TH1>(clone);
            }
        }
    }

    for (const TString& name : names) {
        TH1* clone = FindHistogramRecursiveClone(file,
                                                 name,
                                                 "",
                                                 usedPath,
                                                 cloneName);
        if (clone != nullptr) {
            usedName = name;
            return std::unique_ptr<TH1>(clone);
        }
    }

    return nullptr;
}

void EstimateLegacyWindow(TH1* histogram,
                          double& peakGuess,
                          double& sigmaGuess,
                          double& fitLow,
                          double& fitHigh)
{
    const int numberOfBins = histogram->GetNbinsX();
    const int peakBin = histogram->GetMaximumBin();
    const double xMinimum = histogram->GetXaxis()->GetXmin();
    const double xMaximum = histogram->GetXaxis()->GetXmax();
    const double fullRange = xMaximum - xMinimum;
    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);

    peakGuess = histogram->GetBinCenter(peakBin);
    const double halfMaximum = 0.5 * histogram->GetBinContent(peakBin);

    int leftBin = peakBin;
    int rightBin = peakBin;
    while (leftBin > 1 && histogram->GetBinContent(leftBin) > halfMaximum) --leftBin;
    while (rightBin < numberOfBins
           && histogram->GetBinContent(rightBin) > halfMaximum) ++rightBin;

    double fwhm = histogram->GetBinLowEdge(rightBin + 1)
                - histogram->GetBinLowEdge(leftBin);
    if (!IsFinite(fwhm) || fwhm <= 2.0 * binWidth || fwhm > 0.8 * fullRange) {
        fwhm = 0.10 * fullRange;
    }

    sigmaGuess = fwhm / 2.355;
    if (!IsFinite(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = std::max(histogram->GetRMS(), 5.0 * binWidth);
    }
    if (!IsFinite(sigmaGuess) || sigmaGuess <= 0.0) {
        sigmaGuess = 0.05 * fullRange;
    }

    double window = std::max(3.0 * sigmaGuess, 6.0 * binWidth);
    window = std::min(window, 0.30 * fullRange);
    fitLow = std::max(xMinimum, peakGuess - window);
    fitHigh = std::min(xMaximum, peakGuess + window);

    if (!(fitHigh > fitLow)) {
        fitLow = xMinimum;
        fitHigh = xMaximum;
    }
}

GaussianFit FitLegacyDominantPeak(TH1* histogram,
                                  int run,
                                  const TString& plane,
                                  int paddle,
                                  const TString& type)
{
    GaussianFit result;
    if (histogram == nullptr || histogram->GetEntries() < 10) {
        result.message = "too few entries or missing histogram";
        return result;
    }

    const int peakBin = histogram->GetMaximumBin();
    const double peakHeight = histogram->GetBinContent(peakBin);
    if (!(peakHeight > 0.0)) {
        result.message = "histogram maximum is non-positive";
        return result;
    }

    double peakGuess = 0.0;
    double sigmaGuess = 0.0;
    EstimateLegacyWindow(histogram,
                         peakGuess,
                         sigmaGuess,
                         result.fitLow,
                         result.fitHigh);

    const TString name = Form("legacy_gaus_run%d_%s%02d_%s",
                              run, plane.Data(), paddle, type.Data());
    TF1 gaussian(name.Data(), "gaus", result.fitLow, result.fitHigh);
    gaussian.SetParNames("Amplitude", "Mean", "Sigma");
    gaussian.SetParameters(peakHeight, peakGuess, sigmaGuess);

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    const double fullRange = histogram->GetXaxis()->GetXmax()
                           - histogram->GetXaxis()->GetXmin();
    gaussian.SetParLimits(0, 0.0, std::max(10.0 * peakHeight, 1.0));
    gaussian.SetParLimits(1, result.fitLow, result.fitHigh);
    gaussian.SetParLimits(2, std::max(0.05 * binWidth, 1.0e-9), fullRange);

    TFitResultPtr fit = histogram->Fit(&gaussian, "QR0SN");
    result.status = static_cast<int>(fit);
    result.amplitude = gaussian.GetParameter(0);
    result.amplitudeError = gaussian.GetParError(0);
    result.mean = gaussian.GetParameter(1);
    result.meanError = gaussian.GetParError(1);
    result.sigma = std::fabs(gaussian.GetParameter(2));
    result.sigmaError = gaussian.GetParError(2);
    result.chi2 = gaussian.GetChisquare();
    result.ndf = gaussian.GetNDF();

    result.ok = result.status == 0
             && IsFinite(result.mean)
             && IsFinite(result.meanError)
             && result.meanError > 0.0
             && IsFinite(result.sigma)
             && result.sigma > 0.0;
    result.message = result.ok ? "ok" : Form("fit status %d", result.status);
    return result;
}

std::vector<PeakCandidate> FindLocalMaxima(TH1* histogram,
                                           double relativeThreshold)
{
    std::vector<PeakCandidate> candidates;
    if (histogram == nullptr || histogram->GetNbinsX() < 7) return candidates;

    std::unique_ptr<TH1> smoothed(CloneHistogram(
        histogram, Form("%s_peak_search", histogram->GetName())));
    if (!smoothed) return candidates;
    smoothed->Smooth(3);

    const double maximum = smoothed->GetMaximum();
    if (!(maximum > 0.0)) return candidates;
    const double threshold = std::max(0.0, relativeThreshold) * maximum;

    for (int bin = 3; bin <= smoothed->GetNbinsX() - 2; ++bin) {
        const double center = smoothed->GetBinContent(bin);
        if (center < threshold) continue;
        if (center >= smoothed->GetBinContent(bin - 1)
            && center > smoothed->GetBinContent(bin + 1)) {
            PeakCandidate candidate;
            candidate.bin = bin;
            candidate.x = smoothed->GetBinCenter(bin);
            candidate.height = center;
            candidates.push_back(candidate);
        }
    }
    return candidates;
}

bool SelectTwoSeparatedPeaks(TH1* histogram,
                             PeakCandidate& leftPeak,
                             PeakCandidate& rightPeak)
{
    std::vector<PeakCandidate> candidates = FindLocalMaxima(histogram, 0.03);
    if (candidates.size() < 2) candidates = FindLocalMaxima(histogram, 0.01);
    if (candidates.size() < 2) return false;

    std::sort(candidates.begin(), candidates.end(),
              [](const PeakCandidate& a, const PeakCandidate& b) {
                  return a.height > b.height;
              });

    const double xRange = histogram->GetXaxis()->GetXmax()
                        - histogram->GetXaxis()->GetXmin();
    const double binWidth = histogram->GetXaxis()->GetBinWidth(1);
    const double minimumSeparation = std::max(6.0 * binWidth, 0.03 * xRange);
    const std::size_t maximumCandidates = std::min<std::size_t>(12, candidates.size());

    bool found = false;
    double bestScore = -1.0;
    PeakCandidate bestA;
    PeakCandidate bestB;

    for (std::size_t i = 0; i < maximumCandidates; ++i) {
        for (std::size_t j = i + 1; j < maximumCandidates; ++j) {
            const double separation = std::fabs(candidates[i].x - candidates[j].x);
            if (separation < minimumSeparation) continue;

            const double score = std::min(candidates[i].height, candidates[j].height)
                               + 0.05 * (candidates[i].height + candidates[j].height)
                               + 0.001 * separation;
            if (score > bestScore) {
                bestScore = score;
                bestA = candidates[i];
                bestB = candidates[j];
                found = true;
            }
        }
    }

    if (!found) return false;
    if (bestA.x < bestB.x) {
        leftPeak = bestA;
        rightPeak = bestB;
    } else {
        leftPeak = bestB;
        rightPeak = bestA;
    }
    return true;
}

int FindValleyBin(TH1* histogram, int leftPeakBin, int rightPeakBin)
{
    if (histogram == nullptr || leftPeakBin >= rightPeakBin) return -1;

    std::unique_ptr<TH1> smoothed(CloneHistogram(
        histogram, Form("%s_valley_search", histogram->GetName())));
    if (!smoothed) return -1;
    smoothed->Smooth(2);

    int valleyBin = leftPeakBin + 1;
    double minimum = smoothed->GetBinContent(valleyBin);
    for (int bin = leftPeakBin + 1; bin < rightPeakBin; ++bin) {
        const double value = smoothed->GetBinContent(bin);
        if (value < minimum) {
            minimum = value;
            valleyBin = bin;
        }
    }
    return valleyBin;
}

void EstimateBoundedWindow(TH1* histogram,
                           int peakBin,
                           int lowerBoundaryBin,
                           int upperBoundaryBin,
                           double& fitLow,
                           double& fitHigh,
                           double& sigmaGuess)
{
    const int numberOfBins = histogram->GetNbinsX();
    lowerBoundaryBin = std::max(1, lowerBoundaryBin);
    upperBoundaryBin = std::min(numberOfBins, upperBoundaryBin);

    const double peakHeight = histogram->GetBinContent(peakBin);
    const double halfHeight = 0.5 * peakHeight;
    int halfLeft = peakBin;
    int halfRight = peakBin;

    while (halfLeft > lowerBoundaryBin
           && histogram->GetBinContent(halfLeft) > halfHeight) --halfLeft;
    while (halfRight < upperBoundaryBin
           && histogram->GetBinContent(halfRight) > halfHeight) ++halfRight;

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    double fwhm = histogram->GetBinLowEdge(halfRight + 1)
                - histogram->GetBinLowEdge(halfLeft);
    if (!IsFinite(fwhm) || fwhm < 2.0 * binWidth) fwhm = 6.0 * binWidth;

    sigmaGuess = std::max(fwhm / 2.355, binWidth);
    const double peakX = histogram->GetBinCenter(peakBin);
    const double lowerX = histogram->GetBinLowEdge(lowerBoundaryBin);
    const double upperX = histogram->GetBinLowEdge(upperBoundaryBin + 1);

    fitLow = std::max(lowerX, peakX - 2.8 * sigmaGuess);
    fitHigh = std::min(upperX, peakX + 2.8 * sigmaGuess);

    if (fitHigh - fitLow < 5.0 * binWidth) {
        fitLow = std::max(lowerX, peakX - 3.0 * binWidth);
        fitHigh = std::min(upperX, peakX + 3.0 * binWidth);
    }
}

GaussianFit FitBoundedPeak(TH1* histogram,
                           int run,
                           const TString& plane,
                           int paddle,
                           const TString& label,
                           int peakBin,
                           int lowerBoundaryBin,
                           int upperBoundaryBin)
{
    GaussianFit result;
    if (histogram == nullptr || peakBin < 1 || peakBin > histogram->GetNbinsX()) {
        result.message = "invalid histogram or peak bin";
        return result;
    }

    double sigmaGuess = 0.0;
    EstimateBoundedWindow(histogram,
                          peakBin,
                          lowerBoundaryBin,
                          upperBoundaryBin,
                          result.fitLow,
                          result.fitHigh,
                          sigmaGuess);
    if (!(result.fitHigh > result.fitLow)) {
        result.message = "invalid fit range";
        return result;
    }

    const double peakX = histogram->GetBinCenter(peakBin);
    const double peakHeight = histogram->GetBinContent(peakBin);
    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    const double fullRange = histogram->GetXaxis()->GetXmax()
                           - histogram->GetXaxis()->GetXmin();

    const TString name = Form("%s_gaus_run%d_%s%02d",
                              label.Data(), run, plane.Data(), paddle);
    TF1 gaussian(name.Data(), "gaus", result.fitLow, result.fitHigh);
    gaussian.SetParNames("Amplitude", "Mean", "Sigma");
    gaussian.SetParameters(peakHeight, peakX, sigmaGuess);
    gaussian.SetParLimits(0, 0.0, std::max(10.0 * peakHeight, 1.0));
    gaussian.SetParLimits(1, result.fitLow, result.fitHigh);
    gaussian.SetParLimits(2, std::max(0.25 * binWidth, 1.0e-9), fullRange);

    TFitResultPtr fit = histogram->Fit(&gaussian, "QR0SN");
    result.status = static_cast<int>(fit);
    result.amplitude = gaussian.GetParameter(0);
    result.amplitudeError = gaussian.GetParError(0);
    result.mean = gaussian.GetParameter(1);
    result.meanError = gaussian.GetParError(1);
    result.sigma = std::fabs(gaussian.GetParameter(2));
    result.sigmaError = gaussian.GetParError(2);
    result.chi2 = gaussian.GetChisquare();
    result.ndf = gaussian.GetNDF();

    result.ok = result.status == 0
             && IsFinite(result.mean)
             && IsFinite(result.meanError)
             && result.meanError > 0.0
             && IsFinite(result.sigma)
             && result.sigma > 0.0
             && result.mean >= result.fitLow
             && result.mean <= result.fitHigh;
    result.message = result.ok ? "ok" : Form("fit status %d", result.status);
    return result;
}

HistogramResult AnalyzeHistogram(TH1* histogram,
                                 int run,
                                 const TString& plane,
                                 int paddle,
                                 const TString& type,
                                 const TString& histogramName,
                                 const TString& histogramPath)
{
    HistogramResult result;
    result.run = run;
    result.plane = plane;
    result.paddle = paddle;
    result.type = type;
    result.histogram = histogramName;
    result.histogramPath = histogramPath;
    result.found = histogram != nullptr;

    if (histogram == nullptr) {
        result.legacy.message = "histogram not found";
        result.twoPeakMessage = (type == "QDC")
                              ? "histogram not found"
                              : "not applicable for TDC";
        return result;
    }

    result.entries = histogram->GetEntries();
    result.legacy = FitLegacyDominantPeak(histogram, run, plane, paddle, type);

    if (type != "QDC") {
        result.twoPeakMessage = "not applicable for TDC";
        return result;
    }

    if (result.entries < 20) {
        result.twoPeakMessage = "too few entries for two-peak fit";
        return result;
    }

    PeakCandidate leftPeak;
    PeakCandidate rightPeak;
    result.twoPeaksFound = SelectTwoSeparatedPeaks(histogram, leftPeak, rightPeak);
    if (!result.twoPeaksFound) {
        result.twoPeakMessage = "could not identify two separated peaks";
        return result;
    }

    result.valleyBin = FindValleyBin(histogram, leftPeak.bin, rightPeak.bin);
    if (result.valleyBin <= leftPeak.bin || result.valleyBin >= rightPeak.bin) {
        result.twoPeakMessage = "could not identify valley between peaks";
        return result;
    }
    result.valleyX = histogram->GetBinCenter(result.valleyBin);

    result.pedestal = FitBoundedPeak(histogram,
                                     run,
                                     plane,
                                     paddle,
                                     "pedestal",
                                     leftPeak.bin,
                                     1,
                                     result.valleyBin);
    result.signal = FitBoundedPeak(histogram,
                                   run,
                                   plane,
                                   paddle,
                                   "signal",
                                   rightPeak.bin,
                                   result.valleyBin,
                                   histogram->GetNbinsX());

    result.twoPeakFitOk = result.pedestal.ok
                       && result.signal.ok
                       && result.signal.mean > result.pedestal.mean;

    if (!result.twoPeakFitOk) {
        result.twoPeakMessage = Form("pedestal: %s; signal: %s",
                                     result.pedestal.message.Data(),
                                     result.signal.message.Data());
        return result;
    }

    result.separation = result.signal.mean - result.pedestal.mean;
    result.separationError = std::hypot(result.signal.meanError,
                                        result.pedestal.meanError);
    result.twoPeakMessage = "ok";
    return result;
}

void DrawResult(TCanvas* canvas,
                TH1* histogram,
                const HistogramResult& result,
                const TString& pngPath,
                const TString& pdfPath)
{
    if (canvas == nullptr) return;
    canvas->Clear();
    canvas->cd();
    canvas->SetGrid();

    TLatex text;
    text.SetNDC(kTRUE);
    text.SetTextSize(0.030);

    if (histogram == nullptr) {
        text.DrawLatex(0.12, 0.78,
                       Form("Missing histogram: run %d %s%02d %s",
                            result.run,
                            result.plane.Data(),
                            result.paddle,
                            result.type.Data()));
        canvas->SaveAs(pngPath.Data());
        if (pdfPath.Length() > 0) canvas->Print(pdfPath.Data());
        return;
    }

    histogram->SetLineWidth(2);
    histogram->SetTitle(Form("run %d %s%02d %s: %s;%s channel;Counts",
                             result.run,
                             result.plane.Data(),
                             result.paddle,
                             result.type.Data(),
                             result.histogram.Data(),
                             result.type.Data()));
    histogram->Draw("hist");

    std::vector<std::unique_ptr<TF1>> functions;
    std::unique_ptr<TLine> valleyLine;

    if (result.type == "QDC") {
        if (result.pedestal.ok) {
            functions.emplace_back(new TF1(
                Form("draw_ped_run%d_%s%02d", result.run, result.plane.Data(), result.paddle),
                "gaus", result.pedestal.fitLow, result.pedestal.fitHigh));
            functions.back()->SetParameters(result.pedestal.amplitude,
                                             result.pedestal.mean,
                                             result.pedestal.sigma);
            functions.back()->SetLineColor(kBlue + 1);
            functions.back()->SetLineWidth(3);
            functions.back()->Draw("same");
        }
        if (result.signal.ok) {
            functions.emplace_back(new TF1(
                Form("draw_sig_run%d_%s%02d", result.run, result.plane.Data(), result.paddle),
                "gaus", result.signal.fitLow, result.signal.fitHigh));
            functions.back()->SetParameters(result.signal.amplitude,
                                             result.signal.mean,
                                             result.signal.sigma);
            functions.back()->SetLineColor(kRed + 1);
            functions.back()->SetLineWidth(3);
            functions.back()->Draw("same");
        }

        if (IsFinite(result.valleyX)) {
            valleyLine.reset(new TLine(result.valleyX,
                                       0.0,
                                       result.valleyX,
                                       histogram->GetMaximum()));
            valleyLine->SetLineColor(kGray + 2);
            valleyLine->SetLineStyle(2);
            valleyLine->Draw("same");
        }

        double y = 0.87;
        if (result.pedestal.ok) {
            text.DrawLatex(0.56, y,
                           Form("Pedestal #mu = %s #pm %s",
                                FormatHundredth(result.pedestal.mean).Data(),
                                FormatHundredth(result.pedestal.meanError).Data()));
            y -= 0.045;
        }
        if (result.signal.ok) {
            text.DrawLatex(0.56, y,
                           Form("Signal #mu = %s #pm %s",
                                FormatHundredth(result.signal.mean).Data(),
                                FormatHundredth(result.signal.meanError).Data()));
            y -= 0.045;
        }
        if (result.twoPeakFitOk) {
            text.DrawLatex(0.56, y,
                           Form("Separation = %s #pm %s",
                                FormatHundredth(result.separation).Data(),
                                FormatHundredth(result.separationError).Data()));
        } else {
            text.SetTextColor(kRed + 1);
            text.DrawLatex(0.12, 0.92,
                           Form("TWO-PEAK CHECK: %s", result.twoPeakMessage.Data()));
        }
    } else {
        if (result.legacy.ok) {
            functions.emplace_back(new TF1(
                Form("draw_legacy_run%d_%s%02d_%s",
                     result.run,
                     result.plane.Data(),
                     result.paddle,
                     result.type.Data()),
                "gaus", result.legacy.fitLow, result.legacy.fitHigh));
            functions.back()->SetParameters(result.legacy.amplitude,
                                             result.legacy.mean,
                                             result.legacy.sigma);
            functions.back()->SetLineColor(kRed + 1);
            functions.back()->SetLineWidth(3);
            functions.back()->Draw("same");
        }

        text.DrawLatex(0.58, 0.86,
                       Form("Peak = %s #pm %s",
                            FormatHundredth(result.legacy.mean).Data(),
                            FormatHundredth(result.legacy.meanError).Data()));
        text.DrawLatex(0.58, 0.815,
                       Form("#sigma = %s #pm %s",
                            FormatHundredth(result.legacy.sigma).Data(),
                            FormatHundredth(result.legacy.sigmaError).Data()));
    }

    canvas->SaveAs(pngPath.Data());
    if (pdfPath.Length() > 0) canvas->Print(pdfPath.Data());
    canvas->Clear();
}

void WriteCsvHeader(std::ofstream& csv)
{
    // First 21 columns are exactly the original CSV schema.
    csv << "run,plane,paddle,type,histogram,hist_path,found,fit_ok,fit_status,"
        << "entries,peak,peak_uncertainty,sigma,sigma_uncertainty,"
        << "amplitude,amplitude_uncertainty,chi2,ndf,fit_low,fit_high,message,"

        // Appended QDC two-peak fields.
        << "two_peaks_found,two_peak_fit_ok,"
        << "pedestal_peak,pedestal_peak_uncertainty,pedestal_sigma,"
        << "pedestal_sigma_uncertainty,pedestal_amplitude,"
        << "pedestal_amplitude_uncertainty,pedestal_chi2,pedestal_ndf,"
        << "pedestal_fit_low,pedestal_fit_high,pedestal_fit_status,"
        << "signal_peak,signal_peak_uncertainty,signal_sigma,"
        << "signal_sigma_uncertainty,signal_amplitude,"
        << "signal_amplitude_uncertainty,signal_chi2,signal_ndf,"
        << "signal_fit_low,signal_fit_high,signal_fit_status,"
        << "peak_separation,peak_separation_uncertainty,two_peak_message\n";
}

void WriteFitFields(std::ofstream& csv, const GaussianFit& fit)
{
    csv << FormatHundredth(fit.mean) << ","
        << FormatHundredth(fit.meanError) << ","
        << FormatHundredth(fit.sigma) << ","
        << FormatHundredth(fit.sigmaError) << ","
        << FormatHundredth(fit.amplitude) << ","
        << FormatHundredth(fit.amplitudeError) << ","
        << FormatHundredth(fit.chi2) << ","
        << FormatHundredth(fit.ndf) << ","
        << FormatHundredth(fit.fitLow) << ","
        << FormatHundredth(fit.fitHigh) << ","
        << fit.status;
}

void WriteCsvRow(std::ofstream& csv, const HistogramResult& result)
{
    csv << result.run << ","
        << result.plane << ","
        << std::setw(2) << std::setfill('0') << result.paddle << std::setfill(' ') << ","
        << result.type << ","
        << CsvSafe(result.histogram) << ","
        << CsvSafe(result.histogramPath) << ","
        << (result.found ? 1 : 0) << ","
        << (result.legacy.ok ? 1 : 0) << ","
        << result.legacy.status << ","
        << FormatHundredth(result.entries) << ","
        << FormatHundredth(result.legacy.mean) << ","
        << FormatHundredth(result.legacy.meanError) << ","
        << FormatHundredth(result.legacy.sigma) << ","
        << FormatHundredth(result.legacy.sigmaError) << ","
        << FormatHundredth(result.legacy.amplitude) << ","
        << FormatHundredth(result.legacy.amplitudeError) << ","
        << FormatHundredth(result.legacy.chi2) << ","
        << FormatHundredth(result.legacy.ndf) << ","
        << FormatHundredth(result.legacy.fitLow) << ","
        << FormatHundredth(result.legacy.fitHigh) << ","
        << CsvSafe(result.legacy.message) << ","
        << (result.twoPeaksFound ? 1 : 0) << ","
        << (result.twoPeakFitOk ? 1 : 0) << ",";

    WriteFitFields(csv, result.pedestal);
    csv << ",";
    WriteFitFields(csv, result.signal);
    csv << ","
        << FormatHundredth(result.separation) << ","
        << FormatHundredth(result.separationError) << ","
        << CsvSafe(result.twoPeakMessage)
        << "\n";
}

} // namespace BHPeakFitsTwoPeakDetail

void BHPeakFits_QDC_TDC_memsafe_fixedpaths(const char* inputDir = ".",
                                            const char* outputDir = "BH_peak_output",
                                            int firstRun = 35566,
                                            int lastRun = 35572,
                                            bool savePlots = true,
                                            bool allowRecursiveSearch = false)
{
    using namespace BHPeakFitsTwoPeakDetail;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatFormat("6.2f");
    gStyle->SetFitFormat("6.2f");
    TH1::AddDirectory(kFALSE);

    const TString outputDirectory(outputDir);
    MakeDirectory(outputDirectory);

    const TString csvName = Form("%s/BH_QDC_TDC_peak_summary_runs%d_%d.csv",
                                 outputDirectory.Data(), firstRun, lastRun);
    std::ofstream csv(csvName.Data());
    if (!csv.is_open()) {
        std::cerr << "ERROR: could not create " << csvName << std::endl;
        return;
    }
    WriteCsvHeader(csv);

    const std::vector<TString> planes = {"BHC", "BHD"};
    const std::vector<TString> types = {"QDC", "TDC"};

    int requested = 0;
    int found = 0;
    int legacyFitsOk = 0;
    int qdcTwoPeakFitsOk = 0;

    for (int run = firstRun; run <= lastRun; ++run) {
        const TString fileName = Form("%s/run%d.root", inputDir, run);
        std::unique_ptr<TFile> file(TFile::Open(fileName.Data(), "READ"));

        if (!file || file->IsZombie()) {
            std::cerr << "WARNING: could not open " << fileName << std::endl;
            continue;
        }

        std::cout << "\nProcessing " << fileName << std::endl;
        const TString runDirectory = Form("%s/run%d", outputDirectory.Data(), run);
        MakeDirectory(runDirectory);

        for (const TString& plane : planes) {
            const int numberOfPaddles = (plane == "BHC") ? 16 : 13;
            const TString planeDirectory = Form("%s/%s", runDirectory.Data(), plane.Data());
            MakeDirectory(planeDirectory);

            for (const TString& type : types) {
                const TString typeDirectory = Form("%s/%s",
                                                   planeDirectory.Data(),
                                                   type.Data());
                MakeDirectory(typeDirectory);

                TString pdfPath = "";
                std::unique_ptr<TCanvas> canvas;
                if (savePlots) {
                    pdfPath = Form("%s/run%d_%s_%s_gaussian_fits.pdf",
                                   typeDirectory.Data(),
                                   run,
                                   plane.Data(),
                                   type.Data());
                    canvas.reset(new TCanvas(
                        Form("canvas_run%d_%s_%s", run, plane.Data(), type.Data()),
                        Form("run %d %s %s Gaussian fits", run, plane.Data(), type.Data()),
                        950,
                        700));
                    canvas->Print(Form("%s[", pdfPath.Data()));
                }

                for (int paddle = 0; paddle < numberOfPaddles; ++paddle) {
                    ++requested;

                    TString histogramName;
                    TString histogramPath;
                    std::unique_ptr<TH1> histogram = GetHistogramClone(
                        file.get(),
                        run,
                        plane,
                        paddle,
                        type,
                        histogramName,
                        histogramPath,
                        allowRecursiveSearch);

                    HistogramResult result = AnalyzeHistogram(histogram.get(),
                                                              run,
                                                              plane,
                                                              paddle,
                                                              type,
                                                              histogramName,
                                                              histogramPath);

                    if (result.found) ++found;
                    if (result.legacy.ok) ++legacyFitsOk;
                    if (result.type == "QDC" && result.twoPeakFitOk) ++qdcTwoPeakFitsOk;

                    WriteCsvRow(csv, result);

                    std::cout << "  run " << run
                              << " " << plane << Form("%02d", paddle)
                              << " " << type;
                    if (type == "QDC") {
                        if (result.twoPeakFitOk) {
                            std::cout << " : pedestal = "
                                      << FormatHundredth(result.pedestal.mean)
                                      << " +/- "
                                      << FormatHundredth(result.pedestal.meanError)
                                      << ", signal = "
                                      << FormatHundredth(result.signal.mean)
                                      << " +/- "
                                      << FormatHundredth(result.signal.meanError)
                                      << ", separation = "
                                      << FormatHundredth(result.separation)
                                      << " +/- "
                                      << FormatHundredth(result.separationError);
                        } else {
                            std::cout << " : TWO-PEAK CHECK: " << result.twoPeakMessage;
                        }
                    } else {
                        std::cout << " : peak = "
                                  << FormatHundredth(result.legacy.mean)
                                  << " +/- "
                                  << FormatHundredth(result.legacy.meanError)
                                  << " [" << result.legacy.message << "]";
                    }
                    std::cout << std::endl;

                    if (savePlots && canvas) {
                        const TString status = (type == "QDC")
                                             ? (result.twoPeakFitOk ? "OK" : "CHECK")
                                             : (result.legacy.ok ? "OK" : "CHECK");
                        const TString pngPath = Form("%s/run%d_%s%02d_%s_%s.png",
                                                     typeDirectory.Data(),
                                                     run,
                                                     plane.Data(),
                                                     paddle,
                                                     type.Data(),
                                                     status.Data());
                        DrawResult(canvas.get(),
                                   histogram.get(),
                                   result,
                                   pngPath,
                                   pdfPath);
                    }
                }

                if (savePlots && canvas) {
                    canvas->Print(Form("%s]", pdfPath.Data()));
                    canvas.reset();
                }
            }
        }

        file->Close();
        file.reset();
        gDirectory = nullptr;
    }

    csv.close();

    std::cout << "\nDone.\n"
              << "Requested histograms:       " << requested << "\n"
              << "Found histograms:           " << found << "\n"
              << "Successful legacy fits:     " << legacyFitsOk << "\n"
              << "Successful QDC 2-peak fits: " << qdcTwoPeakFitsOk << "\n"
              << "CSV summary:                " << csvName << "\n"
              << "Plots directory:            " << outputDirectory << "\n"
              << "savePlots:                  " << (savePlots ? "true" : "false") << "\n"
              << "recursive search:           " << (allowRecursiveSearch ? "true" : "false")
              << std::endl;
}

// Wrapper matching the filename, so ROOT can run the macro directly.
void BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded(const char* inputDir = ".",
                                                   const char* outputDir = "BH_peak_output",
                                                   int firstRun = 35566,
                                                   int lastRun = 35572,
                                                   bool savePlots = true,
                                                   bool allowRecursiveSearch = false)
{
    BHPeakFits_QDC_TDC_memsafe_fixedpaths(inputDir,
                                          outputDir,
                                          firstRun,
                                          lastRun,
                                          savePlots,
                                          allowRecursiveSearch);
}

// Backward-compatible wrapper retained from the earlier macro.
void BHPeakFits_QDC_TDC_memsafe(const char* inputDir = ".",
                                const char* outputDir = "BH_peak_output",
                                int firstRun = 35566,
                                int lastRun = 35572,
                                bool savePlots = true,
                                bool allowRecursiveSearch = false)
{
    BHPeakFits_QDC_TDC_memsafe_fixedpaths(inputDir,
                                          outputDir,
                                          firstRun,
                                          lastRun,
                                          savePlots,
                                          allowRecursiveSearch);
}
