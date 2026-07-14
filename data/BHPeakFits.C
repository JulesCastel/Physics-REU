// BHPeakFits.C
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
//   root -l -b -q 'BHPeakFits.C(".", "BH_peak_output", 35566, 35572)'
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
    double baseline = 0.0;
    double prominence = 0.0;
    int leftHalfWidthBins = 0;
    int rightHalfWidthBins = 0;
    int widthBins = 0;
    double area = 0.0;
    double score = 0.0;
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
    int searchStartBin = -1;
    int searchEndBin = -1;
    int pedestalSeedBin = -1;
    int signalSeedBin = -1;
    int valleyBin = -1;
    double valleyX = std::numeric_limits<double>::quiet_NaN();
    GaussianFit pedestal;
    GaussianFit signal;
    double separation = std::numeric_limits<double>::quiet_NaN();
    double separationError = std::numeric_limits<double>::quiet_NaN();

    // Records how the QDC peak seeds were chosen.  Exceptional run-specific
    // overrides therefore remain visible and auditable in the output CSV.
    TString seedMethod = "";
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

void GetPhysicalSearchBounds(TH1* histogram,
                             int& searchStartBin,
                             int& searchEndBin)
{
    searchStartBin = -1;
    searchEndBin = -1;
    if (histogram == nullptr) return;

    const int n = histogram->GetNbinsX();
    const int leftGuard = std::max(8, static_cast<int>(std::lround(0.005 * n)));
    const int rightGuard = std::max(8, static_cast<int>(std::lround(0.020 * n)));

    searchStartBin = std::max(2, 1 + leftGuard);
    searchEndBin = std::min(n - 1, n - rightGuard);

    if (searchEndBin <= searchStartBin + 12) {
        searchStartBin = 2;
        searchEndBin = n - 1;
    }
}

std::vector<double> SmoothBinContents(TH1* histogram, int radius = -1)
{
    std::vector<double> smoothed;
    if (histogram == nullptr) return smoothed;

    const int n = histogram->GetNbinsX();
    smoothed.assign(n + 1, 0.0);
    if (radius < 0) {
        radius = std::max(3, static_cast<int>(std::lround(n / 512.0)));
    }

    for (int bin = 1; bin <= n; ++bin) {
        const int low = std::max(1, bin - radius);
        const int high = std::min(n, bin + radius);
        double weightedSum = 0.0;
        double weightSum = 0.0;

        for (int other = low; other <= high; ++other) {
            const double weight = static_cast<double>(radius + 1 - std::abs(other - bin));
            weightedSum += weight * histogram->GetBinContent(other);
            weightSum += weight;
        }
        smoothed[bin] = (weightSum > 0.0) ? weightedSum / weightSum : 0.0;
    }

    return smoothed;
}

int FindSmoothedMaximumNear(TH1* histogram,
                            const std::vector<double>& smoothed,
                            double targetX,
                            double halfWindowX,
                            int searchStartBin,
                            int searchEndBin)
{
    if (histogram == nullptr || smoothed.empty()) return -1;

    int low = histogram->GetXaxis()->FindBin(targetX - halfWindowX);
    int high = histogram->GetXaxis()->FindBin(targetX + halfWindowX);
    low = std::max(searchStartBin, low);
    high = std::min(searchEndBin, high);
    if (high < low) return -1;

    int bestBin = low;
    double bestValue = smoothed[low];
    for (int bin = low + 1; bin <= high; ++bin) {
        if (smoothed[bin] > bestValue) {
            bestValue = smoothed[bin];
            bestBin = bin;
        }
    }
    return bestBin;
}

void FillSeedCandidate(TH1* histogram,
                       const std::vector<double>& smoothed,
                       int bin,
                       PeakCandidate& candidate)
{
    candidate = PeakCandidate();
    if (histogram == nullptr || smoothed.empty() || bin < 1
        || bin > histogram->GetNbinsX()) return;

    candidate.bin = bin;
    candidate.x = histogram->GetBinCenter(bin);
    candidate.height = smoothed[bin];
}

// The attached run35569 ROOT file shows two pathological QDC histograms:
// BHD04 contains a pedestal shoulder near channel 205, while BHD06 contains
// a very large first-bin spike.  Their physical signal modes are near channels
// 825 and 770.  The override searches locally around those measured regions
// rather than hard-coding an exact ROOT bin.
bool ApplyKnownPeakSeedOverride(TH1* histogram,
                                int run,
                                const TString& plane,
                                int paddle,
                                PeakCandidate& pedestal,
                                PeakCandidate& signal,
                                int& searchStartBin,
                                int& searchEndBin,
                                TString& seedMethod)
{
    const bool matchingRun = (run == 35569 || run == 33569);
    if (!matchingRun || plane != "BHD" || (paddle != 4 && paddle != 6)) {
        return false;
    }

    GetPhysicalSearchBounds(histogram, searchStartBin, searchEndBin);
    const std::vector<double> smoothed = SmoothBinContents(histogram);
    if (smoothed.empty()) return false;

    const double signalTarget = (paddle == 4) ? 825.0 : 770.0;
    const int pedestalBin = FindSmoothedMaximumNear(histogram,
                                                     smoothed,
                                                     140.0,
                                                     70.0,
                                                     searchStartBin,
                                                     searchEndBin);
    const int signalBin = FindSmoothedMaximumNear(histogram,
                                                   smoothed,
                                                   signalTarget,
                                                   180.0,
                                                   searchStartBin,
                                                   searchEndBin);
    if (pedestalBin < 1 || signalBin <= pedestalBin) return false;

    FillSeedCandidate(histogram, smoothed, pedestalBin, pedestal);
    FillSeedCandidate(histogram, smoothed, signalBin, signal);
    seedMethod = Form("known-run override: run %d BHD%02d", run, paddle);
    return true;
}

// Independent fallback used only if broad-peak candidate selection fails.
// It uses separated global windows, so a small local bump on the pedestal's
// descending edge cannot become the signal seed.
bool SelectSeparatedWindowMaxima(TH1* histogram,
                                  PeakCandidate& pedestal,
                                  PeakCandidate& signal,
                                  int& searchStartBin,
                                  int& searchEndBin)
{
    if (histogram == nullptr || histogram->GetNbinsX() < 32) return false;

    GetPhysicalSearchBounds(histogram, searchStartBin, searchEndBin);
    const std::vector<double> smoothed = SmoothBinContents(histogram);
    if (smoothed.empty()) return false;

    const int span = searchEndBin - searchStartBin;
    const int pedestalEnd = std::min(searchEndBin,
        searchStartBin + static_cast<int>(std::lround(0.18 * span)));

    int pedestalBin = searchStartBin;
    for (int bin = searchStartBin + 1; bin <= pedestalEnd; ++bin) {
        if (smoothed[bin] > smoothed[pedestalBin]) pedestalBin = bin;
    }

    const int enforcedGap = std::max(80,
        static_cast<int>(std::lround(0.05 * span)));
    const int signalStart = pedestalBin + enforcedGap;
    if (signalStart >= searchEndBin) return false;

    int signalBin = signalStart;
    for (int bin = signalStart + 1; bin <= searchEndBin; ++bin) {
        if (smoothed[bin] > smoothed[signalBin]) signalBin = bin;
    }

    if (smoothed[pedestalBin] <= 0.0
        || smoothed[signalBin] < std::max(2.0, 0.01 * smoothed[pedestalBin])) {
        return false;
    }

    FillSeedCandidate(histogram, smoothed, pedestalBin, pedestal);
    FillSeedCandidate(histogram, smoothed, signalBin, signal);
    return signal.bin > pedestal.bin;
}

double RangeMinimum(const std::vector<double>& values, int low, int high)
{
    if (values.empty()) return 0.0;
    low = std::max(1, low);
    high = std::min(static_cast<int>(values.size()) - 1, high);
    if (high < low) return 0.0;

    double minimum = values[low];
    for (int i = low + 1; i <= high; ++i) minimum = std::min(minimum, values[i]);
    return minimum;
}

std::vector<PeakCandidate> FindBroadPeakCandidates(TH1* histogram,
                                                    const std::vector<double>& smoothed,
                                                    int searchStartBin,
                                                    int searchEndBin)
{
    std::vector<PeakCandidate> candidates;
    if (histogram == nullptr || smoothed.empty()) return candidates;

    double referenceMaximum = 0.0;
    for (int bin = searchStartBin; bin <= searchEndBin; ++bin) {
        referenceMaximum = std::max(referenceMaximum, smoothed[bin]);
    }
    if (!(referenceMaximum > 0.0)) return candidates;

    const double minimumProminence = std::max(2.0, 0.005 * referenceMaximum);
    const int localRadius = std::max(
        50,
        static_cast<int>(std::lround(0.05 * (searchEndBin - searchStartBin)))
    );

    for (int bin = searchStartBin + 1; bin <= searchEndBin - 1; ++bin) {
        const double center = smoothed[bin];
        if (!(center >= smoothed[bin - 1] && center > smoothed[bin + 1])) continue;

        const double leftMinimum = RangeMinimum(
            smoothed,
            std::max(searchStartBin, bin - localRadius),
            bin
        );
        const double rightMinimum = RangeMinimum(
            smoothed,
            bin,
            std::min(searchEndBin, bin + localRadius)
        );

        const double baseline = std::max(leftMinimum, rightMinimum);
        const double prominence = center - baseline;
        if (prominence < minimumProminence) continue;

        const double halfProminenceLevel = baseline + 0.5 * prominence;
        int left = bin;
        int right = bin;
        while (left > searchStartBin && smoothed[left] > halfProminenceLevel) --left;
        while (right < searchEndBin && smoothed[right] > halfProminenceLevel) ++right;

        const int leftHalfWidth = bin - left;
        const int rightHalfWidth = right - bin;
        const int width = leftHalfWidth + rightHalfWidth;
        if (width < 8 || std::min(leftHalfWidth, rightHalfWidth) < 3) continue;

        const int areaRadius = std::max(8, 2 * std::max(leftHalfWidth, rightHalfWidth));
        const int areaLow = std::max(searchStartBin, bin - areaRadius);
        const int areaHigh = std::min(searchEndBin, bin + areaRadius);
        double excessArea = 0.0;
        for (int areaBin = areaLow; areaBin <= areaHigh; ++areaBin) {
            excessArea += std::max(0.0, smoothed[areaBin] - baseline);
        }

        PeakCandidate candidate;
        candidate.bin = bin;
        candidate.x = histogram->GetBinCenter(bin);
        candidate.height = center;
        candidate.baseline = baseline;
        candidate.prominence = prominence;
        candidate.leftHalfWidthBins = leftHalfWidth;
        candidate.rightHalfWidthBins = rightHalfWidth;
        candidate.widthBins = width;
        candidate.area = excessArea;
        candidate.score = excessArea * std::sqrt(std::max(prominence, 1.0));
        candidates.push_back(candidate);
    }

    return candidates;
}

bool SelectPedestalAndSignal(TH1* histogram,
                             PeakCandidate& pedestal,
                             PeakCandidate& signal,
                             int& searchStartBin,
                             int& searchEndBin)
{
    if (histogram == nullptr || histogram->GetNbinsX() < 32) return false;

    GetPhysicalSearchBounds(histogram, searchStartBin, searchEndBin);
    const std::vector<double> smoothed = SmoothBinContents(histogram);
    const std::vector<PeakCandidate> candidates = FindBroadPeakCandidates(
        histogram,
        smoothed,
        searchStartBin,
        searchEndBin
    );
    if (candidates.size() < 2) return false;

    const int pedestalLimit = searchStartBin
        + static_cast<int>(std::lround(0.30 * (searchEndBin - searchStartBin)));

    bool pedestalFound = false;
    double bestPedestalScore = -1.0;
    for (const PeakCandidate& candidate : candidates) {
        if (candidate.bin > pedestalLimit) continue;
        const double leftPreference = 1.0
            + 0.0005 * static_cast<double>(candidate.bin - searchStartBin);
        const double score = candidate.prominence
                           * std::sqrt(static_cast<double>(candidate.widthBins))
                           / leftPreference;
        if (score > bestPedestalScore) {
            bestPedestalScore = score;
            pedestal = candidate;
            pedestalFound = true;
        }
    }
    if (!pedestalFound) return false;

    const int pedestalCoreHalfWidth = std::max(
        3,
        std::min(pedestal.leftHalfWidthBins, pedestal.rightHalfWidthBins)
    );
    const int signalSearchStart = pedestal.bin
        + std::max(12, 4 * pedestalCoreHalfWidth);

    bool signalFound = false;
    double bestSignalScore = -1.0;
    for (const PeakCandidate& candidate : candidates) {
        if (candidate.bin < signalSearchStart) continue;
        if (candidate.score > bestSignalScore) {
            bestSignalScore = candidate.score;
            signal = candidate;
            signalFound = true;
        }
    }

    return signalFound && signal.bin > pedestal.bin;
}

int FindValleyBin(TH1* histogram,
                  int leftPeakBin,
                  int rightPeakBin,
                  const std::vector<double>& smoothed)
{
    if (histogram == nullptr || smoothed.empty() || leftPeakBin >= rightPeakBin) return -1;

    int valleyBin = leftPeakBin + 1;
    double minimum = smoothed[valleyBin];
    for (int bin = leftPeakBin + 1; bin < rightPeakBin; ++bin) {
        if (smoothed[bin] < minimum) {
            minimum = smoothed[bin];
            valleyBin = bin;
        }
    }
    return valleyBin;
}

bool EstimateSymmetricCoreWindow(TH1* histogram,
                                 int peakBin,
                                 int lowerBoundaryBin,
                                 int upperBoundaryBin,
                                 double& fitLow,
                                 double& fitHigh,
                                 double& sigmaGuess,
                                 double& allowedMeanShift)
{
    if (histogram == nullptr) return false;

    const int n = histogram->GetNbinsX();
    lowerBoundaryBin = std::max(1, lowerBoundaryBin);
    upperBoundaryBin = std::min(n, upperBoundaryBin);
    if (peakBin <= lowerBoundaryBin || peakBin >= upperBoundaryBin) return false;

    const std::vector<double> smoothed = SmoothBinContents(histogram);
    if (smoothed.empty()) return false;

    const double peakHeight = smoothed[peakBin];
    if (!(peakHeight > 0.0)) return false;

    // Use 65% of the peak prominence above a local baseline.  This is a
    // narrower and more symmetric definition of the Gaussian core than the
    // old absolute half-height rule, so long QDC tails have less leverage.
    const double leftBase = RangeMinimum(smoothed, lowerBoundaryBin, peakBin);
    const double rightBase = RangeMinimum(smoothed, peakBin, upperBoundaryBin);
    const double localBaseline = std::max(leftBase, rightBase);
    const double coreLevel = localBaseline + 0.65 * (peakHeight - localBaseline);
    if (!(coreLevel < peakHeight)) return false;

    int halfLeft = peakBin;
    int halfRight = peakBin;
    while (halfLeft > lowerBoundaryBin && smoothed[halfLeft] > coreLevel) --halfLeft;
    while (halfRight < upperBoundaryBin && smoothed[halfRight] > coreLevel) ++halfRight;

    int coreHalfWidthBins = std::min(peakBin - halfLeft, halfRight - peakBin);
    if (coreHalfWidthBins < 3) coreHalfWidthBins = 3;

    int halfRangeBins = std::max(
        6,
        static_cast<int>(std::lround(1.10 * coreHalfWidthBins))
    );
    const int availableHalfRange = std::min(
        peakBin - lowerBoundaryBin,
        upperBoundaryBin - peakBin
    );
    halfRangeBins = std::min(halfRangeBins, availableHalfRange);
    if (halfRangeBins < 5) return false;

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    const double peakX = histogram->GetBinCenter(peakBin);
    fitLow = peakX - halfRangeBins * binWidth;
    fitHigh = peakX + halfRangeBins * binWidth;
    sigmaGuess = std::max(
        (coreHalfWidthBins * binWidth) / std::sqrt(2.0 * std::log(2.0)),
        binWidth
    );
    allowedMeanShift = std::max(2.0 * binWidth, 0.25 * coreHalfWidthBins * binWidth);
    return fitHigh > fitLow;
}

bool IsUsableGaussian(const GaussianFit& result)
{
    return result.status == 0
        && IsFinite(result.amplitude)
        && result.amplitude > 0.0
        && IsFinite(result.mean)
        && IsFinite(result.meanError)
        && result.meanError > 0.0
        && IsFinite(result.sigma)
        && result.sigma > 0.0
        && IsFinite(result.sigmaError)
        && result.sigmaError >= 0.0
        && result.mean >= result.fitLow
        && result.mean <= result.fitHigh;
}

void ReadGaussianResult(TF1& gaussian,
                        int status,
                        double fitLow,
                        double fitHigh,
                        GaussianFit& result)
{
    result.status = status;
    result.fitLow = fitLow;
    result.fitHigh = fitHigh;
    result.amplitude = gaussian.GetParameter(0);
    result.amplitudeError = gaussian.GetParError(0);
    result.mean = gaussian.GetParameter(1);
    result.meanError = gaussian.GetParError(1);
    result.sigma = std::fabs(gaussian.GetParameter(2));
    result.sigmaError = gaussian.GetParError(2);
    result.chi2 = gaussian.GetChisquare();
    result.ndf = gaussian.GetNDF();
    result.ok = IsUsableGaussian(result);
}

GaussianFit FitSymmetricGaussianCore(TH1* histogram,
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
    double allowedMeanShift = 0.0;
    if (!EstimateSymmetricCoreWindow(histogram,
                                     peakBin,
                                     lowerBoundaryBin,
                                     upperBoundaryBin,
                                     result.fitLow,
                                     result.fitHigh,
                                     sigmaGuess,
                                     allowedMeanShift)) {
        result.message = "could not construct symmetric core window";
        return result;
    }

    const double peakX = histogram->GetBinCenter(peakBin);
    double peakHeight = histogram->GetBinContent(peakBin);
    for (int bin = std::max(1, peakBin - 2);
         bin <= std::min(histogram->GetNbinsX(), peakBin + 2);
         ++bin) {
        peakHeight = std::max(peakHeight, histogram->GetBinContent(bin));
    }

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    const double minimumSigma = std::max(0.35 * sigmaGuess, 0.50 * binWidth);
    const double maximumSigma = std::max(2.50 * sigmaGuess, minimumSigma + binWidth);
    const double meanLow = std::max(result.fitLow, peakX - allowedMeanShift);
    const double meanHigh = std::min(result.fitHigh, peakX + allowedMeanShift);

    const TString firstName = Form("gaus_core1_%s_run%d_%s%02d",
                                   label.Data(), run, plane.Data(), paddle);
    TF1 first(firstName.Data(), "gaus", result.fitLow, result.fitHigh);
    first.SetParNames("Amplitude", "Mean", "Sigma");
    first.SetParameters(std::max(peakHeight, 1.0), peakX, sigmaGuess);
    first.SetParLimits(0, 0.0, std::max(10.0 * peakHeight, 1.0));
    first.SetParLimits(1, meanLow, meanHigh);
    first.SetParLimits(2, minimumSigma, maximumSigma);

    TFitResultPtr firstFit = histogram->Fit(&first, "QRL0SN");
    ReadGaussianResult(first,
                       static_cast<int>(firstFit),
                       result.fitLow,
                       result.fitHigh,
                       result);

    if (!result.ok) {
        result.message = Form("initial symmetric-core fit status %d", result.status);
        return result;
    }

    const double originalHalfRange = 0.5 * (result.fitHigh - result.fitLow);
    const double centeredAvailableHalfRange = std::min(
        result.mean - result.fitLow,
        result.fitHigh - result.mean
    );
    if (centeredAvailableHalfRange < 5.0 * binWidth) {
        result.message = "ok (initial symmetric core; insufficient room to refine)";
        return result;
    }

    const double refinedHalfRange = std::min(
        centeredAvailableHalfRange,
        std::max(5.0 * binWidth, 1.40 * result.sigma)
    );
    const double refinedLow = result.mean - refinedHalfRange;
    const double refinedHigh = result.mean + refinedHalfRange;

    const TString secondName = Form("gaus_core2_%s_run%d_%s%02d",
                                    label.Data(), run, plane.Data(), paddle);
    TF1 second(secondName.Data(), "gaus", refinedLow, refinedHigh);
    second.SetParNames("Amplitude", "Mean", "Sigma");
    second.SetParameters(result.amplitude, result.mean, result.sigma);
    second.SetParLimits(0, 0.0, std::max(10.0 * peakHeight, 1.0));
    second.SetParLimits(1, meanLow, meanHigh);
    second.SetParLimits(2, minimumSigma, maximumSigma);

    TFitResultPtr secondFit = histogram->Fit(&second, "QRL0SN");
    GaussianFit refined;
    ReadGaussianResult(second,
                       static_cast<int>(secondFit),
                       refinedLow,
                       refinedHigh,
                       refined);

    if (refined.ok) {
        refined.message = "ok (refined symmetric core)";
        return refined;
    }

    result.message = "ok (initial symmetric core; refinement failed)";
    return result;
}

void SetQdcDiagnosticDisplayRange(TH1* histogram,
                                  int searchStartBin,
                                  int searchEndBin)
{
    if (histogram == nullptr) return;
    if (searchStartBin < 1 || searchEndBin <= searchStartBin) {
        GetPhysicalSearchBounds(histogram, searchStartBin, searchEndBin);
    }

    const std::vector<double> smoothed = SmoothBinContents(histogram);
    double physicalMaximum = 0.0;
    for (int bin = searchStartBin; bin <= searchEndBin; ++bin) {
        physicalMaximum = std::max(physicalMaximum, smoothed[bin]);
    }
    if (physicalMaximum > 0.0) {
        histogram->SetMaximum(1.35 * physicalMaximum);
        histogram->SetMinimum(0.0);
    }
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

    PeakCandidate pedestalSeed;
    PeakCandidate signalSeed;

    result.twoPeaksFound = ApplyKnownPeakSeedOverride(histogram,
                                                      run,
                                                      plane,
                                                      paddle,
                                                      pedestalSeed,
                                                      signalSeed,
                                                      result.searchStartBin,
                                                      result.searchEndBin,
                                                      result.seedMethod);
    if (!result.twoPeaksFound) {
        result.twoPeaksFound = SelectPedestalAndSignal(histogram,
                                                       pedestalSeed,
                                                       signalSeed,
                                                       result.searchStartBin,
                                                       result.searchEndBin);
        if (result.twoPeaksFound) {
            result.seedMethod = "automatic broad-peak search";
        }
    }
    if (!result.twoPeaksFound) {
        result.twoPeaksFound = SelectSeparatedWindowMaxima(histogram,
                                                           pedestalSeed,
                                                           signalSeed,
                                                           result.searchStartBin,
                                                           result.searchEndBin);
        if (result.twoPeaksFound) {
            result.seedMethod = "separated-window fallback";
        }
    }
    if (!result.twoPeaksFound) {
        result.twoPeakMessage = "could not identify pedestal and signal with any seed method";
        return result;
    }

    result.pedestalSeedBin = pedestalSeed.bin;
    result.signalSeedBin = signalSeed.bin;

    const std::vector<double> smoothed = SmoothBinContents(histogram);
    result.valleyBin = FindValleyBin(histogram,
                                     pedestalSeed.bin,
                                     signalSeed.bin,
                                     smoothed);
    if (result.valleyBin <= pedestalSeed.bin || result.valleyBin >= signalSeed.bin) {
        result.twoPeakMessage = "could not identify valley between peaks";
        return result;
    }
    result.valleyX = histogram->GetBinCenter(result.valleyBin);

    result.pedestal = FitSymmetricGaussianCore(histogram,
                                                run,
                                                plane,
                                                paddle,
                                                "pedestal",
                                                pedestalSeed.bin,
                                                result.searchStartBin,
                                                result.valleyBin);
    result.signal = FitSymmetricGaussianCore(histogram,
                                              run,
                                              plane,
                                              paddle,
                                              "signal",
                                              signalSeed.bin,
                                              result.valleyBin,
                                              result.searchEndBin);

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

    if (result.type == "QDC") {
        SetQdcDiagnosticDisplayRange(histogram,
                                     result.searchStartBin,
                                     result.searchEndBin);
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
        << "two_peaks_found,two_peak_fit_ok,seed_method,"
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
        << (result.twoPeakFitOk ? 1 : 0) << ","
        << CsvSafe(result.seedMethod) << ",";

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
                            std::cout << " : [" << result.seedMethod << "] pedestal = "
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


// Wrapper matching BHPeakFits.C, so ROOT can run the renamed macro directly.
void BHPeakFits(const char* inputDir = ".",
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

// Backward-compatible wrapper matching the earlier revised filename.
void BHPeakFits_QDC_TDC_memsafe_fixedpaths_rounded_revised(const char* inputDir = ".",
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
