// BH_QDC_HVScan.C
//
// End-to-end ROOT macro for the BHD/BHC high-voltage scan.
//
// For each requested QDC histogram this macro:
//   1. finds the two dominant, separated peaks;
//   2. fits the left peak (pedestal) and right peak (signal/HV peak)
//      with separate Gaussian functions;
//   3. calculates the peak separation
//          separation = mu_signal - mu_pedestal;
//   4. propagates the fit-mean uncertainties in quadrature
//          separation_uncertainty = sqrt(mu_signal_error^2
//                                        + mu_pedestal_error^2);
//   5. writes a CSV summary;
//   6. makes one BHD summary plot (paddles 4-8) and one BHC summary
//      plot (paddles 4-11), with one sub-pad per physical paddle.
//
// The x axis is the VOLTAGE SETTING RELATIVE TO THE DEFAULT, because
// only the offsets (-0.5, 0.0, +0.5, +1.0 V) were supplied.  Absolute
// high voltage would require the default/readback voltage for each paddle.
//
// The run layout is inferred from a single default-run argument:
//   defaultRun + 0: BHC and BHD at default voltage
//   defaultRun + 1: BHD -0.5 V, paddles 4-8
//   defaultRun + 2: BHD +0.5 V, paddles 4-8
//   defaultRun + 3: BHD +1.0 V, paddles 4-8
//   defaultRun + 4: BHC -0.5 V, paddles 4-11
//   defaultRun + 5: BHC +0.5 V, paddles 4-11
//   defaultRun + 6: BHC +1.0 V, paddles 4-11
//
// Run from a shell, using the run numbers in the current request:
//   root -l -b -q 'BH_QDC_HVScan.C(".", "BH_QDC_HV_scan", 35566, true, true)'
//
// Arguments:
//   inputDir          directory containing runNNNNN.root files
//   outputDir         output directory
//   defaultRun        run with both planes at default voltage
//   saveDiagnostics   save a multipage PDF showing every two-peak fit
//   fitLinear         optionally fit separation versus voltage offset to pol1

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace BH_QDC_HVScan_Detail {

const char* kAlgorithmVersion =
    "HVScan robust separated-region Gaussian+background v6 (2026-07-14)";

struct PeakCandidate {
    int bin = -1;
    double x = std::numeric_limits<double>::quiet_NaN();
    double height = 0.0;
    double baseline = 0.0;
    double prominence = 0.0;
    int leftHalfWidthBins = 0;
    int rightHalfWidthBins = 0;
};

struct GaussianPeakFit {
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
    double backgroundAtPeak = std::numeric_limits<double>::quiet_NaN();
    double backgroundSlope = std::numeric_limits<double>::quiet_NaN();

    TString message = "";
};

struct ScanSetting {
    int run = -1;
    TString plane = "";
    double voltageOffset = 0.0;
    int firstPaddle = -1;
    int lastPaddle = -1;
};

struct ScanPoint {
    int run = -1;
    TString plane = "";
    int paddle = -1;
    double voltageOffset = 0.0;
    TString histogram = "";
    TString histogramPath = "";

    bool found = false;
    bool peaksFound = false;
    bool ok = false;

    int pedestalSeedBin = -1;
    int signalSeedBin = -1;
    double pedestalSeedX = std::numeric_limits<double>::quiet_NaN();
    double signalSeedX = std::numeric_limits<double>::quiet_NaN();
    int signalSearchStartBin = -1;
    int signalSearchEndBin = -1;
    int valleyBin = -1;
    double valleyX = std::numeric_limits<double>::quiet_NaN();
    TString seedMethod = "";

    GaussianPeakFit pedestal;
    GaussianPeakFit signal;

    double separation = std::numeric_limits<double>::quiet_NaN();
    double separationError = std::numeric_limits<double>::quiet_NaN();

    TString message = "";
};

bool IsFinite(double x)
{
    return std::isfinite(x);
}

TString CsvSafe(TString text)
{
    text.ReplaceAll(",", ";");
    text.ReplaceAll("\n", " ");
    text.ReplaceAll("\r", " ");
    return text;
}

std::vector<TString> QdcHistogramNameCandidates(const TString& plane, int paddle)
{
    const TString paddleTag = Form("%s%02d", plane.Data(), paddle);
    std::vector<TString> names;

    // Confirmed layout from the previously inspected ROOT file:
    //   BHC/Paddle00/QDC/qdc_hit_BHC00left
    //   BHD/Paddle00/QDC/qdc_hit_BHD00down
    if (plane == "BHC") {
        names.push_back(Form("qdc_hit_%sleft", paddleTag.Data()));
        names.push_back(Form("qdc_hit_%sright", paddleTag.Data()));
    } else if (plane == "BHD") {
        names.push_back(Form("qdc_hit_%sdown", paddleTag.Data()));
        names.push_back(Form("qdc_hit_%sup", paddleTag.Data()));
    }

    return names;
}

TH1* GetQdcHistogram(TFile* file,
                     const TString& plane,
                     int paddle,
                     TString& usedName,
                     TString& usedPath)
{
    usedName = "";
    usedPath = "";

    if (file == nullptr || file->IsZombie()) return nullptr;

    const std::vector<TString> names = QdcHistogramNameCandidates(plane, paddle);

    for (const TString& name : names) {
        const std::vector<TString> paths = {
            Form("%s/Paddle%02d/QDC/%s", plane.Data(), paddle, name.Data()),
            Form("%s/%s", plane.Data(), name.Data())
        };

        for (const TString& path : paths) {
            TObject* object = file->Get(path.Data());
            if (object != nullptr && object->InheritsFrom(TH1::Class())) {
                usedName = name;
                usedPath = path;
                return dynamic_cast<TH1*>(object);
            }
        }
    }

    return nullptr;
}

// Peak finder used by THIS high-voltage-scan macro.
//
// The older implementation selected the two tallest local maxima after only a
// light TH1::Smooth call.  That allowed fluctuations on the falling pedestal
// shoulder (for example, channel 205 in run 35569 BHD04) to be selected as the
// signal.  It also used the first-bin overflow-like spike when setting the
// detection threshold, which hid the physical peaks in BHD06.
//
// The new method is deliberately different:
//   * build a Gaussian-smoothed copy in memory (the original histogram is never
//     modified and is still used for the fits);
//   * find the pedestal only in the early-QDC region;
//   * estimate the pedestal half-height width;
//   * begin the signal search at least four pedestal half-widths to the right;
//   * choose the strongest broad structure in that separated signal region.
//
// For the supplied 4096-channel histograms this finds approximately
//   run 35569 BHD04: pedestal 138, signal 846
//   run 35569 BHD06: pedestal 144, signal 752.
std::vector<double> BuildGaussianSmoothedCounts(TH1* histogram,
                                                double sigmaBins = 8.0)
{
    std::vector<double> smoothed;
    if (histogram == nullptr) return smoothed;

    const int numberOfBins = histogram->GetNbinsX();
    smoothed.assign(numberOfBins + 1, 0.0); // ROOT-style 1-based indexing.
    if (numberOfBins < 3 || !(sigmaBins > 0.0)) return smoothed;

    const int radius = std::max(3, static_cast<int>(std::ceil(4.0 * sigmaBins)));
    std::vector<double> kernel(2 * radius + 1, 0.0);
    double normalization = 0.0;

    for (int offset = -radius; offset <= radius; ++offset) {
        const double weight = std::exp(-0.5 * offset * offset
                                       / (sigmaBins * sigmaBins));
        kernel[offset + radius] = weight;
        normalization += weight;
    }
    if (!(normalization > 0.0)) return smoothed;
    for (double& weight : kernel) weight /= normalization;

    for (int bin = 1; bin <= numberOfBins; ++bin) {
        double value = 0.0;
        for (int offset = -radius; offset <= radius; ++offset) {
            const int sourceBin = std::max(1, std::min(numberOfBins, bin + offset));
            value += kernel[offset + radius] * histogram->GetBinContent(sourceBin);
        }
        smoothed[bin] = value;
    }

    return smoothed;
}

double PercentileInBinRange(const std::vector<double>& values,
                            int firstBin,
                            int lastBin,
                            double percentile)
{
    if (values.size() <= 1) return 0.0;
    const int maximumBin = static_cast<int>(values.size()) - 1;
    firstBin = std::max(1, std::min(maximumBin, firstBin));
    lastBin = std::max(1, std::min(maximumBin, lastBin));
    if (firstBin > lastBin) std::swap(firstBin, lastBin);

    std::vector<double> sample;
    sample.reserve(lastBin - firstBin + 1);
    for (int bin = firstBin; bin <= lastBin; ++bin) {
        if (IsFinite(values[bin])) sample.push_back(values[bin]);
    }
    if (sample.empty()) return 0.0;

    std::sort(sample.begin(), sample.end());
    percentile = std::max(0.0, std::min(1.0, percentile));
    const double index = percentile * (sample.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(index));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(index));
    const double fraction = index - lower;
    return sample[lower] * (1.0 - fraction) + sample[upper] * fraction;
}

int MaximumBinInRange(const std::vector<double>& values,
                      int firstBin,
                      int lastBin)
{
    if (values.size() <= 1) return -1;
    const int maximumBin = static_cast<int>(values.size()) - 1;
    firstBin = std::max(1, std::min(maximumBin, firstBin));
    lastBin = std::max(1, std::min(maximumBin, lastBin));
    if (firstBin > lastBin) std::swap(firstBin, lastBin);

    int bestBin = firstBin;
    double bestValue = values[firstBin];
    for (int bin = firstBin + 1; bin <= lastBin; ++bin) {
        if (values[bin] > bestValue) {
            bestValue = values[bin];
            bestBin = bin;
        }
    }
    return bestBin;
}

PeakCandidate CharacterizePeak(TH1* histogram,
                                const std::vector<double>& smoothed,
                                int peakBin,
                                int firstBin,
                                int lastBin)
{
    PeakCandidate candidate;
    if (histogram == nullptr || smoothed.size() <= 1) return candidate;

    const int numberOfBins = histogram->GetNbinsX();
    firstBin = std::max(1, std::min(numberOfBins, firstBin));
    lastBin = std::max(1, std::min(numberOfBins, lastBin));
    peakBin = std::max(firstBin, std::min(lastBin, peakBin));

    candidate.bin = peakBin;
    candidate.x = histogram->GetBinCenter(peakBin);
    candidate.height = smoothed[peakBin];
    candidate.baseline = std::max(0.0,
        PercentileInBinRange(smoothed, firstBin, lastBin, 0.10));
    candidate.prominence = std::max(0.0,
        candidate.height - candidate.baseline);

    const double halfProminenceLevel = candidate.baseline
                                     + 0.5 * candidate.prominence;

    int left = peakBin;
    while (left > firstBin && smoothed[left] > halfProminenceLevel) --left;

    int right = peakBin;
    while (right < lastBin && smoothed[right] > halfProminenceLevel) ++right;

    candidate.leftHalfWidthBins = std::max(1, peakBin - left);
    candidate.rightHalfWidthBins = std::max(1, right - peakBin);
    return candidate;
}

bool SelectPedestalAndSignal(TH1* histogram,
                             PeakCandidate& pedestal,
                             PeakCandidate& signal,
                             int& signalSearchStartBin,
                             int& signalSearchEndBin,
                             TString& seedMethod)
{
    if (histogram == nullptr || histogram->GetNbinsX() < 100) return false;

    const int numberOfBins = histogram->GetNbinsX();
    const std::vector<double> smoothed = BuildGaussianSmoothedCounts(histogram, 8.0);
    if (smoothed.size() != static_cast<std::size_t>(numberOfBins + 1)) return false;

    // Ignore edge spikes.  In the supplied files the first QDC bin can be much
    // taller than either physical peak and is not a useful pedestal seed.
    const int edgeGuard = std::max(16,
        static_cast<int>(std::lround(0.004 * numberOfBins)));

    // The pedestal is always in the early-QDC region.  Stopping at 8.5% of the
    // axis prevents the BHC signal peaks near channel 350 from becoming the
    // pedestal while still leaving a very generous window around channels 120-160.
    const int pedestalFirstBin = edgeGuard;
    const int pedestalLastBin = std::min(numberOfBins - edgeGuard,
        std::max(pedestalFirstBin + 40,
                 static_cast<int>(std::lround(0.085 * numberOfBins))));

    const int pedestalBin = MaximumBinInRange(smoothed,
                                               pedestalFirstBin,
                                               pedestalLastBin);
    if (pedestalBin < 1) return false;

    pedestal = CharacterizePeak(histogram,
                                smoothed,
                                pedestalBin,
                                pedestalFirstBin,
                                pedestalLastBin);

    const int pedestalHalfWidth = std::max(pedestal.leftHalfWidthBins,
                                           pedestal.rightHalfWidthBins);

    // The 7%-of-axis floor is about channel 287 for a 4096-bin QDC spectrum.
    // The four-half-width condition is what excludes the BHD04 shoulder at 205.
    signalSearchStartBin = std::max({
        pedestal.bin + 80,
        pedestal.bin + 4 * pedestalHalfWidth,
        static_cast<int>(std::lround(0.070 * numberOfBins))
    });
    signalSearchEndBin = std::min(numberOfBins - edgeGuard,
        static_cast<int>(std::lround(0.600 * numberOfBins)));

    if (signalSearchStartBin >= signalSearchEndBin) return false;

    const int signalBin = MaximumBinInRange(smoothed,
                                            signalSearchStartBin,
                                            signalSearchEndBin);
    if (signalBin < 1) return false;

    signal = CharacterizePeak(histogram,
                              smoothed,
                              signalBin,
                              signalSearchStartBin,
                              signalSearchEndBin);

    // Require a sustained structure rather than a one-bin electronics spike.
    const int signalWidth = signal.leftHalfWidthBins + signal.rightHalfWidthBins;
    const double minimumProminence = std::max(2.0, 0.005 * pedestal.height);
    if (!(signal.prominence > minimumProminence) || signalWidth < 8) return false;

    seedMethod = "Gaussian-smoothed separated-region maxima";
    return signal.bin > pedestal.bin;
}

int FindValleyBin(TH1* histogram, int leftPeakBin, int rightPeakBin)
{
    if (histogram == nullptr || leftPeakBin >= rightPeakBin) return -1;

    const std::vector<double> smoothed = BuildGaussianSmoothedCounts(histogram, 8.0);
    if (smoothed.empty()) return -1;

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

GaussianPeakFit FitGaussianPeak(TH1* histogram,
                                int run,
                                const TString& plane,
                                int paddle,
                                const TString& label,
                                const PeakCandidate& peak,
                                int lowerBoundaryBin,
                                int upperBoundaryBin)
{
    GaussianPeakFit result;
    if (histogram == nullptr || peak.bin < 1
        || peak.bin > histogram->GetNbinsX()) {
        result.message = "invalid histogram or peak seed";
        return result;
    }

    const int numberOfBins = histogram->GetNbinsX();
    lowerBoundaryBin = std::max(1, lowerBoundaryBin);
    upperBoundaryBin = std::min(numberOfBins, upperBoundaryBin);
    if (lowerBoundaryBin >= peak.bin || upperBoundaryBin <= peak.bin) {
        result.message = "peak seed outside fit boundaries";
        return result;
    }

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peak.bin);
    const double peakX = histogram->GetBinCenter(peak.bin);
    const double fullRange = histogram->GetXaxis()->GetXmax()
                           - histogram->GetXaxis()->GetXmin();

    // The geometric mean is a compromise between the short and long sides of a
    // skewed peak.  Using the shorter side alone made many of the previous fits
    // visibly too small; using the longer side alone lets the tail pull the mean.
    const double coreHalfWidthBins = std::max(3.0,
        std::sqrt(static_cast<double>(peak.leftHalfWidthBins)
                * static_cast<double>(peak.rightHalfWidthBins)));

    int halfRangeBins = std::max(8,
        static_cast<int>(std::lround(2.20 * coreHalfWidthBins)));

    // Keep the data interval exactly symmetric around the seed and inside the
    // valley/axis boundaries.
    halfRangeBins = std::min(halfRangeBins, peak.bin - lowerBoundaryBin - 1);
    halfRangeBins = std::min(halfRangeBins, upperBoundaryBin - peak.bin - 1);
    if (halfRangeBins < 5) {
        result.message = "not enough symmetric bins around peak";
        return result;
    }

    result.fitLow = peakX - halfRangeBins * binWidth;
    result.fitHigh = peakX + halfRangeBins * binWidth;

    const double sigmaGuess = std::max(binWidth,
        coreHalfWidthBins * binWidth / std::sqrt(2.0 * std::log(2.0)));
    const double backgroundGuess = std::max(0.0, peak.baseline);
    const double amplitudeGuess = std::max(1.0,
        peak.height - backgroundGuess);
    const double halfRangeX = halfRangeBins * binWidth;

    const TString functionName = Form("gaus_bg_%s_run%d_%s%02d",
                                      label.Data(), run, plane.Data(), paddle);

    // A symmetric Gaussian is fitted together with a local linear background.
    // The background absorbs a slowly varying pedestal/signal tail without
    // forcing the Gaussian mean and width to chase that asymmetric tail.
    TF1 model(functionName.Data(),
              "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[5])",
              result.fitLow,
              result.fitHigh);
    model.SetParNames("Amplitude", "Mean", "Sigma",
                      "Background", "BackgroundSlope", "SeedCenter");
    model.SetParameters(amplitudeGuess,
                        peakX,
                        sigmaGuess,
                        backgroundGuess,
                        0.0,
                        peakX);
    model.FixParameter(5, peakX);

    const double meanFreedom = std::max(2.0 * binWidth,
                                        0.45 * coreHalfWidthBins * binWidth);
    const double slopeLimit = std::max(1.0e-9,
        std::max(peak.height, 1.0) / std::max(halfRangeX, binWidth));

    model.SetParLimits(0, 0.0, std::max(10.0 * peak.height, 10.0));
    model.SetParLimits(1, peakX - meanFreedom, peakX + meanFreedom);
    model.SetParLimits(2,
                       std::max(0.35 * sigmaGuess, 0.5 * binWidth),
                       std::min(3.0 * sigmaGuess, fullRange));
    model.SetParLimits(3, 0.0, std::max(2.0 * peak.height, 10.0));
    model.SetParLimits(4, -slopeLimit, slopeLimit);

    // Fit the ORIGINAL unsmoothed histogram.  Smoothing is used only for seeds.
    // Q: quiet, R: range, 0: no automatic drawing, S: return result,
    // N: do not attach this temporary model to the histogram.
    TFitResultPtr fit = histogram->Fit(&model, "QR0SN");
    result.status = static_cast<int>(fit);

    result.amplitude = model.GetParameter(0);
    result.amplitudeError = model.GetParError(0);
    result.mean = model.GetParameter(1);
    result.meanError = model.GetParError(1);
    result.sigma = std::fabs(model.GetParameter(2));
    result.sigmaError = model.GetParError(2);
    result.backgroundAtPeak = model.GetParameter(3);
    result.backgroundSlope = model.GetParameter(4);
    result.chi2 = model.GetChisquare();
    result.ndf = model.GetNDF();

    result.ok = (result.status == 0
                 && IsFinite(result.mean)
                 && IsFinite(result.meanError)
                 && result.meanError > 0.0
                 && IsFinite(result.sigma)
                 && result.sigma > 0.0
                 && result.mean >= result.fitLow
                 && result.mean <= result.fitHigh);

    result.message = result.ok
        ? "ok (Gaussian plus local linear background)"
        : Form("fit status %d", result.status);
    return result;
}

ScanPoint AnalyzeHistogram(TH1* histogram,
                           int run,
                           const TString& plane,
                           int paddle,
                           double voltageOffset,
                           const TString& histogramName,
                           const TString& histogramPath)
{
    ScanPoint point;
    point.run = run;
    point.plane = plane;
    point.paddle = paddle;
    point.voltageOffset = voltageOffset;
    point.histogram = histogramName;
    point.histogramPath = histogramPath;
    point.found = (histogram != nullptr);

    if (histogram == nullptr) {
        point.message = "histogram not found";
        return point;
    }

    if (histogram->GetEntries() < 20) {
        point.message = "too few entries";
        return point;
    }

    PeakCandidate pedestalSeed;
    PeakCandidate signalSeed;
    point.peaksFound = SelectPedestalAndSignal(histogram,
                                               pedestalSeed,
                                               signalSeed,
                                               point.signalSearchStartBin,
                                               point.signalSearchEndBin,
                                               point.seedMethod);

    if (!point.peaksFound) {
        point.message = "could not identify separated pedestal and signal regions";
        return point;
    }

    point.pedestalSeedBin = pedestalSeed.bin;
    point.signalSeedBin = signalSeed.bin;
    point.pedestalSeedX = pedestalSeed.x;
    point.signalSeedX = signalSeed.x;

    point.valleyBin = FindValleyBin(histogram,
                                    pedestalSeed.bin,
                                    signalSeed.bin);
    if (point.valleyBin <= pedestalSeed.bin
        || point.valleyBin >= signalSeed.bin) {
        point.message = "could not identify valley between peaks";
        return point;
    }

    point.valleyX = histogram->GetBinCenter(point.valleyBin);

    point.pedestal = FitGaussianPeak(histogram,
                                     run,
                                     plane,
                                     paddle,
                                     "pedestal",
                                     pedestalSeed,
                                     1,
                                     point.valleyBin);

    point.signal = FitGaussianPeak(histogram,
                                   run,
                                   plane,
                                   paddle,
                                   "signal",
                                   signalSeed,
                                   point.valleyBin,
                                   histogram->GetNbinsX());

    point.ok = point.pedestal.ok
            && point.signal.ok
            && point.signal.mean > point.pedestal.mean;

    if (!point.ok) {
        point.message = Form("pedestal: %s; signal: %s",
                             point.pedestal.message.Data(),
                             point.signal.message.Data());
        return point;
    }

    point.separation = point.signal.mean - point.pedestal.mean;

    // The professor's requested quadrature rule.  These are the Gaussian
    // MEAN uncertainties, not the Gaussian widths sigma.
    point.separationError = std::hypot(point.signal.meanError,
                                       point.pedestal.meanError);

    point.message = Form("ok; seeds %.1f / %.1f; %s",
                         histogram->GetBinCenter(point.pedestalSeedBin),
                         histogram->GetBinCenter(point.signalSeedBin),
                         point.seedMethod.Data());
    return point;
}

void DrawDiagnostic(TCanvas* canvas,
                    TH1* histogram,
                    const ScanPoint& point,
                    const TString& diagnosticsPdf)
{
    if (canvas == nullptr || diagnosticsPdf.Length() == 0) return;

    canvas->Clear();
    canvas->cd();
    canvas->SetGrid();

    // Diagnostic-page layout only.  Keep the plot labels inside the page and
    // close to their axes without changing any fitting or analysis behavior.
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.13);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.04);

    TLatex text;
    text.SetNDC(kTRUE);
    text.SetTextSize(0.030);

    if (histogram == nullptr) {
        text.DrawLatex(0.12, 0.82,
                       Form("run %d %s%02d: histogram not found",
                            point.run, point.plane.Data(), point.paddle));
        canvas->Print(diagnosticsPdf.Data());
        return;
    }

    histogram->SetLineWidth(2);

    // Set the axis titles separately so the page title can be positioned
    // precisely with TLatex instead of ROOT's high default title box.
    histogram->SetTitle("");
    histogram->GetXaxis()->SetTitle(Form("%s channel", point.histogram.Data()));
    histogram->GetYaxis()->SetTitle("Counts");

    histogram->GetXaxis()->SetTitleSize(0.044);
    histogram->GetYaxis()->SetTitleSize(0.044);
    histogram->GetXaxis()->SetLabelSize(0.032);
    histogram->GetYaxis()->SetLabelSize(0.032);

    // Smaller offsets move each title toward its axis and away from the page
    // boundary.  This prevents clipping while leaving the data area intact.
    histogram->GetXaxis()->SetTitleOffset(1.02);
    histogram->GetYaxis()->SetTitleOffset(1.05);

    // Scale the diagnostic to the sustained physical spectrum rather than an
    // isolated first/last-bin electronics spike.  The histogram contents are
    // not modified; only the displayed y range is changed.
    const std::vector<double> diagnosticSmooth =
        BuildGaussianSmoothedCounts(histogram, 8.0);
    double physicalMaximum = 0.0;
    if (!diagnosticSmooth.empty()) {
        const int firstPhysicalBin = std::max(16,
            static_cast<int>(std::lround(0.004 * histogram->GetNbinsX())));
        const int lastPhysicalBin = std::min(histogram->GetNbinsX(),
            static_cast<int>(std::lround(0.600 * histogram->GetNbinsX())));
        for (int bin = firstPhysicalBin; bin <= lastPhysicalBin; ++bin) {
            physicalMaximum = std::max(physicalMaximum, diagnosticSmooth[bin]);
        }
    }
    if (physicalMaximum > 0.0) histogram->SetMaximum(1.25 * physicalMaximum);
    histogram->Draw("hist");

    TLatex plotTitle;
    plotTitle.SetNDC(kTRUE);
    plotTitle.SetTextAlign(22);
    plotTitle.SetTextSize(0.040);
    plotTitle.DrawLatex(0.50, 0.945,
                        Form("run %d, %s paddle %02d, #DeltaHV = %+.1f V",
                             point.run,
                             point.plane.Data(),
                             point.paddle,
                             point.voltageOffset));

    std::unique_ptr<TF1> pedestalFunction;
    std::unique_ptr<TF1> signalFunction;

    // FitGaussianPeak fits a Gaussian PLUS a local linear background.  The
    // amplitude stored in GaussianPeakFit is therefore the height ABOVE that
    // background, not the full histogram height.  Drawing a bare "gaus" with
    // that amplitude made the curves appear shifted downward by exactly the
    // fitted background.  Draw the same complete model used in the fit so the
    // diagnostic curve lies on the measured spectrum.
    const char* fittedModel =
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[5])";

    if (point.pedestal.ok) {
        const double pedestalSeedX =
            histogram->GetBinCenter(point.pedestalSeedBin);
        pedestalFunction.reset(new TF1(
            Form("draw_ped_run%d_%s%02d", point.run, point.plane.Data(), point.paddle),
            fittedModel,
            point.pedestal.fitLow,
            point.pedestal.fitHigh));
        pedestalFunction->SetParameters(point.pedestal.amplitude,
                                        point.pedestal.mean,
                                        point.pedestal.sigma,
                                        point.pedestal.backgroundAtPeak,
                                        point.pedestal.backgroundSlope,
                                        pedestalSeedX);
        pedestalFunction->SetLineColor(kBlue + 1);
        pedestalFunction->SetLineWidth(3);
        pedestalFunction->Draw("same");
    }

    if (point.signal.ok) {
        const double signalSeedX =
            histogram->GetBinCenter(point.signalSeedBin);
        signalFunction.reset(new TF1(
            Form("draw_sig_run%d_%s%02d", point.run, point.plane.Data(), point.paddle),
            fittedModel,
            point.signal.fitLow,
            point.signal.fitHigh));
        signalFunction->SetParameters(point.signal.amplitude,
                                      point.signal.mean,
                                      point.signal.sigma,
                                      point.signal.backgroundAtPeak,
                                      point.signal.backgroundSlope,
                                      signalSeedX);
        signalFunction->SetLineColor(kRed + 1);
        signalFunction->SetLineWidth(3);
        signalFunction->Draw("same");
    }

    std::unique_ptr<TLine> valleyLine;
    if (IsFinite(point.valleyX)) {
        valleyLine.reset(new TLine(point.valleyX,
                                   0.0,
                                   point.valleyX,
                                   histogram->GetMaximum()));
        valleyLine->SetLineStyle(2);
        valleyLine->SetLineColor(kGray + 2);
        valleyLine->Draw("same");
    }

    std::unique_ptr<TLine> pedestalSeedLine;
    std::unique_ptr<TLine> signalSeedLine;
    if (point.pedestalSeedBin > 0) {
        const double seedX = histogram->GetBinCenter(point.pedestalSeedBin);
        pedestalSeedLine.reset(new TLine(seedX, 0.0, seedX,
                                         histogram->GetMaximum()));
        pedestalSeedLine->SetLineStyle(3);
        pedestalSeedLine->SetLineColor(kBlue + 1);
        pedestalSeedLine->Draw("same");
    }
    if (point.signalSeedBin > 0) {
        const double seedX = histogram->GetBinCenter(point.signalSeedBin);
        signalSeedLine.reset(new TLine(seedX, 0.0, seedX,
                                       histogram->GetMaximum()));
        signalSeedLine->SetLineStyle(3);
        signalSeedLine->SetLineColor(kRed + 1);
        signalSeedLine->Draw("same");
    }

    text.SetTextColor(kBlack);
    text.SetTextSize(0.030);

    double y = 0.87;
    if (point.pedestal.ok) {
        text.DrawLatex(0.56, y,
                       Form("Pedestal #mu = %.2f #pm %.2f",
                            point.pedestal.mean,
                            point.pedestal.meanError));
        y -= 0.045;
    }
    if (point.signal.ok) {
        text.DrawLatex(0.56, y,
                       Form("Signal #mu = %.2f #pm %.2f",
                            point.signal.mean,
                            point.signal.meanError));
        y -= 0.045;
    }
    if (point.ok) {
        text.DrawLatex(0.56, y,
                       Form("Separation = %.2f #pm %.2f",
                            point.separation,
                            point.separationError));
        y -= 0.045;
    }
    if (!point.ok) {
        text.SetTextColor(kRed + 1);
        text.DrawLatex(0.12, 0.92, Form("CHECK: %s", point.message.Data()));
    }

    canvas->Print(diagnosticsPdf.Data());
}

void WriteCsvHeader(std::ofstream& output)
{
    output
        << "run,plane,paddle,voltage_offset_from_default_V,"
        << "histogram,histogram_path,found,two_peaks_found,fit_ok,"
        << "pedestal_seed_bin,pedestal_seed_channel,"
        << "signal_seed_bin,signal_seed_channel,"
        << "signal_search_start_bin,signal_search_end_bin,seed_method,"
        << "pedestal_peak,pedestal_peak_uncertainty,pedestal_sigma,"
        << "pedestal_sigma_uncertainty,pedestal_background,"
        << "pedestal_background_slope,pedestal_fit_status,"
        << "signal_peak,signal_peak_uncertainty,signal_sigma,"
        << "signal_sigma_uncertainty,signal_background,"
        << "signal_background_slope,signal_fit_status,"
        << "peak_separation,peak_separation_uncertainty,message\n";
}

void WriteCsvRow(std::ofstream& output, const ScanPoint& point)
{
    output << point.run << ","
           << point.plane << ","
           << point.paddle << ","
           << std::fixed << std::setprecision(2)
           << point.voltageOffset << ","
           << point.histogram << ","
           << point.histogramPath << ","
           << (point.found ? 1 : 0) << ","
           << (point.peaksFound ? 1 : 0) << ","
           << (point.ok ? 1 : 0) << ","
           << point.pedestalSeedBin << ","
           << point.pedestalSeedX << ","
           << point.signalSeedBin << ","
           << point.signalSeedX << ","
           << point.signalSearchStartBin << ","
           << point.signalSearchEndBin << ","
           << CsvSafe(point.seedMethod) << ","
           << point.pedestal.mean << ","
           << point.pedestal.meanError << ","
           << point.pedestal.sigma << ","
           << point.pedestal.sigmaError << ","
           << point.pedestal.backgroundAtPeak << ","
           << point.pedestal.backgroundSlope << ","
           << point.pedestal.status << ","
           << point.signal.mean << ","
           << point.signal.meanError << ","
           << point.signal.sigma << ","
           << point.signal.sigmaError << ","
           << point.signal.backgroundAtPeak << ","
           << point.signal.backgroundSlope << ","
           << point.signal.status << ","
           << point.separation << ","
           << point.separationError << ","
           << CsvSafe(point.message)
           << "\n";
}

void DrawPlaneSummary(const TString& plane,
                      int firstPaddle,
                      int lastPaddle,
                      const std::map<std::string,
                                     std::map<int, std::vector<ScanPoint>>>& allData,
                      const TString& outputDirectory,
                      bool fitLinear)
{
    auto planeIterator = allData.find(plane.Data());
    if (planeIterator == allData.end()) {
        std::cerr << "No data available for " << plane << std::endl;
        return;
    }

    const int numberOfPaddles = lastPaddle - firstPaddle + 1;
    const int columns = (plane == "BHD") ? 3 : 4;
    const int rows = static_cast<int>(std::ceil(
        numberOfPaddles / static_cast<double>(columns)));

    TCanvas canvas(Form("canvas_%s_hv_scan", plane.Data()),
                   Form("%s QDC peak separation versus high voltage", plane.Data()),
                   1800,
                   1050);
    canvas.Divide(columns, rows, 0.002, 0.002);

    std::vector<TGraphErrors*> graphs;
    std::vector<TF1*> lineFits;

    const TString slopeCsvName = Form("%s/%s_QDC_HV_linear_fit_parameters.csv",
                                      outputDirectory.Data(), plane.Data());
    std::ofstream slopeCsv;
    if (fitLinear) {
        slopeCsv.open(slopeCsvName.Data());
        if (slopeCsv.is_open()) {
            slopeCsv << "plane,paddle,n_points,fit_status,intercept_channels,"
                     << "intercept_uncertainty_channels,slope_channels_per_V,"
                     << "slope_uncertainty_channels_per_V,chi2,ndf\n";
        }
    }

    for (int paddle = firstPaddle; paddle <= lastPaddle; ++paddle) {
        const int padNumber = paddle - firstPaddle + 1;
        canvas.cd(padNumber);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.14);

        auto paddleIterator = planeIterator->second.find(paddle);
        if (paddleIterator == planeIterator->second.end()) {
            TLatex missing;
            missing.SetNDC(kTRUE);
            missing.SetTextSize(0.05);
            missing.DrawLatex(0.20, 0.55,
                              Form("No data for %s paddle %02d",
                                   plane.Data(), paddle));
            if (slopeCsv.is_open()) {
                slopeCsv << plane << "," << paddle
                         << ",0,-999,nan,nan,nan,nan,nan,nan\n";
            }
            continue;
        }

        std::vector<ScanPoint> points = paddleIterator->second;
        std::sort(points.begin(), points.end(),
                  [](const ScanPoint& a, const ScanPoint& b) {
                      return a.voltageOffset < b.voltageOffset;
                  });

        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> xError;
        std::vector<double> yError;
        std::vector<int> runs;

        for (const ScanPoint& point : points) {
            if (!point.ok) continue;
            x.push_back(point.voltageOffset);
            y.push_back(point.separation);
            xError.push_back(0.0); // No supply/readback uncertainty was provided.
            yError.push_back(point.separationError);
            runs.push_back(point.run);
        }

        if (x.empty()) {
            TLatex missing;
            missing.SetNDC(kTRUE);
            missing.SetTextSize(0.05);
            missing.DrawLatex(0.13, 0.55,
                              Form("No successful fits for %s paddle %02d",
                                   plane.Data(), paddle));
            if (slopeCsv.is_open()) {
                slopeCsv << plane << "," << paddle
                         << ",0,-998,nan,nan,nan,nan,nan,nan\n";
            }
            continue;
        }

        TGraphErrors* graph = new TGraphErrors(static_cast<int>(x.size()),
                                               x.data(),
                                               y.data(),
                                               xError.data(),
                                               yError.data());
        graphs.push_back(graph);

        graph->SetName(Form("graph_%s_paddle%02d", plane.Data(), paddle));
        graph->SetTitle(Form("%s paddle %02d;%s voltage setting relative to default (V);"
                             "QDC peak separation, #mu_{signal}-#mu_{pedestal} (channels)",
                             plane.Data(),
                             paddle,
                             plane.Data()));
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.15);
        graph->SetLineWidth(2);

        // Draw points and error bars only.  The trend is represented by the
        // uncertainty-weighted linear fit, not by dot-to-dot segments.
        graph->Draw("AP");
        graph->GetXaxis()->SetLimits(-0.70, 1.20);
        graph->GetXaxis()->SetNdivisions(505);
        graph->GetXaxis()->SetTitleSize(0.047);
        graph->GetYaxis()->SetTitleSize(0.047);
        graph->GetXaxis()->SetLabelSize(0.040);
        graph->GetYaxis()->SetLabelSize(0.040);
        graph->GetYaxis()->SetTitleOffset(1.35);

        int lineStatus = -997;
        double intercept = std::numeric_limits<double>::quiet_NaN();
        double interceptError = std::numeric_limits<double>::quiet_NaN();
        double slope = std::numeric_limits<double>::quiet_NaN();
        double slopeError = std::numeric_limits<double>::quiet_NaN();
        double chi2 = std::numeric_limits<double>::quiet_NaN();
        double ndf = std::numeric_limits<double>::quiet_NaN();

        if (fitLinear && x.size() >= 3) {
            TF1* line = new TF1(Form("line_%s_paddle%02d",
                                     plane.Data(), paddle),
                                "pol1", -0.5, 1.0);
            lineFits.push_back(line);
            line->SetLineStyle(2);
            line->SetLineWidth(3);

            TFitResultPtr lineResult = graph->Fit(line, "QRSN");
            lineStatus = static_cast<int>(lineResult);
            intercept = line->GetParameter(0);
            interceptError = line->GetParError(0);
            slope = line->GetParameter(1);
            slopeError = line->GetParError(1);
            chi2 = line->GetChisquare();
            ndf = line->GetNDF();

            if (lineStatus == 0) line->Draw("same");
        }

        if (slopeCsv.is_open()) {
            slopeCsv << plane << ","
                     << paddle << ","
                     << x.size() << ","
                     << lineStatus << ","
                     << std::fixed << std::setprecision(6)
                     << intercept << ","
                     << interceptError << ","
                     << slope << ","
                     << slopeError << ","
                     << chi2 << ","
                     << ndf << "\n";
        }

        gPad->Modified();
        gPad->Update();

        const double yMinimum = graph->GetYaxis()->GetXmin();
        const double yMaximum = graph->GetYaxis()->GetXmax();
        const double ySpan = std::max(yMaximum - yMinimum, 1.0);

        TLatex runLabel;
        runLabel.SetTextSize(0.025);
        runLabel.SetTextAlign(21);
        for (std::size_t i = 0; i < x.size(); ++i) {
            runLabel.DrawLatex(x[i],
                               y[i] + 0.035 * ySpan,
                               Form("%d", runs[i]));
        }

        TLatex note;
        note.SetNDC(kTRUE);
        note.SetTextSize(0.026);
        note.DrawLatex(0.16, 0.92, "Point labels are run numbers");

        if (fitLinear && lineStatus == 0) {
            note.SetTextSize(0.030);
            note.DrawLatex(0.16, 0.86,
                           Form("slope = %.2f #pm %.2f channels/V",
                                slope, slopeError));
            note.DrawLatex(0.16, 0.81,
                           Form("#chi^{2}/ndf = %.2f/%d",
                                chi2, static_cast<int>(ndf)));
        } else if (fitLinear) {
            note.SetTextColor(kRed + 1);
            note.DrawLatex(0.16, 0.86,
                           Form("linear fit unavailable (status %d)", lineStatus));
        }
    }

    if (slopeCsv.is_open()) slopeCsv.close();

    canvas.cd();
    const TString pngName = Form("%s/%s_QDC_peak_separation_vs_HV.png",
                                 outputDirectory.Data(), plane.Data());
    const TString pdfName = Form("%s/%s_QDC_peak_separation_vs_HV.pdf",
                                 outputDirectory.Data(), plane.Data());
    canvas.SaveAs(pngName.Data());
    canvas.SaveAs(pdfName.Data());
}

} // namespace BH_QDC_HVScan_Detail

void BH_QDC_HVScan(const char* inputDir = ".",
                   const char* outputDir = "BH_QDC_HV_scan",
                   int defaultRun = 35566,
                   bool saveDiagnostics = true,
                   bool fitLinear = true)
{
    using namespace BH_QDC_HVScan_Detail;

    gROOT->SetBatch(kTRUE);
    std::cout << "BH_QDC_HVScan algorithm: " << kAlgorithmVersion << "\n";
    std::cout << "Default run: " << defaultRun
              << "; linear fits: " << (fitLinear ? "enabled" : "disabled")
              << "\n" << std::endl;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetFitFormat("6.2f");
    gStyle->SetEndErrorSize(4.0);
    TH1::AddDirectory(kFALSE);

    const TString outputDirectory(outputDir);
    gSystem->mkdir(outputDirectory.Data(), kTRUE);

    // Sorted in increasing voltage offset so the summary graphs run left to right.
    const std::vector<ScanSetting> settings = {
        {defaultRun + 1, "BHD", -0.5, 4, 8},
        {defaultRun + 0, "BHD",  0.0, 4, 8},
        {defaultRun + 2, "BHD", +0.5, 4, 8},
        {defaultRun + 3, "BHD", +1.0, 4, 8},

        {defaultRun + 4, "BHC", -0.5, 4, 11},
        {defaultRun + 0, "BHC",  0.0, 4, 11},
        {defaultRun + 5, "BHC", +0.5, 4, 11},
        {defaultRun + 6, "BHC", +1.0, 4, 11}
    };

    const TString csvName = Form("%s/BH_QDC_HV_scan_summary_runs%d_%d.csv",
                                 outputDirectory.Data(),
                                 defaultRun,
                                 defaultRun + 6);
    std::ofstream csv(csvName.Data());
    if (!csv.is_open()) {
        std::cerr << "ERROR: could not create " << csvName << std::endl;
        return;
    }
    WriteCsvHeader(csv);

    TString diagnosticsPdf = "";
    std::unique_ptr<TCanvas> diagnosticCanvas;

    if (saveDiagnostics) {
        diagnosticsPdf = Form("%s/BH_QDC_two_peak_fit_diagnostics.pdf",
                              outputDirectory.Data());
        diagnosticCanvas.reset(new TCanvas("qdc_two_peak_diagnostic_canvas",
                                           "QDC two-peak fit diagnostics",
                                           1050,
                                           760));
        diagnosticCanvas->Print(Form("%s[", diagnosticsPdf.Data()));
    }

    std::map<std::string, std::map<int, std::vector<ScanPoint>>> allData;

    int requested = 0;
    int found = 0;
    int successful = 0;

    for (const ScanSetting& setting : settings) {
        const TString fileName = Form("%s/run%d.root", inputDir, setting.run);
        std::unique_ptr<TFile> file(TFile::Open(fileName.Data(), "READ"));

        if (!file || file->IsZombie()) {
            std::cerr << "WARNING: could not open " << fileName << std::endl;

            for (int paddle = setting.firstPaddle;
                 paddle <= setting.lastPaddle;
                 ++paddle) {
                ++requested;
                ScanPoint missing;
                missing.run = setting.run;
                missing.plane = setting.plane;
                missing.paddle = paddle;
                missing.voltageOffset = setting.voltageOffset;
                missing.message = "ROOT file not found";
                WriteCsvRow(csv, missing);
                allData[setting.plane.Data()][paddle].push_back(missing);

                if (saveDiagnostics) {
                    DrawDiagnostic(diagnosticCanvas.get(),
                                   nullptr,
                                   missing,
                                   diagnosticsPdf);
                }
            }
            continue;
        }

        std::cout << "\nProcessing " << fileName
                  << " for " << setting.plane
                  << " at relative voltage setting "
                  << std::showpos << std::fixed << std::setprecision(1)
                  << setting.voltageOffset << " V"
                  << std::noshowpos << std::endl;

        for (int paddle = setting.firstPaddle;
             paddle <= setting.lastPaddle;
             ++paddle) {
            ++requested;

            TString histogramName;
            TString histogramPath;
            TH1* histogram = GetQdcHistogram(file.get(),
                                             setting.plane,
                                             paddle,
                                             histogramName,
                                             histogramPath);

            ScanPoint point = AnalyzeHistogram(histogram,
                                               setting.run,
                                               setting.plane,
                                               paddle,
                                               setting.voltageOffset,
                                               histogramName,
                                               histogramPath);

            if (point.found) ++found;
            if (point.ok) ++successful;

            WriteCsvRow(csv, point);
            allData[setting.plane.Data()][paddle].push_back(point);

            std::cout << "  " << setting.plane
                      << Form("%02d", paddle)
                      << ": ";
            if (point.ok) {
                std::cout << "seeds = " << std::fixed << std::setprecision(1)
                          << point.pedestalSeedX << " / " << point.signalSeedX
                          << "; pedestal = " << std::setprecision(2)
                          << point.pedestal.mean << " +/- " << point.pedestal.meanError
                          << ", signal = " << point.signal.mean
                          << " +/- " << point.signal.meanError
                          << ", separation = " << point.separation
                          << " +/- " << point.separationError;
            } else {
                std::cout << "CHECK: " << point.message;
            }
            std::cout << std::endl;

            if (saveDiagnostics) {
                DrawDiagnostic(diagnosticCanvas.get(),
                               histogram,
                               point,
                               diagnosticsPdf);
            }
        }

        file->Close();
        file.reset();
        gDirectory = nullptr;
    }

    if (saveDiagnostics && diagnosticCanvas) {
        diagnosticCanvas->Print(Form("%s]", diagnosticsPdf.Data()));
        diagnosticCanvas.reset();
    }

    csv.close();

    DrawPlaneSummary("BHD", 4, 8, allData, outputDirectory, fitLinear);
    DrawPlaneSummary("BHC", 4, 11, allData, outputDirectory, fitLinear);

    std::cout << "\nDone.\n"
              << "Requested QDC histograms: " << requested << "\n"
              << "Found histograms:         " << found << "\n"
              << "Successful two-peak fits: " << successful << "\n"
              << "CSV summary:              " << csvName << "\n"
              << "BHD summary:              "
              << outputDirectory << "/BHD_QDC_peak_separation_vs_HV.png\n"
              << "BHC summary:              "
              << outputDirectory << "/BHC_QDC_peak_separation_vs_HV.png\n";

    if (fitLinear) {
        std::cout << "BHD slope CSV:            "
                  << outputDirectory << "/BHD_QDC_HV_linear_fit_parameters.csv\n"
                  << "BHC slope CSV:            "
                  << outputDirectory << "/BHC_QDC_HV_linear_fit_parameters.csv\n";
    }

    if (saveDiagnostics) {
        std::cout << "Fit diagnostics:          " << diagnosticsPdf << "\n";
    }

    std::cout << std::endl;
}
