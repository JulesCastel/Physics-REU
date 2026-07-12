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
//   root -l -b -q 'BH_QDC_HVScan.C(".", "BH_QDC_HV_scan", 33566, true, false)'
//
// The earlier Gaussian-fit work used run numbers beginning with 35566.
// If those are the actual file names, use:
//   root -l -b -q 'BH_QDC_HVScan.C(".", "BH_QDC_HV_scan", 35566, true, false)'
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

struct PeakCandidate {
    int bin = -1;
    double x = std::numeric_limits<double>::quiet_NaN();
    double height = 0.0;
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

    int valleyBin = -1;
    double valleyX = std::numeric_limits<double>::quiet_NaN();

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

std::vector<PeakCandidate> FindLocalMaxima(TH1* histogram, double relativeThreshold)
{
    std::vector<PeakCandidate> candidates;
    if (histogram == nullptr || histogram->GetNbinsX() < 7) return candidates;

    std::unique_ptr<TH1> smoothed(dynamic_cast<TH1*>(histogram->Clone(
        Form("%s_smoothed_for_peak_search", histogram->GetName()))));
    if (!smoothed) return candidates;

    smoothed->SetDirectory(nullptr);
    smoothed->Smooth(3);

    const int numberOfBins = smoothed->GetNbinsX();
    const double maximum = smoothed->GetMaximum();
    if (!(maximum > 0.0)) return candidates;

    const double threshold = std::max(0.0, relativeThreshold) * maximum;

    for (int bin = 3; bin <= numberOfBins - 2; ++bin) {
        const double center = smoothed->GetBinContent(bin);
        const double left1 = smoothed->GetBinContent(bin - 1);
        const double right1 = smoothed->GetBinContent(bin + 1);

        if (center < threshold) continue;
        if (!(center >= left1 && center > right1)) continue;

        PeakCandidate candidate;
        candidate.bin = bin;
        candidate.x = smoothed->GetBinCenter(bin);
        candidate.height = center;
        candidates.push_back(candidate);
    }

    return candidates;
}

bool SelectTwoSeparatedPeaks(TH1* histogram,
                             PeakCandidate& leftPeak,
                             PeakCandidate& rightPeak,
                             double relativeThreshold = 0.03)
{
    if (histogram == nullptr) return false;

    std::vector<PeakCandidate> candidates = FindLocalMaxima(histogram, relativeThreshold);
    if (candidates.size() < 2 && relativeThreshold > 0.01) {
        candidates = FindLocalMaxima(histogram, 0.01);
    }
    if (candidates.size() < 2) return false;

    std::sort(candidates.begin(), candidates.end(),
              [](const PeakCandidate& a, const PeakCandidate& b) {
                  return a.height > b.height;
              });

    const double xRange = histogram->GetXaxis()->GetXmax()
                        - histogram->GetXaxis()->GetXmin();
    const double binWidth = histogram->GetXaxis()->GetBinWidth(1);
    const double minimumSeparation = std::max(6.0 * binWidth, 0.03 * xRange);

    bool pairFound = false;
    double bestScore = -1.0;
    PeakCandidate bestA;
    PeakCandidate bestB;

    const std::size_t maximumCandidates = std::min<std::size_t>(candidates.size(), 12);

    for (std::size_t i = 0; i < maximumCandidates; ++i) {
        for (std::size_t j = i + 1; j < maximumCandidates; ++j) {
            const double separation = std::fabs(candidates[i].x - candidates[j].x);
            if (separation < minimumSeparation) continue;

            // Favor a pair in which BOTH peaks are prominent.  A small separation
            // term breaks near-ties in favor of clearly distinct peaks.
            const double smallerHeight = std::min(candidates[i].height,
                                                  candidates[j].height);
            const double score = smallerHeight
                               + 0.05 * (candidates[i].height + candidates[j].height)
                               + 0.001 * separation;

            if (score > bestScore) {
                bestScore = score;
                bestA = candidates[i];
                bestB = candidates[j];
                pairFound = true;
            }
        }
    }

    if (!pairFound) return false;

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
    if (histogram == nullptr) return -1;
    if (leftPeakBin >= rightPeakBin) return -1;

    std::unique_ptr<TH1> smoothed(dynamic_cast<TH1*>(histogram->Clone(
        Form("%s_smoothed_for_valley", histogram->GetName()))));
    if (!smoothed) return -1;

    smoothed->SetDirectory(nullptr);
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

void EstimateFitWindow(TH1* histogram,
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
           && histogram->GetBinContent(halfLeft) > halfHeight) {
        --halfLeft;
    }
    while (halfRight < upperBoundaryBin
           && histogram->GetBinContent(halfRight) > halfHeight) {
        ++halfRight;
    }

    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);
    double fwhm = histogram->GetBinLowEdge(halfRight + 1)
                - histogram->GetBinLowEdge(halfLeft);

    if (!IsFinite(fwhm) || fwhm < 2.0 * binWidth) {
        fwhm = 6.0 * binWidth;
    }

    sigmaGuess = std::max(fwhm / 2.355, binWidth);

    const double peakX = histogram->GetBinCenter(peakBin);
    const double lowerBoundaryX = histogram->GetBinLowEdge(lowerBoundaryBin);
    const double upperBoundaryX = histogram->GetBinLowEdge(upperBoundaryBin + 1);

    fitLow = std::max(lowerBoundaryX, peakX - 2.8 * sigmaGuess);
    fitHigh = std::min(upperBoundaryX, peakX + 2.8 * sigmaGuess);

    // Guarantee enough width for a stable three-parameter Gaussian fit.
    if (fitHigh - fitLow < 5.0 * binWidth) {
        fitLow = std::max(lowerBoundaryX, peakX - 3.0 * binWidth);
        fitHigh = std::min(upperBoundaryX, peakX + 3.0 * binWidth);
    }
}

GaussianPeakFit FitGaussianPeak(TH1* histogram,
                                int run,
                                const TString& plane,
                                int paddle,
                                const TString& label,
                                int peakBin,
                                int lowerBoundaryBin,
                                int upperBoundaryBin)
{
    GaussianPeakFit result;

    if (histogram == nullptr || peakBin < 1 || peakBin > histogram->GetNbinsX()) {
        result.message = "invalid histogram or peak bin";
        return result;
    }

    double sigmaGuess = 0.0;
    EstimateFitWindow(histogram,
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
    const double fullRange = histogram->GetXaxis()->GetXmax()
                           - histogram->GetXaxis()->GetXmin();
    const double binWidth = histogram->GetXaxis()->GetBinWidth(peakBin);

    const TString functionName = Form("gaus_%s_run%d_%s%02d",
                                      label.Data(), run, plane.Data(), paddle);
    TF1 gaussian(functionName.Data(), "gaus", result.fitLow, result.fitHigh);
    gaussian.SetParNames("Amplitude", "Mean", "Sigma");
    gaussian.SetParameters(peakHeight, peakX, sigmaGuess);
    gaussian.SetParLimits(0, 0.0, std::max(10.0 * peakHeight, 1.0));
    gaussian.SetParLimits(1, result.fitLow, result.fitHigh);
    gaussian.SetParLimits(2, std::max(0.25 * binWidth, 1.0e-9), fullRange);

    // Q: quiet, R: use range, 0: do not draw, S: return fit result,
    // N: do not attach the fit function to the histogram.
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

    result.ok = (result.status == 0
                 && IsFinite(result.mean)
                 && IsFinite(result.meanError)
                 && result.meanError > 0.0
                 && IsFinite(result.sigma)
                 && result.sigma > 0.0
                 && result.mean >= result.fitLow
                 && result.mean <= result.fitHigh);

    result.message = result.ok ? "ok" : Form("fit status %d", result.status);
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

    PeakCandidate leftPeak;
    PeakCandidate rightPeak;
    point.peaksFound = SelectTwoSeparatedPeaks(histogram, leftPeak, rightPeak);

    if (!point.peaksFound) {
        point.message = "could not identify two separated peaks";
        return point;
    }

    point.valleyBin = FindValleyBin(histogram, leftPeak.bin, rightPeak.bin);
    if (point.valleyBin <= leftPeak.bin || point.valleyBin >= rightPeak.bin) {
        point.message = "could not identify valley between peaks";
        return point;
    }

    point.valleyX = histogram->GetBinCenter(point.valleyBin);

    point.pedestal = FitGaussianPeak(histogram,
                                     run,
                                     plane,
                                     paddle,
                                     "pedestal",
                                     leftPeak.bin,
                                     1,
                                     point.valleyBin);

    point.signal = FitGaussianPeak(histogram,
                                   run,
                                   plane,
                                   paddle,
                                   "signal",
                                   rightPeak.bin,
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

    point.message = "ok";
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
    histogram->SetTitle(Form("run %d, %s paddle %02d, #DeltaHV = %+.1f V;%s channel;Counts",
                             point.run,
                             point.plane.Data(),
                             point.paddle,
                             point.voltageOffset,
                             point.histogram.Data()));
    histogram->Draw("hist");

    std::unique_ptr<TF1> pedestalFunction;
    std::unique_ptr<TF1> signalFunction;

    if (point.pedestal.ok) {
        pedestalFunction.reset(new TF1(
            Form("draw_ped_run%d_%s%02d", point.run, point.plane.Data(), point.paddle),
            "gaus",
            point.pedestal.fitLow,
            point.pedestal.fitHigh));
        pedestalFunction->SetParameters(point.pedestal.amplitude,
                                        point.pedestal.mean,
                                        point.pedestal.sigma);
        pedestalFunction->SetLineColor(kBlue + 1);
        pedestalFunction->SetLineWidth(3);
        pedestalFunction->Draw("same");
    }

    if (point.signal.ok) {
        signalFunction.reset(new TF1(
            Form("draw_sig_run%d_%s%02d", point.run, point.plane.Data(), point.paddle),
            "gaus",
            point.signal.fitLow,
            point.signal.fitHigh));
        signalFunction->SetParameters(point.signal.amplitude,
                                      point.signal.mean,
                                      point.signal.sigma);
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
        << "pedestal_peak,pedestal_peak_uncertainty,pedestal_sigma,"
        << "pedestal_sigma_uncertainty,pedestal_fit_status,"
        << "signal_peak,signal_peak_uncertainty,signal_sigma,"
        << "signal_sigma_uncertainty,signal_fit_status,"
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
           << point.pedestal.mean << ","
           << point.pedestal.meanError << ","
           << point.pedestal.sigma << ","
           << point.pedestal.sigmaError << ","
           << point.pedestal.status << ","
           << point.signal.mean << ","
           << point.signal.meanError << ","
           << point.signal.sigma << ","
           << point.signal.sigmaError << ","
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
    const int rows = static_cast<int>(std::ceil(numberOfPaddles / static_cast<double>(columns)));

    TCanvas canvas(Form("canvas_%s_hv_scan", plane.Data()),
                   Form("%s QDC peak separation versus high voltage", plane.Data()),
                   1800,
                   1050);
    canvas.Divide(columns, rows, 0.002, 0.002);

    std::vector<TGraphErrors*> graphs;
    std::vector<TF1*> lineFits;

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
                              Form("No data for %s paddle %02d", plane.Data(), paddle));
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
        graph->Draw("APL");
        graph->GetXaxis()->SetLimits(-0.70, 1.20);
        graph->GetXaxis()->SetNdivisions(505);
        graph->GetXaxis()->SetTitleSize(0.047);
        graph->GetYaxis()->SetTitleSize(0.047);
        graph->GetXaxis()->SetLabelSize(0.040);
        graph->GetYaxis()->SetLabelSize(0.040);
        graph->GetYaxis()->SetTitleOffset(1.35);

        if (fitLinear && x.size() >= 3) {
            TF1* line = new TF1(Form("line_%s_paddle%02d", plane.Data(), paddle),
                                "pol1", -0.5, 1.0);
            lineFits.push_back(line);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            graph->Fit(line, "QRN");
            line->Draw("same");
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
            runLabel.DrawLatex(x[i], y[i] + 0.035 * ySpan, Form("%d", runs[i]));
        }

        TLatex note;
        note.SetNDC(kTRUE);
        note.SetTextSize(0.028);
        note.DrawLatex(0.16, 0.91, "Point labels are run numbers");
    }

    canvas.cd();
    const TString pngName = Form("%s/%s_QDC_peak_separation_vs_HV.png",
                                 outputDirectory.Data(), plane.Data());
    const TString pdfName = Form("%s/%s_QDC_peak_separation_vs_HV.pdf",
                                 outputDirectory.Data(), plane.Data());
    canvas.SaveAs(pngName.Data());
    canvas.SaveAs(pdfName.Data());

    // The canvas/pads retain the drawn ROOT objects until this function returns.
    // Their total number is small (at most eight graphs and eight optional lines).
}

} // namespace BH_QDC_HVScan_Detail

void BH_QDC_HVScan(const char* inputDir = ".",
                   const char* outputDir = "BH_QDC_HV_scan",
                   int defaultRun = 33566,
                   bool saveDiagnostics = true,
                   bool fitLinear = false)
{
    using namespace BH_QDC_HVScan_Detail;

    gROOT->SetBatch(kTRUE);
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
                std::cout << "pedestal = " << std::fixed << std::setprecision(2)
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

    if (saveDiagnostics) {
        std::cout << "Fit diagnostics:          " << diagnosticsPdf << "\n";
    }

    std::cout << std::endl;
}
