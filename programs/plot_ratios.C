#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

// ------------------------------------------------------------
// Constants: use GeV-based natural units.
// Cross section units cancel in the mu/e ratio.
// ------------------------------------------------------------
const double alpha = 1.0 / 137.035999084;     // fine structure constant
const double Mp    = 0.9382720813;            // proton mass [GeV]
const double me    = 0.00051099895;           // electron mass [GeV]
const double mmu   = 0.1056583755;            // muon mass [GeV]
const double mup   = 2.79284734463;           // proton magnetic moment
const double hbarc = 0.1973269804;            // [GeV fm]

// Convert Q^2 from GeV^2 to fm^-2 for the Borkowski fits.
double GeV2_to_fm2(double Q2_GeV2)
{
    return Q2_GeV2 / (hbarc * hbarc);
}

// ------------------------------------------------------------
// Helper: single-pole sum used by Borkowski-style fits.
// G(q^2) = sum_i a_i m_i^2 / (m_i^2 + q^2)
// with q^2 and m_i^2 in fm^-2.
// ------------------------------------------------------------
double PoleSum4(double q2_fm2, const double a[4], const double m2[4])
{
    double sum = 0.0;
    for (int i = 0; i < 4; i++) {
        sum += a[i] * m2[i] / (m2[i] + q2_fm2);
    }
    return sum;
}

// ------------------------------------------------------------
// 1. Borkowski 1974
// Returns GEp and GMp, with GMp normalized so GMp(0)=mu_p.
// ------------------------------------------------------------
void FF_Borkowski1974(double Q2, double &GE, double &GM)
{
    const double q2 = GeV2_to_fm2(Q2);

    const double aE[4]  = { 0.250,  1.230, -0.442, -0.038 };
    const double m2E[4] = { 3.53,  15.02,  44.08, 154.2  };

    const double aM[4]  = { 0.843,  0.522, -0.366,  0.001 };
    const double m2M[4] = { 8.73,  15.02,  44.08, 355.4  };

    GE = PoleSum4(q2, aE, m2E);
    GM = mup * PoleSum4(q2, aM, m2M);
}

// ------------------------------------------------------------
// 2. Borkowski 1975
// Returns GEp and GMp, with GMp normalized so GMp(0)=mu_p.
// ------------------------------------------------------------
void FF_Borkowski1975(double Q2, double &GE, double &GM)
{
    const double q2 = GeV2_to_fm2(Q2);

    const double aE[4]  = { 0.219,  1.371, -0.634,  0.044 };
    const double m2E[4] = { 353.0, 15.02,  44.08, 154.20 };

    const double aM[4]  = { 0.794,  0.594, -0.393,  0.005 };
    const double m2M[4] = { 873.0, 15.02,  44.08, 355.40 };

    GE = PoleSum4(q2, aE, m2E);
    GM = mup * PoleSum4(q2, aM, m2M);
}

// ------------------------------------------------------------
// 3. Bosted 1995
// Bosted uses Q = sqrt(Q^2) in GeV.
// ------------------------------------------------------------
void FF_Bosted1995(double Q2, double &GE, double &GM)
{
    const double Q = std::sqrt(Q2);

    GE = 1.0 / (
        1.0
        + 0.62 * Q
        + 0.68 * std::pow(Q, 2)
        + 2.80 * std::pow(Q, 3)
        + 0.83 * std::pow(Q, 4)
    );

    const double GM_over_mup = 1.0 / (
        1.0
        + 0.35 * Q
        + 2.44 * std::pow(Q, 2)
        + 0.50 * std::pow(Q, 3)
        + 1.04 * std::pow(Q, 4)
        + 0.34 * std::pow(Q, 5)
    );

    GM = mup * GM_over_mup;
}

// ------------------------------------------------------------
// 4. Arrington 2004 Rosenbluth form factors
// G = 1 / (1 + p2 Q^2 + p4 Q^4 + ...)
// Q^2 in GeV^2.
// ------------------------------------------------------------
void FF_Arrington2004(double Q2, double &GE, double &GM)
{
    const double Q4  = Q2 * Q2;
    const double Q6  = Q4 * Q2;
    const double Q8  = Q4 * Q4;
    const double Q10 = Q8 * Q2;
    const double Q12 = Q6 * Q6;

    GE = 1.0 / (
        1.0
        + 3.226  * Q2
        + 1.508  * Q4
        - 0.3773 * Q6
        + 0.611  * Q8
        - 0.1853 * Q10
        + 1.596e-2 * Q12
    );

    const double GM_over_mup = 1.0 / (
        1.0
        + 3.19      * Q2
        + 1.355     * Q4
        + 0.151     * Q6
        - 1.143e-2  * Q8
        + 5.333e-4  * Q10
        - 9.003e-6  * Q12
    );

    GM = mup * GM_over_mup;
}

// ------------------------------------------------------------
// 5. Kelly 2004
// G = (1 + a1 tau) / (1 + b1 tau + b2 tau^2 + b3 tau^3)
// tau = Q^2 / (4 Mp^2).
// ------------------------------------------------------------
void FF_Kelly2004(double Q2, double &GE, double &GM)
{
    const double tau = Q2 / (4.0 * Mp * Mp);

    GE = (1.0 - 0.24 * tau) /
         (1.0 + 10.98 * tau + 12.82 * tau * tau + 21.97 * tau * tau * tau);

    const double GM_over_mup =
        (1.0 + 0.12 * tau) /
        (1.0 + 10.97 * tau + 18.86 * tau * tau + 6.55 * tau * tau * tau);

    GM = mup * GM_over_mup;
}

// ------------------------------------------------------------
// 6. Arrington, Melnitchouk, Tjon 2007
// G = (1 + a1 tau + a2 tau^2 + a3 tau^3)
//     / (1 + b1 tau + ... + b5 tau^5)
// ------------------------------------------------------------
void FF_AMT2007(double Q2, double &GE, double &GM)
{
    const double t  = Q2 / (4.0 * Mp * Mp);
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t2 * t2;
    const double t5 = t4 * t;

    const double GE_num = 1.0 + 3.439 * t - 1.602 * t2 + 0.068 * t3;
    const double GE_den = 1.0 + 15.055 * t + 48.061 * t2 + 99.304 * t3
                              + 0.012 * t4 + 8.650 * t5;

    const double GM_num = 1.0 - 1.465 * t + 1.260 * t2 + 0.262 * t3;
    const double GM_den = 1.0 + 9.627 * t + 0.000 * t2 + 0.000 * t3
                              + 11.179 * t4 + 13.245 * t5;

    GE = GE_num / GE_den;
    GM = mup * (GM_num / GM_den);
}

// ------------------------------------------------------------
// Elastic kinematics at fixed incoming lepton momentum.
// Solve for outgoing lepton momentum p' from energy conservation:
// E + Mp = E' + Ep_recoil.
// ------------------------------------------------------------
double SolvePprime(double p, double theta, double m)
{
    const double E = std::sqrt(p * p + m * m);
    const double c = std::cos(theta);

    auto f = [&](double pp) {
        const double EpLep = std::sqrt(pp * pp + m * m);
        const double pRec2 = p * p + pp * pp - 2.0 * p * pp * c;
        const double EpRec = std::sqrt(Mp * Mp + pRec2);
        return E + Mp - EpLep - EpRec;
    };

    double lo = 0.0;
    double hi = p;

    double flo = f(lo);
    double fhi = f(hi);

    if (flo * fhi > 0.0) {
        return NAN;
    }

    for (int i = 0; i < 100; i++) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = f(mid);

        if (flo * fmid > 0.0) {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
            fhi = fmid;
        }
    }

    return 0.5 * (lo + hi);
}

// ------------------------------------------------------------
// Preedom-Tegen massive-lepton unpolarized cross section.
// pBeam is in GeV/c.
// theta is in radians.
// ff is one of the form factor functions above.
// ------------------------------------------------------------
typedef void (*FFunc)(double Q2, double &GE, double &GM);

double CrossSection(double pBeam, double theta, double m, FFunc ff)
{
    const double E  = std::sqrt(pBeam * pBeam + m * m);
    const double pp = SolvePprime(pBeam, theta, m);

    if (!std::isfinite(pp)) return NAN;

    const double Ep = std::sqrt(pp * pp + m * m);

    // Positive Q^2 = -q^2. For elastic scattering from a proton at rest:
    const double Q2 = 2.0 * Mp * (E - Ep);

    if (Q2 <= 0.0) return NAN;

    const double sinHalf = std::sin(theta / 2.0);
    const double sinHalf2 = sinHalf * sinHalf;

    const double x = Q2 / (4.0 * E * Ep);  // x = -q^2 / (4 E E')
    if (x <= 0.0 || x >= 1.0) return NAN;

    const double d = std::sqrt(
        (1.0 - (m * m) / (E * E)) /
        (1.0 - (m * m) / (Ep * Ep))
    );

    double GE, GM;
    ff(Q2, GE, GM);

    const double eta = Q2 / (4.0 * Mp * Mp);

    const double sigma_ns =
        (alpha * alpha / (4.0 * E * E))
        * (1.0 - x) / (x * x)
        * (1.0 / d)
        / (1.0 + (2.0 * E * d / Mp) * sinHalf2 + (E / Mp) * (1.0 - d));

    const double R =
        (GE * GE + eta * GM * GM) / (1.0 + eta)
        + (2.0 * eta - (m * m) / (Mp * Mp)) * GM * GM * x / (1.0 - x);

    return sigma_ns * R;
}

double MuOverElectronRatio(double pBeam, double theta, FFunc ff)
{
    const double sigMu = CrossSection(pBeam, theta, mmu, ff);
    const double sigEl = CrossSection(pBeam, theta, me,  ff);

    if (!std::isfinite(sigMu) || !std::isfinite(sigEl) || sigEl == 0.0) {
        return NAN;
    }

    return sigMu / sigEl;
}

// ------------------------------------------------------------
// Main ROOT function.
// Creates three plots:
//   ratio_p210.png
//   ratio_p160.png
//   ratio_p115.png
// ------------------------------------------------------------
void plot_ratios()
{
    gStyle->SetOptStat(0);

    struct Model {
        const char *name;
        FFunc ff;
        int color;
        int style;
    };

    std::vector<Model> models = {
        {"Borkowski 1974", FF_Borkowski1974, kRed + 1,     1},
        {"Borkowski 1975", FF_Borkowski1975, kOrange + 7,  2},
        {"Bosted 1995",    FF_Bosted1995,    kGreen + 2,   3},
        {"Arrington 2004", FF_Arrington2004, kBlue + 1,    4},
        {"Kelly 2004",     FF_Kelly2004,     kViolet + 1,  5},
        {"AMT 2007",       FF_AMT2007,       kBlack,       6}
    };

    // Beam momenta in GeV/c.
    std::vector<double> momenta = {0.210, 0.160, 0.115};

    // Avoid exactly 0 and 180 degrees.
    const int nPts = 171;
    const double thetaMinDeg = 5.0;
    const double thetaMaxDeg = 102.0;

    for (double pBeam : momenta) {
        TCanvas *c = new TCanvas(
            Form("c_p%.0f", pBeam * 1000.0),
            Form("p = %.0f MeV/c", pBeam * 1000.0),
            1000, 700
        );

        TMultiGraph *mg = new TMultiGraph();
        TLegend *leg = new TLegend(0.14, 0.62, 0.42, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        double maxSpreadPercent = 0.0;
        double angleAtMaxSpread = 0.0;

        for (const auto &model : models) {
            std::vector<double> x;
            std::vector<double> y;

            for (int i = 0; i < nPts; i++) {
                const double thetaDeg =
                    thetaMinDeg
                    + (thetaMaxDeg - thetaMinDeg) * i / double(nPts - 1);

                const double thetaRad = thetaDeg * TMath::DegToRad();

                const double ratio =
                    MuOverElectronRatio(pBeam, thetaRad, model.ff);

                if (std::isfinite(ratio)) {
                    x.push_back(thetaDeg);
                    y.push_back(ratio);
                }
            }

            TGraph *gr = new TGraph((int)x.size(), x.data(), y.data());
            gr->SetLineColor(model.color);
            gr->SetLineStyle(model.style);
            gr->SetLineWidth(3);

            mg->Add(gr, "L");
            leg->AddEntry(gr, model.name, "l");
        }

        // Also compute the maximum spread among form-factor choices.
        for (int i = 0; i < nPts; i++) {
            const double thetaDeg =
                thetaMinDeg
                + (thetaMaxDeg - thetaMinDeg) * i / double(nPts - 1);

            const double thetaRad = thetaDeg * TMath::DegToRad();

            std::vector<double> ratios;
            for (const auto &model : models) {
                const double r = MuOverElectronRatio(pBeam, thetaRad, model.ff);
                if (std::isfinite(r)) ratios.push_back(r);
            }

            if (ratios.size() >= 2) {
                const auto [rminIt, rmaxIt] =
                    std::minmax_element(ratios.begin(), ratios.end());

                const double rmin = *rminIt;
                const double rmax = *rmaxIt;
                const double mean = 0.5 * (rmin + rmax);

                const double spreadPercent = 100.0 * (rmax - rmin) / mean;

                if (spreadPercent > maxSpreadPercent) {
                    maxSpreadPercent = spreadPercent;
                    angleAtMaxSpread = thetaDeg;
                }
            }
        }

        mg->SetTitle(Form(
            "#mu p / e p cross-section ratio at p = %.0f MeV/c;#theta_{lab} [deg];(d#sigma_{#mu p}/d#Omega)/(d#sigma_{e p}/d#Omega)",
            pBeam * 1000.0
        ));

        mg->Draw("A");
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();

        leg->Draw();

        c->SetGrid();
        c->SaveAs(Form("ratio_p%.0f.png", pBeam * 1000.0));

        std::cout << std::fixed << std::setprecision(3)
                  << "p = " << pBeam * 1000.0 << " MeV/c: "
                  << "max model spread = " << maxSpreadPercent
                  << "% near theta = " << angleAtMaxSpread
                  << " deg" << std::endl;
    }
}