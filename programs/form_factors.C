// form_factors.C
// ROOT macro for comparing proton form-factor parametrizations.
//
// The code keeps each paper's parametrization in its own clearly labeled
// function block.  Convention in this file:
//   Q2/qsq is positive Q^2 in GeV^2.
//   GEp is normalized so GEp(0) = 1.
//   GMp functions return the physical Sachs GMp, so GMp(0) = mu_p.

#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TAxis.h"

#include <vector>

constexpr double kHbarC = 0.197327;   // GeV*fm
constexpr double kMp    = 0.938272;   // proton mass, GeV
constexpr double kMuP   = 2.793;      // proton magnetic moment in nuclear magnetons

// -----------------------------------------------------------------------------
// Baseline dipole form factor
// -----------------------------------------------------------------------------
double gdipole(double qsq)
{
    // Standard dipole, Q^2 in GeV^2.
    return TMath::Power(1.0 + qsq/0.71, -2.0);
}

// -----------------------------------------------------------------------------
// Shared helper for the Borkowski pole fits
// -----------------------------------------------------------------------------
double qsq_gev2_to_fm2(double qsq)
{
    // Borkowski uses q^2 in fm^-2.  Since hbar*c = 0.197327 GeV*fm,
    // q^2[fm^-2] = Q^2[GeV^2] / (hbar*c)^2.
    return qsq / (kHbarC*kHbarC);
}

double borkowski_pole_sum(double qsq, const double a[4], const double m2[4])
{
    const double q2fm = qsq_gev2_to_fm2(qsq);
    double g = 0.0;
    for (int i = 0; i < 4; ++i) {
        g += a[i] * m2[i] / (m2[i] + q2fm);
    }
    return g;
}

// -----------------------------------------------------------------------------
// Borkowski et al. 1974, Nucl. Phys. A222, 269
// Four-pole fit; GEp(0)=1 and GMp(0)/mu_p=1.
// -----------------------------------------------------------------------------
double gep_borkowski1974(double qsq)
{
    const double a[4]  = { 0.250,  1.230, -0.442, -0.038 };
    const double m2[4] = { 3.530, 15.020, 44.080, 154.200 };
    return borkowski_pole_sum(qsq, a, m2);
}

double gmp_borkowski1974(double qsq)
{
    const double a[4]  = { 0.843,  0.522, -0.366,  0.001 };
    const double m2[4] = { 8.730, 15.020, 44.080, 355.400 };
    return kMuP * borkowski_pole_sum(qsq, a, m2);
}

// -----------------------------------------------------------------------------
// Borkowski et al. 1975, Nucl. Phys. B93, 461
// Updated four-pole fit; GEp(0)=1 and GMp(0)/mu_p=1.
// -----------------------------------------------------------------------------
double gep_borkowski1975(double qsq)
{
    const double a[4]  = { 0.219,  1.371, -0.634,  0.044 };
    const double m2[4] = { 3.530, 15.020, 44.080, 154.200 };
    return borkowski_pole_sum(qsq, a, m2);
}

double gmp_borkowski1975(double qsq)
{
    const double a[4]  = { 0.794,  0.594, -0.393,  0.005 };
    const double m2[4] = { 8.730, 15.020, 44.080, 355.400 };
    return kMuP * borkowski_pole_sum(qsq, a, m2);
}

// -----------------------------------------------------------------------------
// Bosted 1995, Phys. Rev. C 51, 409
// Polynomial in Q = sqrt(Q^2), with Q in GeV.
// -----------------------------------------------------------------------------
double gep_bosted1995(double qsq)
{
    const double q = TMath::Sqrt(qsq);
    return 1.0 / (1.0 + q*(0.62 + q*(0.68 + q*(2.80 + q*0.83))));
}

double gmp_bosted1995(double qsq)
{
    const double q = TMath::Sqrt(qsq);
    return kMuP / (1.0 + q*(0.35 + q*(2.44 + q*(0.50 + q*(1.04 + q*0.34)))));
}

// -----------------------------------------------------------------------------
// Arrington 2004, Phys. Rev. C 69, 022201(R)
// Rosenbluth fit.  Polynomial in Q^2, Q^2 in GeV^2.
// -----------------------------------------------------------------------------
double gep_arrington2004_rosenbluth(double qsq)
{
    const double den = 1.0 + qsq*(3.226 + qsq*(1.508 + qsq*(-0.3773 + qsq*(0.611 + qsq*(-0.1853 + qsq*0.01596)))));
    return 1.0 / den;
}

double gmp_arrington2004_rosenbluth(double qsq)
{
    const double den = 1.0 + qsq*(3.190 + qsq*(1.355 + qsq*(0.151 + qsq*(-0.0114 + qsq*(0.000533 + qsq*(-0.00000900))))));
    return kMuP / den;
}

// -----------------------------------------------------------------------------
// Arrington 2004, Phys. Rev. C 69, 022201(R)
// Polarization-transfer-compatible fit.  Polynomial in Q^2, Q^2 in GeV^2.
// -----------------------------------------------------------------------------
double gep_arrington2004_polarization(double qsq)
{
    const double den = 1.0 + qsq*(2.940 + qsq*(3.040 + qsq*(-2.255 + qsq*(2.002 + qsq*(-0.5338 + qsq*0.04875)))));
    return 1.0 / den;
}

double gmp_arrington2004_polarization(double qsq)
{
    const double den = 1.0 + qsq*(3.000 + qsq*(1.390 + qsq*(0.122 + qsq*(-0.00834 + qsq*(0.000425 + qsq*(-0.00000779))))));
    return kMuP / den;
}

// -----------------------------------------------------------------------------
// Kelly 2004, Phys. Rev. C 70, 068202
// Rational function in tau = Q^2/(4 Mp^2).
// -----------------------------------------------------------------------------
double gep_kelly2004(double tau)
{
    const double num = 1.0 - 0.24*tau;
    const double den = 1.0 + tau*(10.98 + tau*(12.82 + tau*21.97));
    return num / den;
}

double gmp_kelly2004(double tau)
{
    const double num = kMuP * (1.0 + 0.12*tau);
    const double den = 1.0 + tau*(10.97 + tau*(18.86 + tau*6.55));
    return num / den;
}

// -----------------------------------------------------------------------------
// Arrington, Melnitchouk, Tjon 2007, Phys. Rev. C 76, 035205
// TPE-corrected rational function in tau = Q^2/(4 Mp^2).
// -----------------------------------------------------------------------------
double gep_amt2007(double tau)
{
    const double num = 1.0 + tau*(3.439 + tau*(-1.602 + tau*0.068));
    const double den = 1.0 + tau*(15.055 + tau*(48.061 + tau*(99.304 + tau*(0.012 + tau*8.650))));
    return num / den;
}

double gmp_amt2007(double tau)
{
    const double num = kMuP * (1.0 + tau*(-1.465 + tau*(1.260 + tau*0.262)));
    const double den = 1.0 + tau*(9.627 + tau*(0.000 + tau*(0.000 + tau*(11.179 + tau*13.245))));
    return num / den;
}

// -----------------------------------------------------------------------------
// Small drawing helper
// -----------------------------------------------------------------------------
void add_curve(TMultiGraph* mg,
               TLegend* legend,
               const int n,
               const double* x,
               const double* y,
               const int color,
               const int width,
               const char* label)
{
    TGraph* gr = new TGraph(n, x, y);
    gr->SetLineColor(color);
    gr->SetLineWidth(width);
    mg->Add(gr, "L");
    legend->AddEntry(gr, label, "l");
}

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------
void form_factors()
{
    gROOT->Reset();
    gSystem->Load("libPhysics");

    const int n = 101;
    double q2[n];

    // GE/GD arrays
    double ge_b74[n], ge_b75[n], ge_bosted[n], ge_arr_r[n], ge_arr_p[n], ge_kelly[n], ge_amt[n];

    // GM/(mu_p GD) arrays
    double gm_b74[n], gm_b75[n], gm_bosted[n], gm_arr_r[n], gm_arr_p[n], gm_kelly[n], gm_amt[n];

    for (int iq = 0; iq < n; ++iq) {
        const double qsq = 0.01 * iq;       // Q^2 in GeV^2, from 0 to 1 GeV^2
        const double tau = qsq/(4.0*kMp*kMp);
        const double gd  = gdipole(qsq);

        q2[iq] = qsq;

        ge_b74[iq]    = gep_borkowski1974(qsq)              / gd;
        ge_b75[iq]    = gep_borkowski1975(qsq)              / gd;
        ge_bosted[iq] = gep_bosted1995(qsq)                 / gd;
        ge_arr_r[iq]  = gep_arrington2004_rosenbluth(qsq)   / gd;
        ge_arr_p[iq]  = gep_arrington2004_polarization(qsq) / gd;
        ge_kelly[iq]  = gep_kelly2004(tau)                  / gd;
        ge_amt[iq]    = gep_amt2007(tau)                    / gd;

        gm_b74[iq]    = gmp_borkowski1974(qsq)              / (kMuP*gd);
        gm_b75[iq]    = gmp_borkowski1975(qsq)              / (kMuP*gd);
        gm_bosted[iq] = gmp_bosted1995(qsq)                 / (kMuP*gd);
        gm_arr_r[iq]  = gmp_arrington2004_rosenbluth(qsq)   / (kMuP*gd);
        gm_arr_p[iq]  = gmp_arrington2004_polarization(qsq) / (kMuP*gd);
        gm_kelly[iq]  = gmp_kelly2004(tau)                  / (kMuP*gd);
        gm_amt[iq]    = gmp_amt2007(tau)                    / (kMuP*gd);
    }

    // Canvas 1: electric form factors
    TCanvas* c1 = new TCanvas("GE_wrt_dipole", "GE/GD comparison", 1200, 800);
    TMultiGraph* ge_mg = new TMultiGraph();
    TLegend* ge_leg = new TLegend(0.62, 0.58, 0.90, 0.90);
    ge_leg->SetHeader("Form-factor paper", "C");
    ge_leg->SetFillColor(kWhite);
    ge_leg->SetBorderSize(1);
    ge_leg->SetTextSize(0.032);

    add_curve(ge_mg, ge_leg, n, q2, ge_b74,    kOrange+7, 2, "Borkowski 1974");
    add_curve(ge_mg, ge_leg, n, q2, ge_b75,    kOrange+1, 2, "Borkowski 1975");
    add_curve(ge_mg, ge_leg, n, q2, ge_bosted, kGreen+2,  2, "Bosted 1995");
    add_curve(ge_mg, ge_leg, n, q2, ge_arr_r,  kBlue+1,   2, "Arrington 2004 Ros.");
    add_curve(ge_mg, ge_leg, n, q2, ge_arr_p,  kCyan+2,   2, "Arrington 2004 Pol.");
    add_curve(ge_mg, ge_leg, n, q2, ge_kelly,  kRed+1,    2, "Kelly 2004");
    add_curve(ge_mg, ge_leg, n, q2, ge_amt,    kMagenta+1,2, "AMT 2007");

    ge_mg->SetTitle("Proton electric form factor comparisons");
    ge_mg->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    ge_mg->GetYaxis()->SetTitle("G_{E}^{p} / G_{D}");
    ge_mg->GetXaxis()->CenterTitle();
    ge_mg->GetYaxis()->CenterTitle();
    ge_mg->GetYaxis()->SetTitleOffset(1.2);
    ge_mg->Draw("AL");
    ge_leg->Draw();
    c1->Modified();
    c1->Update();

    // Canvas 2: magnetic form factors
    TCanvas* c2 = new TCanvas("GM_wrt_dipole", "GM/(mu_p GD) comparison", 1200, 800);
    TMultiGraph* gm_mg = new TMultiGraph();
    TLegend* gm_leg = new TLegend(0.62, 0.58, 0.90, 0.90);
    gm_leg->SetHeader("Form-factor paper", "C");
    gm_leg->SetFillColor(kWhite);
    gm_leg->SetBorderSize(1);
    gm_leg->SetTextSize(0.032);

    add_curve(gm_mg, gm_leg, n, q2, gm_b74,    kOrange+7, 2, "Borkowski 1974");
    add_curve(gm_mg, gm_leg, n, q2, gm_b75,    kOrange+1, 2, "Borkowski 1975");
    add_curve(gm_mg, gm_leg, n, q2, gm_bosted, kGreen+2,  2, "Bosted 1995");
    add_curve(gm_mg, gm_leg, n, q2, gm_arr_r,  kBlue+1,   2, "Arrington 2004 Ros.");
    add_curve(gm_mg, gm_leg, n, q2, gm_arr_p,  kCyan+2,   2, "Arrington 2004 Pol.");
    add_curve(gm_mg, gm_leg, n, q2, gm_kelly,  kRed+1,    2, "Kelly 2004");
    add_curve(gm_mg, gm_leg, n, q2, gm_amt,    kMagenta+1,2, "AMT 2007");

    gm_mg->SetTitle("Proton magnetic form factor comparisons");
    gm_mg->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    gm_mg->GetYaxis()->SetTitle("G_{M}^{p} / (#mu_{p} G_{D})");
    gm_mg->GetXaxis()->CenterTitle();
    gm_mg->GetYaxis()->CenterTitle();
    gm_mg->GetYaxis()->SetTitleOffset(1.2);
    gm_mg->Draw("AL");
    gm_leg->Draw();
    c2->Modified();
    c2->Update();
}

// Keep the original macro name available too.
void gratios()
{
    form_factors();
}
