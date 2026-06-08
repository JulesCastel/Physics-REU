#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <string>

long double ep_cross_section_preedom_tegen_p(
    long double p_beam_MeVc = 200.0L,  // MUSE beam momentum in MeV/c.
    long double theta_deg = 30.0L,    // Scattering angle in degrees, should be between 20 and 100 degrees.
    std::string lepton = "electron"    // Choose "electron" or "muon".
) {
    const long double alpha = 1.0L / 137.035999084L;       // Fine-structure constant alpha.

    const long double electron_mass = 0.0L;   // Electron rest mass in GeV/c^2.

    const long double muon_mass = 0.1056583755L;           // Muon rest mass in GeV/c^2.

    const long double Mp = 0.93827208816L;                 // Proton rest mass in GeV/c^2.

    const long double mu_p = 2.79284734463L;               // Proton magnetic moment in nuclear magnetons; used for GM = mu_p * GE.

    const long double dipole_scale = 0.71L;                // Dipole form-factor scale in GeV^2.

    const long double gev2_to_mb = 0.389379338L;           // Unit conversion: 1 GeV^(-2) = 0.389379338 mb.

    const long double mev_to_gev = 0.001L;                 // Converts MeV to GeV.

    const long double pi = acosl(-1.0L);                   // Pi, used for converting degrees to radians.

    long double ml = electron_mass;                        // Lepton rest mass m in GeV/c^2; defaults to electron.

    if (lepton == "electron" || lepton == "e" || lepton == "Electron") {
        ml = electron_mass;
    }
    else if (lepton == "muon" || lepton == "mu" || lepton == "Muon") {
        ml = muon_mass;
    }
    else {
        std::cout << "Unknown lepton type: " << lepton << "\n";
        std::cout << "Use \"electron\", \"e\", \"muon\", or \"mu\".\n";
        return std::numeric_limits<long double>::quiet_NaN();
    }

    const long double p = p_beam_MeVc * mev_to_gev;        // Incoming lepton momentum magnitude |k| in GeV/c.

    const long double E = sqrtl(p * p + ml * ml);          // Incoming lepton total energy E in GeV.

    const long double theta = theta_deg * pi / 180.0L;     // Scattering angle theta in radians.

    const long double c = cosl(theta);                     // cos(theta).

    const long double s2 = sinl(theta / 2.0L);             // sin(theta/2).

    const long double k = p;                               // Incoming lepton 3-momentum magnitude |k| in GeV/c.

    // Solve elastic kinematics for outgoing total lepton energy E'.
    // Uses q^2 = 2M(E' - E) and q^2 = 2m^2 - 2EE' + 2|k||k'|cos(theta).
    const long double B = Mp + E;                          // Algebraic coefficient from elastic energy-momentum conservation.

    const long double C = Mp * E + ml * ml;                // Algebraic coefficient from elastic energy-momentum conservation.

    const long double A_quad = B * B - k * k * c * c;      // Quadratic coefficient multiplying E'^2.

    const long double B_quad = -2.0L * B * C;              // Quadratic coefficient multiplying E'.

    const long double C_quad = C * C + k * k * c * c * ml * ml; // Constant term in the E' quadratic.

    const long double disc = B_quad * B_quad - 4.0L * A_quad * C_quad; // Quadratic discriminant.

    if (disc < 0.0L) {
        std::cout << "No physical elastic-scattering solution.\n";
        return std::numeric_limits<long double>::quiet_NaN();
    }

    const long double Eprime_1 = (-B_quad + sqrtl(disc)) / (2.0L * A_quad); // First possible outgoing lepton energy in GeV.

    const long double Eprime_2 = (-B_quad - sqrtl(disc)) / (2.0L * A_quad); // Second possible outgoing lepton energy in GeV.

    long double Eprime = Eprime_1;                         // Outgoing lepton total energy E' in GeV.

    if (Eprime_1 > E || Eprime_1 < ml) Eprime = Eprime_2;  // Choose physical root: ml <= E' <= E.

    if (Eprime > E || Eprime < ml) {
        std::cout << "No physical elastic-scattering root was found.\n";
        return std::numeric_limits<long double>::quiet_NaN();
    }

    const long double kprime = sqrtl(Eprime * Eprime - ml * ml); // Outgoing lepton 3-momentum magnitude |k'| in GeV/c.

    // Eq. (3): d factor.
    const long double d = sqrtl(
        (1.0L - ml * ml / (E * E)) /
        (1.0L - ml * ml / (Eprime * Eprime))
    );

    const long double beta_in = sqrtl(1.0L - ml * ml / (E * E));             // Incoming lepton velocity factor |k|/E.

    const long double beta_out = sqrtl(1.0L - ml * ml / (Eprime * Eprime)); // Outgoing lepton velocity factor |k'|/E'.

    // Eq. (3): invariant four-momentum transfer q^2 in GeV^2.
    const long double q2 =
        -4.0L * E * Eprime * s2 * s2 * beta_in * beta_out
        + 2.0L * ml * ml
        - 2.0L * E * Eprime * (1.0L - beta_in * beta_out);

    const long double Q2 = -q2;                             // Positive momentum transfer Q^2 = -q^2 in GeV^2.

    const long double eta = Q2 / (4.0L * Mp * Mp);          // eta = -q^2 / (4M^2), dimensionless.

    const long double x = Q2 / (4.0L * E * Eprime);         // x = -q^2 / (4EE'), dimensionless.

    // Dipole Sachs form factors used by Preedom & Tegen for the proton.
    const long double GE = powl(1.0L + Q2 / dipole_scale, -2.0L); // Proton electric Sachs form factor G_E(q^2).

    const long double GM = mu_p * GE;                       // Proton magnetic Sachs form factor G_M(q^2).

    // Eq. (4): structure factor R.
    const long double R =
        (GE * GE + eta * GM * GM) / (1.0L + eta)
        + (2.0L * eta - ml * ml / (Mp * Mp)) * GM * GM * (x / (1.0L - x));

    // Eq. (2): no-structure cross section.
    const long double sigma_ns =
        (alpha * alpha / (4.0L * E * E))
        * ((1.0L + x) / (x * x))
        * (1.0L / d)
        / (
            1.0L
            + (2.0L * E * d / Mp) * s2 * s2
            + (E / Mp) * (1.0L - d)
        );

    // Eq. (2): full differential cross section dσ/dΩ.
    const long double dsigma_dOmega_GeV2 = sigma_ns * R;    // Cross section in GeV^(-2)/sr.

    const long double dsigma_dOmega_mb = dsigma_dOmega_GeV2 * gev2_to_mb; // Cross section in mb/sr.

    std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10);

    std::cout << "lepton = " << lepton << "\n";

    std::cout << "p_beam = " << p_beam_MeVc << " MeV/c\n";

    std::cout << "m_lepton = " << ml << " GeV/c^2\n";

    std::cout << "E = " << E << " GeV\n";

    std::cout << "E' = " << Eprime << " GeV\n";

    std::cout << "p' = " << kprime / mev_to_gev << " MeV/c\n";

    std::cout << "Q^2 = " << Q2 << " GeV^2\n";

    std::cout << "dσ/dΩ = " << dsigma_dOmega_mb << " mb/sr\n";

    return dsigma_dOmega_mb;
}