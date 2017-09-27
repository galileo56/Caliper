
#include <iostream>
#include <string>
#include "radRet.h"

using namespace std;
// using namespace radiative;

int main () //can I comment here?
{
    //Default options of the xsection_calculator
    double photon_energy, x; // fraction of c.o.m. energy carried away by the photon
    int nf = 6; // number of flavors for Top quark
    int order = 1, ordNRQCD = 3; //Only 1 loop implemented at fixed order, N3LL for NRQCD
    int orderAlp = 4, runAlp = 4, runMass = 4, ordMass = 4, ord1S = 4; // running and matching for alphaS and masses
    double R1S = 35., lambda = 1.;// numbers related to the MSRn mass
    double v1 = 0.1, v2 = 0.5;    // parameters related to the matching procedure
    //string scheme = "MSbar";
    double energy_cm = 500; // c.o.m. energy
    double mZ = 91.187, gammaZ = 2.5042, sin2ThetaW = 0.2312; // ElectroWeak parameters
    double mT = 160., mB = 4.2, mC = 1.3; // quark masses
    double gt = 1.5; // top width
    double eH = 1, hnu = 1; // ratio between renormalization scales and energy of the hadronic system
    double alphaQED = 1/127.925; // e.m. alpha at the Z pole
    double alphamZ = 0.1181;     // strong coupling constant at the Z pole
    //Theta is defined from 7º to 90º (excluding the beampipe), Phi is fully integrated
    double angle = 3.141592 * 48.5 / 180; //Pointing angle (theta) of the cone 48.5º
    double half_aperture =  3.141592 * 41.5 / 180;//Half aperture of the cone 41.5º

    cout << "Fixed Order distribution, MS-bar scheme, one loop \n\n";

    cout << "x, E_photon (GeV), E_hadrons (GeV), differential, cone \n\n";

    for (int i = 1; i < 52; i++) {
    photon_energy = energy_cm * 0.018 + (i - 1)/25.;
    x = photon_energy/energy_cm;

    cout << x << "  " << photon_energy << "  " << energy_cm * (1 - 2 * x) << " "
    << " " << sigmamassrad("MSbar", "all", orderAlp, runAlp, runMass, order, nf,
    mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED, mT, mT, mB, mB, mC, mC, eH,
    energy_cm, x, angle) <<
    "  " << sigmamassradcone("MSbar", "all", orderAlp, runAlp, runMass, order,
    nf, mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED, mT, mT, mB, mB, mC, mC, eH,
    energy_cm, x, angle, half_aperture) << "\n\n" ;

    }

    cout << "\n integration over the cone and E_photon between 9 and 11 GeV \n\n";

   cout << sigmamassradconecum("MSbar", "all", orderAlp, runAlp, runMass, order,
   nf, mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED, mT, mT, mB, mB, mC, mC, eH,
   energy_cm, 0.018, 0.022, angle, half_aperture) << "\n";

    cout << "NRQCD matched to Fixed Order distribution, 1S-MSRb scheme, highest order \n\n";
    for (int i = 1; i < 52; i++) {
    photon_energy = energy_cm * 0.018 + (i - 1)/25.;
    x = photon_energy/energy_cm;

    cout << x << "  " << photon_energy << "  " << energy_cm * (1 - 2 * x) << " "
    << " " << sigmamatchedrad("MSRb", runAlp, runMass, ordMass, ordNRQCD,
    ord1S, R1S, "numeric", lambda, gt, mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED,
    mT, eH, hnu, v1, v2, energy_cm, x, angle) <<
    "  " << sigmamatchedradcone("MSRb", runAlp, runMass, ordMass, ordNRQCD,
    ord1S, R1S, "numeric", lambda, gt, mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED,
    mT, eH, hnu, v1, v2, energy_cm, x, angle, half_aperture) << "\n\n" ;

    }

    cout << "\n integration over the cone and E_photon between 9 and 11 GeV \n\n";

   cout << sigmamatchedradconecum("MSRb", runAlp, runMass, ordMass, ordNRQCD,
   ord1S, R1S, "numeric", lambda, gt, mZ, gammaZ, sin2ThetaW, alphamZ, alphaQED,
   mT, eH, hnu, v1, v2, energy_cm, 0.018, 0.022, angle, half_aperture) << "\n";

    return 0;
}
