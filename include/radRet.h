#ifndef RADRET_H
#define RADRET_H

double sigmahad(char const* str, char const* curr, int orderAlp, int runAlp,
int order, int nf, double mZ, double gammaZ, double sin2ThetaW, double amZ,
double amZQED, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double Q);

double sigmarad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x, double theta);

double sigmaradcum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x0, double x1, double theta);

double sigmaradcone(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x, double theta,
double deltaTheta);

double sigmaradconecum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x0, double x1, double theta,
double deltaTheta);

double sigmamass(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double mu, double Q);

double sigmamassrad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x, double theta);

double sigmamassradcum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x0, double x1,
double theta);

double sigmamassradcone(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x, double theta,
double deltaTheta);

double sigmamassradconecum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x0, double x1,
double theta, double deltaTheta);

double sigmamatched(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double gammaZ, double sinW, double aMz, double aMzQED, double mT,
double mu, double nu, double v1, double v2, double Q);

double sigmamatchedrad(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x,
double theta);

double sigmamatchedradcum(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x0,
double x1, double theta);

double sigmamatchedradcone(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x,
double theta, double deltatheta);

double sigmamatchedradconecum(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x0,
double x1, double theta, double deltatheta);

#endif // RADRET_H
