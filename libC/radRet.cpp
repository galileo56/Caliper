extern "C" void f90sigmahad_(char const* str, char const* curr, int* orderAlp,
int* runAlp, int* order, int* nf, double* mZ, double* gammaZ, double* sin2ThetaW,
double* amZ, double* amZQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double* Q, double* res);

extern "C" void f90sigmarad_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x, double* theta,
double* res);

extern "C" void f90sigmaradcum_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x0, double* x1, double* theta,
double* res);

extern "C" void f90sigmaradcone_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x, double* theta,
double* deltaTheta, double* res);

extern "C" void f90sigmaradconecum_(char const* str, char const* curr,
int* orderAlpha, int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ,
double* thetaW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* theta,
double* x0, double* x1, double* deltaTheta, double* res);

extern "C" void f90sigmamass_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* mu, double * Q, double* res);

extern "C" void f90sigmamassrad_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x,
double* theta, double* res);

extern "C" void f90sigmamassradcum_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x0,
double* x1, double* theta, double* res);

extern "C" void f90sigmamassradcone_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x,
double* theta, double* deltaTheta, double* res);

extern "C" void f90sigmamassradconecum_(char const* str, char const* curr,
int* orderAlpha, int* runAlpha, int* runMass, int* order, int* nf, double* Mz,
double* gammaZ, double* sinW, double* aMz, double* aMzQED, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* eH,
double * Q, double* x0, double* x1, double* theta, double* deltaTheta,
double* res);

extern "C" void f90sigmamatched_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* res);

extern "C" void f90sigmamatchedrad_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x, double* theta, double* res);

extern "C" void f90sigmamatchedradcum_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x0, double* x1, double* theta, double* res);

extern "C" f90sigmamatchedradcone_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x, double* theta, double* deltatheta, double* res);

extern "C"f90sigmamatchedradconecum_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x0, double* x1, double* theta, double* deltatheta,
double* res);

double sigmahad(char const* str, char const* curr, int orderAlp, int runAlp,
int order, int nf, double mZ, double gammaZ, double sin2ThetaW, double amZ,
double amZQED, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double Q){

  double res;

  f90sigmahad_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
  &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

  return res;
}

double sigmarad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x, double theta){

  double res;

  f90sigmarad_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
  &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
  &Q, &x, &theta, &res);

  return res;
}

double sigmaradcum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x0, double x1, double theta){

   double res;

   f90sigmaradcum_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
   &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
   &Q, &x0, &x1, &theta, &res);

  return res;
}

double sigmaradcone(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x, double theta,
double deltaTheta){

   double res;

   f90sigmaradcone_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
   &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
   &Q, &x, &theta, &deltaTheta, &res);

  return res;
}

double sigmaradconecum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x0, double x1, double theta,
double deltaTheta){

   double res;

   f90sigmaradconecum_(str, curr, &orderAlpha, &runAlpha, &order, &nf,
   &Mz, &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
   &Q, &x0, &x1, &theta, &deltaTheta, &res);

  return res;
}

double sigmamass(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double mu, double Q){

  double res;

  f90sigmamass_(str, curr, &orderAlpha, &runAlpha, &runMass, &order, &nf, &Mz,
  &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

 return res;
}

static double sigmamassrad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x, double theta){

  double res;

  f90sigmamassrad_(str, curr, &orderAlpha, &runAlpha, &runMass, &order,
  &nf, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
  &Q, &x, &theta, &res);

 return res;
}

double sigmamassradcum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x0, double x1,
double theta){

  double res;

  f90sigmamassradcum_(str, curr, &orderAlpha, &runAlpha, &runMass, &order,
  &nf, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
  &Q, &x0, &x1, &theta, &res);

 return res;
}

double sigmamassradcone(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x, double theta,
double deltaTheta){

  double res;

  f90sigmamassradcone_(str, curr, &orderAlpha, &runAlpha, &runMass, &order,
  &nf, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
  &Q, &x, &theta, &deltaTheta, &res);

 return res;
}

double sigmamassradconecum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double eH, double Q, double x0, double x1,
double theta, double deltaTheta){

  double res;

  f90sigmamassradconecum_(str, curr, &orderAlpha, &runAlpha, &runMass, &order,
  &nf, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
  &Q, &x0, &x1, &theta, &deltaTheta, &res);

 return res;
}

double sigmamatched(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double gammaZ, double sinW, double aMz, double aMzQED, double mT,
double mu, double nu, double v1, double v2, double Q){

  double res;

  f90sigmamatched_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S,
  method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &res);

 return res;
}

double sigmamatchedrad(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x,
double theta){

  double res;

  f90sigmamatchedrad_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S,
  &R1S, method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &x, &theta, &res);

 return res;
}

double sigmamatchedradcum(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x0,
double x1, double theta){

  double res;

  f90sigmamatchedradcum_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S,
  &R1S, method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &x0, &x1, &theta, &res);

 return res;
}

double sigmamatchedradcone(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x,
double theta, double deltatheta){

  double res;

  f90sigmamatchedradcone_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S,
  &R1S, method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &x, &theta, &deltatheta, &res);

 return res;
}

double sigmamatchedradconecum(char const* str, int runAlpha, int runMass,
int ordMass, int order, int ord1S, double R1S, char const* method, double lambda,
double gt, double Mz, double gammaZ, double sinW, double aMz, double aMzQED,
double mT, double mu, double nu, double v1, double v2, double Q, double x0,
double x1, double theta, double deltatheta){

  double res;

  f90sigmamatchedradconecum_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S,
  &R1S, method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &x0, &x1, &theta, &deltatheta, &res);

 return res;
}
