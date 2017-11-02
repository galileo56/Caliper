:Evaluate:   BeginPackage["Caliper`"]

:Evaluate:   Print["     Package for Massive and Massless Event Shapes "]
:Evaluate:   Print["     Author:            Vicent Mateu               "]
:Evaluate:   Print["     Last modification: 12 - 05 - 2017             "]
:Evaluate:   Print["     Version:           test 1                     "]

:Evaluate:  amZdef           = 0.1181
:Evaluate:  amZdelt          = 0.0011
:Evaluate:  mCdef            = 1.28
:Evaluate:  mCdelt           = 0.03
:Evaluate:  mZdef            = 91.187
:Evaluate:  Gammadef         = 1553.0647546066
:Evaluate:  gammaZdef        = 2.4952
:Evaluate:  sin2ThetaWdef    = 0.23119
:Evaluate:  gammaZPythia     = 2.5042
:Evaluate:  sin2ThetaWPythia = 0.2312
:Evaluate:  aQEDdef          = 0.00781751
:Evaluate:  MUpsilonexp      = {9.399, 9.46030, 9.85944, 9.89278, 9.8993, 9.91221, 9.999, 10.02326, 10.1637, 10.2325, 10.25546, 10.2598, 10.26865, 10.3552}
:Evaluate:  ErrUpsilon       = {0.0023, 0.00026, 0.00042, 0.00031, 0.0008, 0.00031, 0.004, 0.00031, 0.0014, 0.0005, 0.0005, 0.0012, 0.0005, 0.0005}
:Evaluate:  MJexp            = {2.98034, 3.0969}
:Evaluate:  MJexpErr         = {0.0005, 0.000006}
:Evaluate:  MJexpErr         = {0.0005, 0.000006}
:Evaluate:  MJexpErrdataJ    = {MJexp, MJexpErr} // Transpose
:Evaluate:  Narr             = {1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};
:Evaluate:  Larr             = {0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 0};
:Evaluate:  Sarr             = {0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1};
:Evaluate:  Jarr             = {0, 1, 0, 1, 1, 2, 0, 1, 2, 0, 1, 1, 2, 1};
:Evaluate:  Quantum          = {Narr, Larr, Sarr, Jarr} // Transpose

:Evaluate:  QSwitch::usage = "Delta1S[nl, orderAlpha, runAlpha, orderMass, runMass, ord1S, muLam, xLam, method, mZ, aMz, mt, gt, R]"
:Evaluate:  Delta1S::usage = "Delta1S[nl, orderAlpha, runAlpha, orderMass, runMass, muLam, xLam, method, mZ, aMz, mt, R]"
:Evaluate:  rNRQCD::usage = "rNRQCD[nl, order, scheme, method, orderAlpha, runAlpha, orderMass, runMass, ord1S, R1S, muLam, xLam, mZ, aMz, Q, mtpole, gt, h, nu]"
:Evaluate:  A1Pole::usage = "A1Pole[nl, order, En, mtpole, gamtop, asoft, VcsNNLL, musoft]"
:Evaluate:  TTbar::usage = "ttbar[energy, topmass, topgamma, alphas0, mue0, cutn, cutv, c0, c1, c2, cdeltapotc, cdeltapot1, cfullc, cfull1, crm2, kincm, kinca, ijknflg, ijgcflg, kincv, ijvflg] cross section"
:Evaluate:  TTbarList::usage = "ttbarList[energy, topmass, topgamma, alphas0, mue0, cutn, cutv, c0, c1, c2, cdeltapotc, cdeltapot1, cfullc, cfull1, crm2, kincm, kinca, ijknflg, ijgcflg, kincv, ijvflg] cross section and distribution list"
:Evaluate:  CdiGamma::usage = "CdiGamma[x]"
:Evaluate:  CtriGamma::usage = "CtriGamma[x]"
:Evaluate:  HypGeo::usage = "HypGeo[a,b,c,z]"
:Evaluate:  QFromV::usage = "QFromV[v, m, gt]"
:Evaluate:  VC::usage = "VC[q, m, gt]"
:Evaluate:  VStar::usage = "VStar[q, m, gt]"
:Evaluate:  VRootStar::usage = "VRootStar[q, m, gt]"
:Evaluate:  SwitchOff::usage = "SwitchOff[q, m, gt, v0, v1]"
:Evaluate:  VssLL::usage = "VssLL[nl, ah, as]"
:Evaluate:  Vk1sLL::usage = "Vk1sLL[nl, as, au]"
:Evaluate:  Vk2sLL::usage = "Vk2sLL[nl, ah, as]"
:Evaluate:  VkeffsLL::usage = "VkeffsLL[nl, ah, as]"
:Evaluate:  VcsLL::usage = "VcsLL[as]"
:Evaluate:  VrsLL::usage = "VrsLL[nl, as, au]"
:Evaluate:  V2sLL::usage = "V2sLL[nl, ah, as, au]"
:Evaluate:  XiNLL::usage = "XiNLL[nl, ah, as, au]"
:Evaluate:  VceffsNNLL::usage = "VceffsNNLL[nl, asNNLL, as, au]"
:Evaluate:  XiNNLLmixUsoft::usage = "XiNNLLmixUsoft[nl, ah, as]"
:Evaluate:  MLLc2::usage = "MLLc2[nl, ah, au]"
:Evaluate:  MNLLc1::usage = "MNLLc1[nl, ah, as, au]"
:Evaluate:  MNLLplusNNLLnonmixc1::usage = "MNLLplusNNLLnonmixc1[nl, ah, as, au]"
:Evaluate:  MNNLLAllc1InclSoftMixLog::usage = "MNNLLAllc1InclSoftMixLog[nl, ah, as, au, nu, hh, ss]"
:Evaluate:  XiNNLLSoftMixLogc1::usage = "XiNNLLSoftMixLogc1[ah, nu, hh]"
:Evaluate:  XiNNLLnonmix::usage = "XiNNLLnonmix[nl, ah, as, au, hh, ss]"
:Evaluate:  DeltaBottomCharm::usage = "DeltaBottomCharm[r1,r2] double massive bubble"
:Evaluate:  GammaRBottomCharm::usage = "GammaRBottomCharm[r1,r2] R-anomalous dimension from the double massive bubble"
:Evaluate:  Pi0::usage = "Pi0[z] tree-level massive vacuum polarization function"
:Evaluate:  Pi0Der::usage = "Pi0Der[i,z] i-th derivative of the tree-level massive vacuum polarization function"
:Evaluate:  Pi1Der::usage = "Pi1Der[i,z] i-th derivative of the one-loop massive vacuum polarization function"
:Evaluate:  Pi1::usage = "Pi1[z] one-loop massive vacuum polarization function"
:Evaluate:  Pi3::usage = "Pi3[z] three-loop massive vacuum polarization function"
:Evaluate:  Pi2::usage = "Pi2[z] two-loop massive vacuum polarization function"
:Evaluate:  Pi2Der::usage = "Pi2Der[z] derivative of the two-loop massive vacuum polarization function"
:Evaluate:  P2::usage = "P2[z] integrand for massive bubble"
:Evaluate:  P2Double::usage = "P2Double[r1,r2] integral for double massive bubble"
:Evaluate:  P2Int::usage = "P2Int[r] integral for massive bubble"
:Evaluate:  DeltaCharmExact::usage = "DeltaCharmExact[charm, type, scheme, average, n, l, j, s, nl, mH, mL, mu, alp] computes the subleading massive charm corrections to quarkonium masses"
:Evaluate:  UpsilonDeltaCharmBin::usage = "UpsilonDeltaCharmBin[n, l, alp, mb, mc] computes the massive charm corrections to quarkonium masses"
:Evaluate:  UpsilonDeltaCharm::usage = "UpsilonDeltaCharm[n, l, alp, mb, mc] computes the massive charm corrections to quarkonium masses"
:Evaluate:  GammaRCharm2::usage = "GammaRCharm2[z] computes the massive charm corrections to the R-anomalous dimension at two loops"
:Evaluate:  GammaRCharm3::usage = "GammaRCharm3[nl, nh, z] computes the massive charm corrections to the R-anomalous dimension at three loops"
:Evaluate:  DeltaCharmGlue::usage = "DeltaCharmGlue[z] computes the massive charm corrections to the nm and nm^2 piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharmGlueDer::usage = "DeltaCharmGlueDer[z] computes the derivative of the massive charm corrections to the nm and nm^2 piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharmNl::usage = "DeltaCharmNl[z] computes the massive charm corrections to the nl piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharmNh::usage = "DeltaCharmNh[z] computes the massive charm corrections to the nh piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharmNhDer::usage = "DeltaCharmNhDer[z] computes the derivative of the massive charm corrections to the nh piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharmNlDer::usage = "DeltaCharmNlDer[z] computes the derivative of the massive charm corrections to the nl piece of the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharm2::usage = "DeltaCharm2[z] computes the massive charm corrections to the MSbar-pole mass relation at two loops"
:Evaluate:  DeltaCharm3::usage = "DeltaCharm3[nl, nh, z] computes the massive charm corrections to the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharm3Der::usage = "DeltaCharm3Der[nl, nh, z] computes the derivative of the massive charm corrections to the MSbar-pole mass relation at three loops"
:Evaluate:  DeltaCharm2Der::usage = "DeltaCharm2Der[z] computes the derivative of the massive charm corrections to the MSbar-pole mass relation at two loops"
:Evaluate:  LegendreList::usage = "LegendreList[n, k, x] computes the of the first n + 1 k-th derivative of the Legendre Polynomials"
:Evaluate:  QLegendreList::usage = "QLegendreList[n, x] computes the of the first n + 1 Legendre Polynomial"
:Evaluate:  BreitUnstable::usage = "BreitUnstable[shape, mt, Q, gamma, n, k, x] computes the LO distribution for unstable tops convoluted with a BreitWigner"
:Evaluate:  MCtop::usage = "MCtop[shape, mt, Q, n, k, x] computes the LO distribution for unstable tops"
:Evaluate:  DeltaMCtop::usage = "DeltaMCtop[shape, mt, Q] computes the LO distribution for unstable tops"
:Evaluate:  pFq::usage = "pFq[a,b,z] computes the generalized hypergeometric function"
:Evaluate:  FindOrigin::usage = "FindOrigin[shape, gap, orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eR, R0, muR0, delta0, h], finds the origin for massless Event Shapes"
:Evaluate:  MassOrigin::usage = "MassOrigin[shape, EShape, gap, scheme, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, R0, muR0, del0, h], finds the origin of Massive Event Shapes"
:Evaluate:  MasslessPieceBin::usage = "MasslessPieceBin[terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda, R0, muR0, delta0, h, tauList], computes bins including profiles for massless cross section"
:Evaluate:  MasslessBinList::usage = "MasslessBinList[terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0, muR0, delta0, h, tauList], computes bins including profiles for massless cross section"
:Evaluate:  MasslessProfList::usage = "MasslessProfList[terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0, muR0, delta0, h, tauList], computes the cross section including profiles for massless cross section"
:Evaluate:  MassiveBinList::usage = "MassiveBinList[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList], computes bins including profiles for massive cross section"
:Evaluate:  MassiveProfList::usage = "MassiveProfList[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList], computes the cross section including profiles for massive cross section"
:Evaluate:  MassiveProf::usage = "MassiveProf[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau], computes the cross section including profiles for massive cross section"
:Evaluate:  MassiveProfPiece::usage = "MassiveProfPiece[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau], computes the cross section including profiles for massive cross section"
:Evaluate:  MassiveProfPieceList::usage = "MassiveProfPieceList[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau], computes the cross section including profiles for massive cross section"
:Evaluate:  MassivePieceBin::usage = "MassivePieceBin[terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList], computes the cross section including profiles for massive cross section"
:Evaluate:  MassiveMoment::usage = "MassiveMoment[terms, hard, shape, EShape, setup, gap, space, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, pow], computes the moment including profiles for massive cross section"
:Evaluate:  MasslessMoment::usage = "MasslessMoment[terms, hard, shape, setup, gap, space, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0, muR0, delta0, h, tau, tau2, pow], computes moments including profiles for massless cross section"
:Evaluate:  MasslessProf::usage = "MasslessProf[terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0, muR0, delta0, h, tau], computes the cross section including profiles for massless cross section"
:Evaluate:  MasslessProfPiece::usage = "MasslessProfPiece[terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda, R0, muR0, delta0, h, tau], computes the cross section including profiles for massless cross section"
:Evaluate:  MasslessProfPieceList::usage = "MasslessProfPieceList[terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda, R0, muR0, delta0, h, tau], computes the cross section including profiles for massless cross section"
:Evaluate:  MassIntersection::usage = "MassIntersection[Q, beta, mu0, delLamb, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, def, EShape], computes the Renormalization scales for massive event shapes"
:Evaluate:  ProfilesMass::usage = "ProfilesMass[Q, beta, mu0, delLamb, R0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape, tau], computes the Renormalization scales for massive event shapes"
:Evaluate:  Profiles::usage = "Profiles[Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, tau], computes the Renormalization scales for massless event shapes"
:Evaluate:  GammaR::usage = "GammaR[str, nf] computes the soft R anomalous Dimension"
:Evaluate:  DiLog::usage = "DiLog[z] computes the polylogarithm"
:Evaluate:  Elliptic3::usage = "Elliptic3[psi, k, c] computes the elliptic function of the third kind, using Carlson forms"
:Evaluate:  Polylog::usage = "Polylog[n,z] computes the polylogarithm"
:Evaluate:  NGLFunction::usage = "NGLFunction[n,z] computes the non-global part of the soft function"
:Evaluate:  ComplexPolylog::usage = "ComplexPolylog[n,z] computes the polylogarithm"
:Evaluate:  NGLSoft::usage = "NGLSoft[nf,z] computes the NGL function in complex space"
:Evaluate:  CLi2::usage = "Cli2[z] computes the complex dilogarithm"
:Evaluate:  CLi3::usage = "Cli3[z] computes the complex trilogarithm"
:Evaluate:  sCoef::usage = "sCoef[str, nf] computes the soft R anomalous Dimension"
:Evaluate:  sCoefGamma::usage = "sCoefGamma[gamma, n, nf] computes the soft R anomalous Dimension"
:Evaluate:  sCoefLambda::usage = "sCoefLambda[str, nf, lambda] computes the MSR R-anomalous Dimension"
:Evaluate:  Ql::usage = "Ql[str, nf, order, lambda] computes the Ql coefficients"
:Evaluate:  aNasy::usage = "aNasy[str, nf, order, lambda] computes the asymptotic aN coefficients"
:Evaluate:  aFromS::usage = "aFromS[str, nf, lambda] computes the a coefficients from the sCoefs"
:Evaluate:  anLambda::usage = "anLambda[str, nf, lambda] computes the lambda dependence of the a coefficients"
:Evaluate:  AnomDim::usage = "AnomDim[str, nf, G4] computes the QCD anomalous dimension"
:Evaluate:  MSbarDeltaPiece::usage = "MSbarDeltaPiece[nl, nh] computes the pole to MS-bar relation"
:Evaluate:  AlphaMatchingLog::usage = "AlphaMatchingLog[str, direction, nf] computes the alpha threshold matching"
:Evaluate:  AlphaMatchingInverse::usage = "AlphaMatchingInverse[str, nf] computes the inverse alpha threshold matching"
:Evaluate:  cCoef::usage = "cCoef[nf, order, m] computes the inverse of the QCD anomalous dimension"
:Evaluate:  PSCoef::usage = "PSCoef[nf, lg] computes the PS mass series coefficients"
:Evaluate:  N12Generic::usage = "N12Generic[aCoef, order, nf, lambda] computes the renormalon sum rule"
:Evaluate:  N12::usage = "N12[str, order, nf, lambda, err] computes the renormalon sum rule"
:Evaluate:  N12Ratio::usage = "N12Ratio[str, order, nf, lambda, err] computes the renormalon residue by the ratio method"
:Evaluate:  N12Residue::usage = "N12Residue[str, order, nf, lambda, err] computes the renormalon residue by the residue method"
:Evaluate:  N12RS::usage = "N12RS[str, order, n, nf, lambda, err] computes the renormalon residue by the MSR-RS mass consistency"
:Evaluate:  Delta::usage = "Delta[str, nf, mu, R] computes the soft renormalon subtractions"
:Evaluate:  DeltaGap::usage = "DeltaGap[str, orderAlpha, runAlpha, runMass, nf, mZ, aMz, mT, muT, mB, muB, mC, muC, mu, R] computes the soft renormalon subtractions"
:Evaluate:  PSDelta::usage = "PSDelta[orderAlpha, runAlpha, nf, mZ, aMz, mT, muT, mB, muB, mC, muC, mu, R, lg] computes the PS mass subtractions"
:Evaluate:  CoefMat::usage = "CoefMat[str, nf, s3] computes the hard, soft and jet matrix elements"
:Evaluate:  wTilde::usage = "wTilde[order, nf, gamma, a0, a1] computes wTilde for a given anomalous dimension gamma"
:Evaluate:  kTilde::usage = "kTilde[order, nf, gamma, a0, a1] computes kTilde for a given anomalous dimension gamma"
:Evaluate:  AlphaQED::usage = "AlphaQED[nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the electromagnetic coupling with flavor matching."
:Evaluate:  AlphaQCD::usage = "AlphaQCD[scheme, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of alpha with flavor matching."
:Evaluate:  AlphaComplex::usage = "AlphaComplex[scheme, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of alpha with flavor matching."
:Evaluate:  MSbarMass::usage = "MSbarMass[order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the quark masses with flavor matching."
:Evaluate:  PoleMass::usage = "PoleMass[orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the quark masses with flavor matching."
:Evaluate:  MSbarMassLow::usage = "MSbarMassLow[order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the quark masses with flavor matching below the mass."
:Evaluate:  RSMass::usage = "RSMass[type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, lambda, R] computes the RS-scheme running of the quark masses."
:Evaluate:  MSRMass::usage = "MSRMass[type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, lambda, R] computes the MSR running of the quark masses."
:Evaluate:  MSRVFNS::usage = "MSRVFNS[up, type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, R] computes the MSR running of the quark masses with flavor matching."
:Evaluate:  MSRTop::usage = "MSRTop[up, type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, mu3, R] computes the MSR running of the top masses for nonzero bottom and charm quark masses."
:Evaluate:  NRQCD::usage = "NRQCD[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] computes the quarkonium energy levels."
:Evaluate:  NRQCDDerCharm::usage = "NRQCDDerCharm[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R, eps] computes the derivative of quarkonium energy levels wrt the charm mass."
:Evaluate:  NRQCDDerAlpha::usage = "NRQCDDerAlpha[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R, eps] computes the derivative of quarkonium energy levels wrt the amZ."
:Evaluate:  MassIter::usage = "MassIter[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu, R] computes the bottom mass from quarkonium energy levels."
:Evaluate:  MassExpand::usage = "MassExpand[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu, R] computes the bottom mass from quarkonium energy levels."
:Evaluate:  FindMass::usage = "FindMass[ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu, R] fits the quark mass from the quarkonium energy levels."
:Evaluate:  FindEnergy::usage = "FindEnergy[ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] finds the quarkonium energy levels."
:Evaluate:  MassError::usage = "MassError[ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x] fits the quark mass from the quarkonium energy levels, including perturbative error."
:Evaluate:  MassList::usage = "MassList[ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR] makes a list of the quark mass from the quarkonium energy levels in a grid of mu-R values."
:Evaluate:  NRQCDList::usage = "NRQCDList[n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR] makes a list of the NRQCD prediction for the quarkonium energy levels in a grid of mu-R values."
:Evaluate:  UpsilonList::usage = "UpsilonList[n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm] makes a list of the NRQCD prediction for the quarkonium energy levels and their derivatives wrt alpha(mZ) and mC, in a grid of mu-R values."
:Evaluate:  CorrMat::usage = "CorrMat[qnlist, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm] Computes the average values of the masses and derivatives wrt alpha and mc, perturbative uncertainties and covariance matrix."
:Evaluate:  Chi2NRQCD::usage = "Chi2NRQCD[qnlist, datalist, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] Computes the chi2 value with respect of the heavy quark mass for quarkonium."
:Evaluate:  Chi2MinNRQCD::usage = "Chi2MinNRQCD[qnlist, datalist, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] minimizes the chi2 value with respect of the heavy quark mass and alphaS(mZ) for quarkonium."
:Evaluate:  Chi2MinAlphaNRQCD::usage = "Chi2MinAlphaNRQCD[qnlist, datalist, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] minimizes the chi2 value with respect to alphaS(mZ) for quarkonium."
:Evaluate:  Chi2MinAlphaMbNRQCD::usage = "Chi2MinAlphaMbNRQCD[qnlist, datalist, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu, R] minimizes the chi2 value for quarkonium."
:Evaluate:  ErrMat::usage = "ErrMat[qnlist, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm] Computes the average values of the masses and derivatives wrt alpha and mc, perturbative uncertainties and covariance matrix."
:Evaluate:  ErrMatrices::usage = "ErrMatrices[qnlist, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm] Computes the average values of the masses and derivatives wrt alpha and mc, perturbative uncertainties and covariance matrix."
:Evaluate:  NRQCDError::usage = "NRQCDError[n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x] computes the quarkonium energy levels, including perturbative error."
:Evaluate:  OptimalRVFNS::usage = "OptimalRVFNS[up, type, n, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2] computes the Optimal R scale for quarkonium."
:Evaluate:  OptimalR::usage = "OptimalR[type, n, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, lambda] computes the Optimal R scale for quarkonium."
:Evaluate:  OptimalR2::usage = "OptimalR2[n, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mass] computes the Optimal R scale for quarkonium."
:Evaluate:  mmfromMSR::usage = "mmfromMSR[type, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R] computes the MSR practical definition running of the quark masses with flavor matching."
:Evaluate:  Rhad::usage = "Rhad[scheme, orderAlpha, runAlpha, order, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu, Q] computes the massless total hadronic cross section."
:Evaluate:  SigmaHad::usage = "SigmaHad[scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, mu, Q] computes the massless total hadronic cross section."
:Evaluate:  SigmaRad::usage = "SigmaRad[scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x, theta] computes the ISR massless total hadronic cross section."
:Evaluate:  SigmaRadCum::usage = "SigmaRadCum[scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, theta] computes the ISR massless total hadronic cross section."
:Evaluate:  SigmaRadCone::usage = "SigmaRadCone[scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x, theta, deltaTheta] computes the ISR massless total hadronic cross section."
:Evaluate:  SigmaRadConeCum::usage = "SigmaRadConeCum[scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, theta, deltaTheta] computes the ISR massless total hadronic cross section."
:Evaluate:  RhadCoefs::usage = "RhadCoefs[nf] computes the massless total hadronic cross section series coefficients."
:Evaluate:  RhadMass::usage = "RhadMass[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, mu, Q] computes the massive total hadronic cross section."
:Evaluate:  SigmaMass::usage = "SigmaMass[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, mu, Q] computes the massive total hadronic cross section."
:Evaluate:  SigmaMassRad::usage = "SigmaMassRad[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, eH, Q, x, theta] computes the ISR massive total hadronic cross section."
:Evaluate:  SigmaMassRadCum::usage = "SigmaMassRadCum[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, theta] computes the ISR massive total hadronic cross section."
:Evaluate:  SigmaMassRadCone::usage = "SigmaMassRadCone[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, eH, Q, x, theta, deltaTheta] computes the ISR massive total hadronic cross section."
:Evaluate:  SigmaMassRadConeCum::usage = "SigmaMassRadConeCum[scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, theta, deltaTheta] computes the ISR massive total hadronic cross section."
:Evaluate:  RQCD::usage = "RQCD[scheme, runAlpha, runMass, ordMass, ord1S, R1S, order, method, lambda, gt, Mz, aMz, mT, mu, Q] computes the massive total hadronic cross section for an unstable top quark."
:Evaluate:  RExp::usage = "RExp[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, aMz, mT, mu, nu, Q] computes the threshold-expaded massive total hadronic cross section for an unstable top quark."
:Evaluate:  SigmaMatched::usage = "Rmatched[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, hnu, v1, v2, Q] computes the matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  SigmaMatchedRadCum::usage = "SigmaMatchedRadCum[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, hnu, v1, v2, Q, x0, x1, theta] computes the ISR matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  SigmaMatchedRad::usage = "SigmaMatchedRad[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, hnu, v1, v2, Q, x, theta] computes the ISR matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  SigmaMatchedRadCone::usage = "SigmaMatchedRadCone[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, hnu, v1, v2, Q, x, theta, deltaTheta] computes the ISR matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  SigmaMatchedRadConeCum::usage = "SigmaMatchedRadConeCum[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, hnu, v1, v2, Q, x0, x1, theta, deltaTheta] computes the ISR matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  Rmatched::usage = "Rmatched[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, aMz, mT, mu, nu, v1, v2, Q] computes the matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  RmatchedList::usage = "RmatchedList[scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method, lambda, gt, Mz, aMz, mT, h, hnu, v1, v2, Q0, Q1, deltaQ] computes the matched massive total hadronic cross section for an unstable top quark."
:Evaluate:  LambdaQCD::usage = "LambdaQCD[scheme, order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the quark masses with flavor matching."
:Evaluate:  Hyper2F1::usage="Hyper2F1[a, b, c, x] Hypergeometric Function in Fortran"
:Evaluate:  HyperF32Exact::usage="HyperF32Exact[w, x] Hypergeometric Function in Fortran"
:Evaluate:  InteCorre::usage="InteCorre[b, x0, x1] Incomplete Gamma Function"
:Evaluate:  DiffDeltaGap::usage = "DiffDeltaGap[gap, scheme, order, R0, R1, mu0, mu1, muLambda, orderAlpha, runAlpha, nf, Mz, aMz, mT, muT, mB, muB, mC, muC] computes the running of the gap parameter."
:Evaluate:  DiffDeltaGapMass::usage = "DiffDeltaGapMass[gap, order, R0, R1, mu0, mu1, muM, muLambda1, muLambda2, orderAlpha, runAlpha, runMass, nf, Mz, aMz, mT, muT, mB, muB, mC, muC] computes the running of the gap parameter with flavor matching."
:Evaluate:  Kernel::usage = "Kernel[n, width, w, mu, p] computes the first n+1 kernels needed for resummation"
:Evaluate:  GammaDerList::usage = "GammaDerList[n, w] computes the first n+1 derivatives of 1/Gamma"
:Evaluate:  polyGamma::usage = "polyGamma[n, w] computes the first n+1 derivatives of Gamma"
:Evaluate:  NGLKernel::usage = "NGLKernel[n, n1, n2, width, w, mu, p] computes the first 2*n kernels needed for NGL resummation"
:Evaluate:  Taylor::usage = "Taylor[c, lambda, k] computes the Taylor expansion of the shape function"
:Evaluate:  Model::usage = "Model[c, lambda, k, l] computes the shape function"
:Evaluate:  ModelUnstable::usage = "ModelUnstable[shape, mt, Q, c, lambda, n, k, l] computes the shape function convoluted with the unstable distribution"
:Evaluate:  BreitModelUnstable::usage = "BreitModelUnstable[shape, mt, Q, gamma, c, lambda, n, k, l] computes the shape function convoluted with the unstable distribution plus a BreitWigner"
:Evaluate:  BreitModel::usage = "BreitModel[c, lambda, width, k, l] computes the shape function convoluted with a Breit Wigner"
:Evaluate:  MomentModel::usage = "MomentModel[c, lambda, k] computes the shape function"
:Evaluate:  ModelPiece::usage = "ModelPiece[c, lambda, k, l] computes the shape function"
:Evaluate:  TaylorPiece::usage = "TaylorPiece[c, lambda, k] computes the Taylor coefficients of the shape function"
:Evaluate:  DeltaMSbar::usage = "DeltaMSbar[order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of the quark masses with flavor matching."
:Evaluate:  JetMass::usage = "JetMass[orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R, mu] computes the Jet Mass running of the quark masses with flavor matching."
:Evaluate:  mmFromJetMass::usage = "mmFromJetMass[orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R, mu] computes the Jet Mass running of the quark masses with flavor matching."
:Evaluate:  Singular::usage = "Singular[hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularList::usage = "SingularList[hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, clen, lambda, R0, mu0, delta0, h, tau] computes the pieces of the Singular Thrust and C-parameter distrubution"
:Evaluate:  MassNonDist::usage = "MassNonDist[hard, shape, Eshape, setup, gap, space, cum, scheme, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Non distributional part of the Singular Massive Thrust and C-parameter distrubution"
:Evaluate:  MassNonDistPiece::usage = "MassNonDistPiece[hard, shape, Eshape, gap, space, cum, scheme, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Non distributional part of the Singular Massive Thrust and C-parameter distrubution"
:Evaluate:  SingularMass::usage = "SingularMass[hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau] computes the Singular Massive Thrust and C-parameter distrubution"
:Evaluate:  SingularMassPiece::usage = "SingularMassPiece[hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau] computes the Singular Massive Thrust and C-parameter distrubution"
:Evaluate:  SingularPiece::usage = "SingularPiece[hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularHJM::usage = "SingularHJM[hard, setup, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularHJMPiece::usage = "SingularHJMPiece[hard, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularDouble::usage = "SingularDouble[hard, setup, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularDoublePiece::usage = "SingularDoublePiece[hard, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2] computes the Singular Thrust and C-parameter distrubution"
:Evaluate:  SingularHJM1D::usage = "SingularHJM1D[hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular HJM distrubution with a 1D shape function"
:Evaluate:  SingularHJM1DPiece::usage = "SingularHJM1DPiece[hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singular HJM distrubution with the piece of a 1D shape function"
:Evaluate:  EWFactors::usage = "EWFactors[nf, Q, Mz, GammaZ, sin2ThetaW] electroweak factors"
:Evaluate:  NGLIntegral::usage = "NGLIntegral[nf, pow, w1, w2] computes the exact NGL integral"
:Evaluate:  NGLDoubleIntegral::usage = "NGLDoubleIntegral[nf, pow, w1, w2, r] computes the exact NGL integral"
:Evaluate:  ThrustNS1loop::usage = "ThrustNS1loop[tau] computes the 1-loop thrust NS function and its first two derivatives"
:Evaluate:  FOMass::usage = "FOMass[shape, current, m, Q, Mz, gammaZ, sin2ThetaW, tau] computes the 1-loop massive thrust NS function"
:Evaluate:  ThrustNS2loop::usage = "ThrustNS2loop[er, tau] computes the 2-loop thrust NS function and its first two derivatives"
:Evaluate:  NSMass::usage = "NSMass[shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, Mz, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive Thrust, HJM, SJM and C-parameter distribution, with a 1D model."
:Evaluate:  NSMassPiece::usage = "NSMassPiece[shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, Mz, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive Thrust, HJM, SJM, and C-parameter distribution, with a 1D model."
:Evaluate:  HJMNSMass::usage = "HJMNSMass[setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, Mz, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive HJM distribution, with a 2-D model."

:Evaluate:  Begin["`Private`"]

:Evaluate:  Print["You can access the complete function list typing '?Caliper`*' "]

:Begin:
:Function:      hypgeo
:Pattern:       HypGeo[a_, b_, c_, z_]
:Arguments:     {Re[a], Im[a], Re[b], Im[b], Re[c], Im[c], Re[z], Im[z]}
:ArgumentTypes: {Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      cdigamma
:Pattern:       CdiGamma[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ctrigamma
:Pattern:       CtriGamma[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      xinnllnonmix
:Pattern:       XiNNLLnonmix[nl_, ah_, as_, au_, hh_, ss_]
:Arguments:     {nl, ah, as, au, hh, ss}
:ArgumentTypes: {Integer, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      xinnllsoftmixlogc1
:Pattern:       XiNNLLSoftMixLogc1[ah_, nu_, hh_]
:Arguments:     {ah, nu, hh}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      mnnllallc1inclsoftmixlog
:Pattern:       MNNLLAllc1InclSoftMixLog[nl_, ah_, as_, au_, nu_, hh_, ss_]
:Arguments:     {nl, ah, as, au, nu, hh, ss}
:ArgumentTypes: {Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      mnllplusnnllnonmixc1
:Pattern:       MNLLplusNNLLnonmixc1[nl_, ah_, as_, au_]
:Arguments:     {nl, ah, as, au}
:ArgumentTypes: {Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      mnllc1
:Pattern:       MNLLc1[nl_, ah_, as_, au_]
:Arguments:     {nl, ah, as, au}
:ArgumentTypes: {Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vceffsnnll
:Pattern:       VceffsNNLL[nl_, asNNLL_, ah_, as_]
:Arguments:     {nl, asNNLL, ah, as}
:ArgumentTypes: {Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rnrqcd
:Pattern:       rNRQCD[nl_, order_, scheme_, method_, orderAlpha_, runAlpha_,
                orderMass_, runMass_, ord1S_, R1S_, muLam_, xLam_, mZ_, aMz_,
                Q_, mtpole_, gt_, h_, nu_]
:Arguments:     {nl, order, scheme, method, orderAlpha, runAlpha, orderMass,
                runMass, ord1S, R1S, muLam, xLam, mZ, aMz, Q, mtpole, gt, h, nu}
:ArgumentTypes: {Integer, Integer, String, String, Integer, Integer, Integer,
                Integer, Integer, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      qswitch
:Pattern:       QSwitch[nl_, orderAlpha_, runAlpha_, orderMass_, runMass_, ord1S_,
                muLam_, xLam_, method_, mZ_, aMz_, mt_, gt_, R_]
:Arguments:     {nl, orderAlpha, runAlpha, orderMass, runMass, ord1S, muLam, xLam,
                method, mZ, aMz, mt, gt, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, Integer, Real, Real,
                String, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      delta1s
:Pattern:       Delta1S[nl_, orderAlpha_, runAlpha_, orderMass_, runMass_,
                muLam_, xLam_, method_, mZ_, aMz_, mt_, R_]
:Arguments:     {nl, orderAlpha, runAlpha, orderMass, runMass, muLam, xLam,
                method, mZ, aMz, mt, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, Real, Real,
                String, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      a1pole
:Pattern:       A1Pole[nl_, order_, En_, mtpole_, gamtop_, asoft_, VcsNNLL_, musoft_]
:Arguments:     {nl, order, En, mtpole, gamtop, asoft, VcsNNLL, musoft}
:ArgumentTypes: {Integer, Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      xinnllmixusoft
:Pattern:       XiNNLLmixUsoft[nl_, ah_, au_]
:Arguments:     {nl, ah, au}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      mllc2
:Pattern:       MLLc2[nl_, ah_, as_]
:Arguments:     {nl, ah, as}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vssll
:Pattern:       VssLL[nl_, ah_, as_]
:Arguments:     {nl, ah, as}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vcsll
:Pattern:       VcsLL[as_]
:Arguments:     {as}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vrsll
:Pattern:       VrsLL[nl_, as_, au_]
:Arguments:     {nl, as, au}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      v2sll
:Pattern:       V2sLL[nl_, ah_, au_, as_]
:Arguments:     {nl, ah, au, as}
:ArgumentTypes: {Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      xinll
:Pattern:       XiNLL[nl_, ah_, au_, as_]
:Arguments:     {nl, ah, au, as}
:ArgumentTypes: {Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vk1sll
:Pattern:       Vk1sLL[nl_, ah_, as_]
:Arguments:     {nl, ah, as}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vk2sll
:Pattern:       Vk2sLL[nl_, ah_, as_]
:Arguments:     {nl, ah, as}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vkeffsll
:Pattern:       VkeffsLL[nl_, ah_, as_]
:Arguments:     {nl, ah, as}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      qfromv
:Pattern:       QFromV[v_, m_, gt_]
:Arguments:     {v, m, gt}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      switchoff
:Pattern:       SwitchOff[q_, m_, gt_, v0_, v1_]
:Arguments:     {q, m, gt, v0, v1}
:ArgumentTypes: {Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vc
:Pattern:       VC[v_, m_, gt_]
:Arguments:     {v, m, gt}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      vstar
:Pattern:       VStar[v_, m_, gt_]
:Arguments:     {v, m, gt}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      vrootstar
:Pattern:       VRootStar[v_, m_, gt_]
:Arguments:     {v, m, gt}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      ttbar
:Pattern:       TTbar[energy_, topmass_, topgamma_, alphas0_, mue0_, cutn_,
                cutv_,  c0_, c1_, c2_, cdeltapotc_, cdeltapot1_, cfullc_,
                cfull1_, crm2_, kincm_, kinca_, ijknflg_, ijgcflg_, kincv_,
                ijvflg_]

:Arguments:     {energy, topmass, topgamma, alphas0, mue0, cutn, cutv, c0, c1,
                 c2, cdeltapotc, cdeltapot1, cfullc, cfull1, crm2, kincm, kinca,
                 ijknflg, ijgcflg, kincv, ijvflg}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Integer, Integer,
                 Real, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ttbarlist
:Pattern:       TTbarList[energy_, topmass_, topgamma_, alphas0_, mue0_, cutn_,
                cutv_,  c0_, c1_, c2_, cdeltapotc_, cdeltapot1_, cfullc_,
                cfull1_, crm2_, kincm_, kinca_, ijknflg_, ijgcflg_, kincv_,
                ijvflg_]

:Arguments:     {energy, topmass, topgamma, alphas0, mue0, cutn, cutv, c0, c1,
                 c2, cdeltapotc, cdeltapot1, cfullc, cfull1, crm2, kincm, kinca,
                 ijknflg, ijgcflg, kincv, ijvflg}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Integer, Integer,
                 Real, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ewfactors
:Pattern:       EWFactors[nf_, Q_, Mz_, GammaZ_, sin2ThetaW_]
:Arguments:     {nf, Q, Mz, GammaZ, sin2ThetaW}
:ArgumentTypes: {Integer, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      legendrelist
:Pattern:       LegendreList[n_, k_, x_]
:Arguments:     {n, k, x}
:ArgumentTypes: {Integer, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      qlegendrelist
:Pattern:       QLegendreList[n_, x_]
:Arguments:     {n, x}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      gammar
:Pattern:       GammaR[str_, nf_]
:Arguments:     {str, nf}
:ArgumentTypes: {String, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      masslessprof
:Pattern:       MasslessProf[terms_, hard_, shape_, setup_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,
                delta0_, h_, tau_]
:Arguments:     {terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,
                 R0, muR0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, RealList, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      findorigin
:Pattern:       FindOrigin[shape_, gap_, orderAlpha_, runAlpha_, order_, run_, nf_,
                mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_,
                n0_, n1_, t2_, tR_, ts_, slope_, cnt_, eH_, eS_, eR_, R0_, muR0_, delta0_,
                h_]
:Arguments:     {shape, gap, orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, muT, mB,
                 muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH,
                 eS, eR, R0, muR0, delta0, h}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      masslessprofpiece
:Pattern:       MasslessProfPiece[terms_, hard_, shape_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,
                muR0_, delta0_, h_, tau_]
:Arguments:     {terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,
                 lambda, R0, muR0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, Integer, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslessprofpiecelist
:Pattern:       MasslessProfPieceList[terms_, hard_, shape_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,
                muR0_, delta0_, h_, tauList_]
:Arguments:     {terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,
                 lambda, R0, muR0, delta0, h, tauList}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, Integer, Real, Real, Real,
                 Real, Real, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslesspiecebin
:Pattern:       MasslessPieceBin[terms_, hard_, shape_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,
                muR0_, delta0_, h_, tauList_]
:Arguments:     {terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,
                 lambda, R0, muR0, delta0, h, Flatten[tauList], Length[tauList]}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, Integer, Real, Real, Real,
                 Real, Real, RealList, Integer}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslessprofdiffpiece
:Pattern:       MasslessProfPiece[terms_, hard_, shape_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,
                muR0_, delta0_, h_, tau_, tau2_]
:Arguments:     {terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,
                 lambda, R0, muR0, delta0, h, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslessproflist
:Pattern:       MasslessProfList[terms_, hard_, shape_, setup_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,
                delta0_, h_, tauList_]
:Arguments:     {terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,
                 R0, muR0, delta0, h, tauList}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, RealList, Real, Real, Real,
                 Real, Real, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslessbinlist
:Pattern:       MasslessBinList[terms_, hard_, shape_, setup_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,
                delta0_, h_, tauList_]
:Arguments:     {terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,
                 R0, muR0, delta0, h, Flatten[tauList], Length[tauList]}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, RealList, Real, Real, Real,
                 Real, Real, RealList, Integer}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslessprofdiff
:Pattern:       MasslessProf[terms_, hard_, shape_, setup_, gap_, space_, cum_,
                orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,
                tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,
                delta0_, h_, tau_, tau2_]
:Arguments:     {terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,
                 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,
                 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,
                 R0, muR0, delta0, h, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, RealList, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      masslessmoment
:Pattern:       MasslessMoment[terms_, hard_, shape_, setup_, gap_, space_, orderAlpha_,
                runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_, tR_, ts_,
                slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_, delta0_,
                h_, tau_, tau2_, pow_]
:Arguments:     {terms, hard, shape, setup, gap, space, orderAlpha, runAlpha, order, run,
                 nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0,
                 Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0,
                 muR0, delta0, h, tau, tau2, pow}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, RealList, Real, Real, Real,
                 Real, Real, Real, Real, Integer}
:ReturnType:     Real
:End:

:Begin:
:Function:      profiles
:Pattern:       Profiles[Q_, mu0_, R0_, n0_, n1_, t2_, tR_, ts_, slope_, cnt_, eH_, eS_,
                eJ_, eR_, ns_, tau_]
:Arguments:     {Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, tau}
:ArgumentTypes: {Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      profilesmass
:Pattern:       ProfilesMass[Q_, beta_, mu0_, delLamb_, R0_, n0_, delta0_, n1_,
                delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_, mass_, muM_, ns_,
                def_, EShape_, tau_]
:Arguments:     {Q, beta, mu0, delLamb, R0, n0, delta0, n1, delta1, t2, ts, slope,
                 cnt, eH, eS, eJ, mass, muM, ns, def, EShape, tau}
:ArgumentTypes: {Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Integer, String, String, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      mctop
:Pattern:       MCtop[shape_, mt_, Q_, n_, k_, x_]
:Arguments:     {shape, mt, Q, n, k, x}
:ArgumentTypes: {String, Real, Real, Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      breitunstable
:Pattern:       BreitUnstable[shape_, mt_, Q_, gamma_, n_, k_, x_]
:Arguments:     {shape, mt, Q, gamma, n, k, x}
:ArgumentTypes: {String, Real, Real, Real, Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltamctop
:Pattern:       DeltaMCtop[shape_, mt_, Q_]
:Arguments:     {shape, mt, Q}
:ArgumentTypes: {String, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      massintersection
:Pattern:       MassIntersection[Q_, beta_, mu0_, delLamb_, n0_, delta0_, n1_,delta1_,
                t2_, ts_, slope_, cnt_, eH_, eS_, eJ_, mass_, muM_, def_, EShape_]
:Arguments:     {Q, beta, mu0, delLamb, n0, delta0, n1, delta1, t2, ts, slope,
                 cnt, eH, eS, eJ, mass, muM, def, EShape}
:ArgumentTypes: {Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, String, String}
:ReturnType:     Manual
:End:

:Begin:
:Function:      nsmasspiece
:Pattern:       NSMassPiece[shape_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_ runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_]
:Arguments:     {shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,
                orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,
                muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0,
                delta0, h, gammaZ, sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real , Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, IntegerList, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      nsmassdiffpiece
:Pattern:       NSMassPiece[shape_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_ runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_]
:Arguments:     {shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,
                orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,
                muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0,
                delta0, h, gammaZ, sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real , Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, IntegerList, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      nsmass
:Pattern:       NSMass[shape_, setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_]
:Arguments:     {shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha,
                order, run, orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC,
                muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c,
                lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      nsmassdiff
:Pattern:       NSMass[shape_, setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_]
:Arguments:     {shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha,
                order, run, orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC,
                muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c,
                lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      hjmnsmass
:Pattern:       HJMNSMass[setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_]
:Arguments:     {setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,
                orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,
                muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, Flatten[Transpose[c]],
                Length[c], lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, RealList, Integer, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      hjmnsmassdiff
:Pattern:       HJMNSMass[setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,
                runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,
                mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,
                muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_]
:Arguments:     {setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,
                orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,
                muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, Flatten[Transpose[c]],
                Length[c], lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, RealList, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singular
:Pattern:       Singular[hard_, shape_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_,
                delta0_, h_, tau_]
:Arguments:     {hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,
                 s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, RealList, Real, Real, Real, Real,
                 Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularlist
:Pattern:       SingularList[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, clen_, lambda_, R0_, mu0_,
                delta0_, h_, tau_]
:Arguments:     {hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,
                 s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, clen, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Integer, Real, Real, Real, Real,
                 Real, Real}
:ReturnType:     Manual
:End:


:Begin:
:Function:      singulardiff
:Pattern:       Singular[hard_, shape_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_,
                delta0_, h_, tau1_, tau2_]
:Arguments:     {hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,
                 s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, c, lambda, R0, mu0, delta0, h, tau1, tau2}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, RealList, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massiveprof
:Pattern:       MassiveProf[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tau_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massorigin
:Pattern:       MassOrigin[shape_, Eshape_, gap_, scheme_, orderAlpha_, runAlpha_,
                orderMass_, runMass_, order_, run_, nf_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, R0_, muR0_, del0_, h_]
:Arguments:     {shape, Eshape, gap, scheme, orderAlpha, runAlpha, orderMass, runMass,
                 order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2,
                 Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope,
                 cnt, eH, eS, eJ, mass, muM, R0, muR0, del0, h}
:ArgumentTypes: {String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massiveprofpiece
:Pattern:       MassiveProfPiece[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tau_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, Integer, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massiveprofpiecelist
:Pattern:       MassiveProfPieceList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tauList_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tauList}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, Integer, Real, Real, Real, Real, Real, Real, Real, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massivepiecebin
:Pattern:       MassivePieceBin[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tauList_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, Flatten[tauList], Length[tauList]}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, Integer, Real, Real, Real, Real, Real, Real, Real, RealList,
                 Integer}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massiveprofdiffpiece
:Pattern:       MassiveProfPiece[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massiveprofdiff
:Pattern:       MassiveProf[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massiveproflist
:Pattern:       MassiveProfList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tauList_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, tauList}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massivebinlist
:Pattern:       MassiveBinList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,
                sin2ThetaW_, tauList_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,
                 sin2ThetaW, Flatten[tauList], Length[tauList]}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, RealList,
                 Integer}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massivemoment
:Pattern:       MassiveMoment[terms_, hard_, shape_, Eshape_, setup_, gap_, space_,
                scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,
                runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,
                Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,
                mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, delta0_, h_, gammaZ_,
                sin2ThetaW_, tau_, tau2_, pow_]
:Arguments:     {terms, hard, shape, Eshape, setup, gap, space, scheme, abs, current,
                 xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,
                 s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,
                 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,
                 eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, mu0, delta0, h, gammaZ,
                 sin2ThetaW, tau, tau2, pow}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, Real, Real, Integer, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Real, RealList, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Integer}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularmass
:Pattern:       SingularMass[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,
                abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,
                width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi,
                 xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3,
                 mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ,
                 muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,
                 sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, Real, Real, Integer, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 RealList, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularmassdiff
:Pattern:       SingularMass[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,
                abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,
                width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_,
                tau2_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi,
                 xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,
                 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH,
                 muJ, muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h,
                 gammaZ, sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String,
                 String, Real, Real, Integer, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 RealList, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massnondist
:Pattern:       MassNonDist[hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, orderAlpha_, runAlpha_, orderMass_, runMass_, order_,
                run_, nf_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_,
                muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_,
                mu_, c_, lambda_, R0_, mu0_, delta0_, h_, tau_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, orderAlpha,
                 runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,
                 mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,
                 muM, mu, c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, String,
                 Integer, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, RealList, Real,
                 Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massnondistdiff
:Pattern:       MassNonDist[hard_, shape_, Eshape_, setup_, gap_, space_, cum_,
                scheme_, orderAlpha_, runAlpha_, orderMass_, runMass_, order_,
                run_, nf_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_,
                muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_,
                mu_, c_, lambda_, R0_, mu0_, delta0_, h_, tau_, tau2_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, orderAlpha,
                 runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,
                 mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,
                 muM, mu, c, lambda, R0, mu0, delta0, h, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, String,
                 Integer, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, RealList, Real,
                 Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massnondistpiece
:Pattern:       MassNonDistPiece[hard_, shape_, Eshape_, gap_, space_, cum_, scheme_,
                 orderAlpha_, runAlpha_, orderMass_, runMass_, order_, run_, nf_, G3_,
                 mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda1_, muLambda2_,
                 Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_, c_, lambda_, R0_, mu0_,
                 delta0_, h_, tau_]
:Arguments:     {hard, shape, Eshape, gap, space, cum, scheme, orderAlpha, runAlpha,
                 orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB,
                 mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,
                 c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer,
                 Integer, Integer, Integer, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, IntegerList, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      massnondistdiffpiece
:Pattern:       MassNonDistPiece[hard_, shape_, Eshape_, gap_, space_, cum_, scheme_,
                 orderAlpha_, runAlpha_, orderMass_, runMass_, order_, run_, nf_, G3_,
                 mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda1_, muLambda2_,
                 Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_, c_, lambda_, R0_, mu0_,
                 delta0_, h_, tau_, tau2_]
:Arguments:     {hard, shape, Eshape, gap, space, cum, scheme, orderAlpha, runAlpha,
                 orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB,
                 mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,
                 c, lambda, R0, mu0, delta0, h, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, Integer,
                 Integer, Integer, Integer, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, IntegerList, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularmasspiece
:Pattern:       SingularMassPiece[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,
                abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,
                width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB,
                 orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ,
                 amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS,
                 R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,
                 sin2ThetaW, tau}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String, String,
                 Real, Real, Integer, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, IntegerList,
                 Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularmassdiffpiece
:Pattern:       SingularMassPiece[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,
                abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,
                order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,
                muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,
                width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_,
                tau2_]
:Arguments:     {hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB,
                 orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ,
                 amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS,
                 R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,
                 sin2ThetaW, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, String, String, String, String, String,
                 Real, Real, Integer, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, IntegerList,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularhjm
:Pattern:       SingularHJM[hard_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,
                run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,
                R0_, mu0_, delta0_, h_, tau_]
:Arguments:     {hard, setup, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3,
                 s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH,
                 muJ, muS, R, mu, Flatten[Transpose[c]], Length[c], lambda, R0, mu0,
                 delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 RealList, Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularhjm1d
:Pattern:       SingularHJM1D[hard_, gap_, cum_, orderAlpha_, runAlpha_, order_,
                run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_,
                muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,
                R0_, mu0_, delta0_, h_, tau_]
:Arguments:     {hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32,
                 G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, RealList, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularhjm1dpiece
:Pattern:       SingularHJM1DPiece[hard_, gap_, cum_, orderAlpha_, runAlpha_, order_, run_,
                isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_,
                mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_,
                mu0_, delta0_, h_, tau_]
:Arguments:     {hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3,
                 s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 IntegerList, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singulardouble
:Pattern:       SingularDouble[hard_, setup_, gap_, space_, cum1_, cum2_, orderAlpha_, runAlpha_,
                order_, run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,
                R0_, mu0_, delta0_, h_, rho1_, rho2_]
:Arguments:     {hard, setup, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft,
                 nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH,
                 muJ, muS, R, mu, Flatten[Transpose[c]], Length[c], lambda, R0, mu0,
                 delta0, h, rho1, rho2}
:ArgumentTypes: {String, String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 RealList, Integer, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularhjmpiece
:Pattern:       SingularHJMPiece[hard_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_, run_,
                isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,
                muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,
                h_, tau_]
:Arguments:     {hard, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31,
                 s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R,
                 mu, c, lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, Integer, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, IntegerList, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singulardoublepiece
:Pattern:       SingularDoublePiece[hard_, gap_, space_, cum1_, cum2_, orderAlpha_, runAlpha_,
                order_, run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,
                R0_, mu0_, delta0_, h_, rho1_, rho2_]
:Arguments:     {hard, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3,
                 s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,
                 R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2}
:ArgumentTypes: {String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 IntegerList, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singularpiece
:Pattern:       SingularPiece[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,
                run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,
                muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,
                h_, tau_]
:Arguments:     {hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3,
                 mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c,
                 lambda, R0, mu0, delta0, h, tau}
:ArgumentTypes: {String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, IntegerList, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      singulardiffpiece
:Pattern:       SingularPiece[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,
                run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,
                muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,
                h_, tau_, tau2_]
:Arguments:     {hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3,
                 mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c,
                 lambda, R0, mu0, delta0, h, tau, tau2}
:ArgumentTypes: {String, String, String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, IntegerList, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      diffdeltagap
:Pattern:       DiffDeltaGap[str_, scheme_, order_, R0_, R1_, mu0_, mu1_, muLambda_,
                orderAlpha_, runAlpha_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_]
:Arguments:     {str, scheme, order, R0, R1, mu0, mu1, muLambda, orderAlpha, runAlpha,
                 nf, Mz, aMz, mT, muT, mB, muB, mC, muC}
:ArgumentTypes: {String, String, Integer, Real, Real, Real, Real, Real, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      diffdeltagapmass
:Pattern:       DiffDeltaGapMass[str_, order_, R0_, R1_, mu0_, mu1_, muM_, muLambda1_,
                muLambda2_, orderAlpha_, runAlpha_, runMass_, nf_, Mz_, aMz_, mT_, muT_,
                mB_, muB_, mC_, muC_]
:Arguments:     {str, order, R0, R1, mu0, mu1, muM, muLambda1, muLambda2, orderAlpha,
                 runAlpha, runMass, nf, Mz, aMz, mT, muT, mB, muB, mC, muC}
:ArgumentTypes: {String, Integer, Real, Real, Real, Real, Real, Real, Real, Integer,
                 Integer, Integer, Integer, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      hyperf32exact
:Pattern:       HyperF32Exact[w_, x_]
:Arguments:     {w, x}
:ArgumentTypes: {Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      model
:Pattern:       Model[c_, lambda_, k_, l_]
:Arguments:     {c, lambda, k, l}
:ArgumentTypes: {RealList, Real, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      modelunstable
:Pattern:       ModelUnstable[shape_, mt_, Q_, c_, lambda_, n_, k_, l_]
:Arguments:     {shape, mt, Q, c, lambda, n, k, l}
:ArgumentTypes: {String, Real, Real, RealList, Real, Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      breitmodelunstable
:Pattern:       BreitModelUnstable[shape_, mt_, Q_, gamma_, c_, lambda_, n_, k_, l_]
:Arguments:     {shape, mt, Q, gamma, c, lambda, n, k, l}
:ArgumentTypes: {String, Real, Real, Real, RealList, Real, Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      modelunstablediff
:Pattern:       ModelUnstable[shape_, mt_, Q_, c_, lambda_, n_, k_, l_, l2_]
:Arguments:     {shape, mt, Q, c, lambda, n, k, l, l2}
:ArgumentTypes: {String, Real, Real, RealList, Real, Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      modeldiff
:Pattern:       Model[c_, lambda_, k_, l_, l2_]
:Arguments:     {c, lambda, k, l, l2}
:ArgumentTypes: {RealList, Real, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      breitmodel
:Pattern:       BreitModel[c_, lambda_, width_, k_, l_]
:Arguments:     {c, lambda, width, k, l}
:ArgumentTypes: {RealList, Real, Real, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      breitmodeldiff
:Pattern:       BreitModel[c_, lambda_, width_, k_, l_, l2_]
:Arguments:     {c, lambda, width, k, l, l2}
:ArgumentTypes: {RealList, Real, Real, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      taylor
:Pattern:       Taylor[c_, lambda_, k_]
:Arguments:     {c, lambda, k}
:ArgumentTypes: {RealList, Real, Integer}
:ReturnType:    Real
:End:

:Begin:
:Function:      momentmodel
:Pattern:       MomentModel[c_, lambda_, k_]
:Arguments:     {c, lambda, k}
:ArgumentTypes: {RealList, Real, Integer}
:ReturnType:    Real
:End:

:Begin:
:Function:      modelpiece
:Pattern:       ModelPiece[c_, lambda_, k_, l_]
:Arguments:     {c, lambda, k, l}
:ArgumentTypes: {IntegerList, Real, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      taylorpiece
:Pattern:       TaylorPiece[c_, lambda_, k_]
:Arguments:     {c, lambda, k}
:ArgumentTypes: {IntegerList, Real, Integer}
:ReturnType:    Real
:End:

:Begin:
:Function:      hyper2f1
:Pattern:       Hyper2F1[a_, b_, c_, x_]
:Arguments:     {a, b, c, x}
:ArgumentTypes: {Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      intecorre
:Pattern:       InteCorre[b_, x0_, x1_]
:Arguments:     {b, x0, x1}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      anomdim
:Pattern:       AnomDim[str_, nf_, G4_]
:Arguments:     {str, nf, G4}
:ArgumentTypes: {String, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      afroms
:Pattern:       aFromS[str_, nf_, G4_]
:Arguments:     {str, nf, G4}
:ArgumentTypes: {String, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ql
:Pattern:       Ql[str_, nf_, order_, G4_]
:Arguments:     {str, nf, order, G4}
:ArgumentTypes: {String, Integer, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      anasy
:Pattern:       aNasy[str_, nf_, order_, G4_]
:Arguments:     {str, nf, order, G4}
:ArgumentTypes: {String, Integer, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      anlambda
:Pattern:       anLambda[str_, nf_, G4_]
:Arguments:     {str, nf, G4}
:ArgumentTypes: {String, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      msbardeltapiece
:Pattern:       MSbarDeltaPiece[nl_, nh_]
:Arguments:     {nl, nh}
:ArgumentTypes: {Integer, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      alphamatchinglog
:Pattern:       AlphaMatchingLog[str_, direction_, nf_]
:Arguments:     {str, direction, nf}
:ArgumentTypes: {String, String, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      alphamatchinginverse
:Pattern:       AlphaMatchingInverse[str_, nf_]
:Arguments:     {str, nf}
:ArgumentTypes: {String, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ccoef
:Pattern:       cCoef[nf_, order_, n_]
:Arguments:     {nf, order, n}
:ArgumentTypes: {Integer, Integer, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pscoef
:Pattern:       PSCoef[nf_, lg_]
:Arguments:     {nf, lg}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      n12
:Pattern:       N12[str_, order_, nf_, Lambda_, err_]
:Arguments:     {str, order, nf, Lambda, err}
:ArgumentTypes: {String, Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      n12ratio
:Pattern:       N12Ratio[str_, order_, nf_, Lambda_, err_]
:Arguments:     {str, order, nf, Lambda, err}
:ArgumentTypes: {String, Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      n12residue
:Pattern:       N12Residue[str_, order_, nf_, Lambda_, err_]
:Arguments:     {str, order, nf, Lambda, err}
:ArgumentTypes: {String, Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      n12rs
:Pattern:       N12RS[str_, order_, n_, nf_, Lambda_, err_]
:Arguments:     {str, order, n, nf, Lambda, err}
:ArgumentTypes: {String, Integer, Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      n12generic
:Pattern:       N12Generic[aCoef_, order_, nf_, Lambda_]
:Arguments:     {aCoef, order, nf, Lambda}
:ArgumentTypes: {RealList, Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      scoef
:Pattern:       sCoef[str_, nf_]
:Arguments:     {str, nf}
:ArgumentTypes: {String, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      scoefgamma
:Pattern:       sCoefGamma[gama_, nf_]
:Arguments:     {gama, nf}
:ArgumentTypes: {RealList, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      scoeflambda
:Pattern:       sCoefLambda[str_, nf_, lambda_]
:Arguments:     {str, nf, lambda}
:ArgumentTypes: {String, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      polylog
:Pattern:       Polylog[n_, z_]
:Arguments:     {n, z}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      dilog
:Pattern:       DiLog[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      pfq
:Pattern:       pFq[a_, b_, z_]
:Arguments:     {a, b, z}
:ArgumentTypes: {RealList, RealList, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      elliptic3
:Pattern:       Elliptic3[psi_, k_, c_]
:Arguments:     {psi, k, c}
:ArgumentTypes: {Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      nglfunction
:Pattern:       NGLFunction[n_, z_]
:Arguments:     {n, z}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      nglsoft
:Pattern:       NGLSoft[n_, z_]
:Arguments:     {n, Re[z], Im[z]}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      complexpolylog
:Pattern:       ComplexPolylog[n_, z_]
:Arguments:     {n, Re[z], Im[z]}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      cli2
:Pattern:       CLi2[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      upsilondeltacharm
:Pattern:       UpsilonDeltaCharm[n_, l_, alpha_, mb_, mc_]
:Arguments:     {n, l, alpha, mb, mc}
:ArgumentTypes: {Integer, Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmexact
:Pattern:       DeltaCharmExact[charm_, type_, scheme_, average_, n_, l_, j_,
                s_, nl_, mH_, mL_, mu_, alp_]
:Arguments:     {charm, type, scheme, average, n, l, j, s, nl, mH, mL, mu, alp}
:ArgumentTypes: {String, String, String, String, Integer, Integer, Integer,
                Integer, Integer, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      upsilondeltacharmbin
:Pattern:       UpsilonDeltaCharmBin[n_, l_, alpha_, mb_, mc_]
:Arguments:     {n, l, alpha, mb, mc}
:ArgumentTypes: {Integer, Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharm3
:Pattern:       DeltaCharm3[nl_, nh_, z_]
:Arguments:     {nl, nh, z}
:ArgumentTypes: {Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharm3der
:Pattern:       DeltaCharm3Der[nl_, nh_, z_]
:Arguments:     {nl, nh, z}
:ArgumentTypes: {Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      gammarcharm3
:Pattern:       GammaRCharm3[nl_, nh_, z_]
:Arguments:     {nl, nh, z}
:ArgumentTypes: {Integer, Integer, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharm2
:Pattern:       DeltaCharm2[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      p2
:Pattern:       P2[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      pi0
:Pattern:       Pi0[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi0der
:Pattern:       Pi0Der[i_, z_]
:Arguments:     {i, Re[z], Im[z]}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi1der
:Pattern:       Pi1Der[i_, z_]
:Arguments:     {i, Re[z], Im[z]}
:ArgumentTypes: {Integer, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi1
:Pattern:       Pi1[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi3
:Pattern:       Pi3[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi2
:Pattern:       Pi2[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      pi2der
:Pattern:       Pi2Der[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      p2int
:Pattern:       P2Int[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltabottomcharm
:Pattern:       DeltaBottomCharm[z1_, z2_]
:Arguments:     {z1, z2}
:ArgumentTypes: {Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      gammarbottomcharm
:Pattern:       GammaRBottomCharm[z1_, z2_]
:Arguments:     {z1, z2}
:ArgumentTypes: {Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      p2double
:Pattern:       P2Double[z1_, z2_]
:Arguments:     {z1, z2}
:ArgumentTypes: {Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmnh
:Pattern:       DeltaCharmNh[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmglue
:Pattern:       DeltaCharmGlue[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmglueder
:Pattern:       DeltaCharmGlueDer[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmnl
:Pattern:       DeltaCharmNl[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmnhder
:Pattern:       DeltaCharmNhDer[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharmnlder
:Pattern:       DeltaCharmNlDer[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      gammarcharm2
:Pattern:       GammaRCharm2[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltacharm2der
:Pattern:       DeltaCharm2Der[z_]
:Arguments:     {z}
:ArgumentTypes: {Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      thrustns1loop
:Pattern:       ThrustNS1loop[t_]
:Arguments:     {t}
:ArgumentTypes: {Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      fomass
:Pattern:       FOMass[shape_, current_, m_, Q_, Mz_, gammaZ_, sin2ThetaW_, t_]
:Arguments:     {shape, current, m, Q, Mz, gammaZ, sin2ThetaW, t}
:ArgumentTypes: {String, String, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      thrustns2loop
:Pattern:       ThrustNS2loop[er_, t_]
:Arguments:     {er, t}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      cli3
:Pattern:       CLi3[z_]
:Arguments:     {Re[z], Im[z]}
:ArgumentTypes: {Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      delta
:Pattern:       Delta[str_, nf_, mu_, R_]
:Arguments:     {str, nf, mu, R}
:ArgumentTypes: {String, Integer, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      deltagap
:Pattern:       DeltaGap[str_, orderAlpha_, runAlpha_, runMass_, nf_, mZ_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, mu_, R_]
:Arguments:     {str, orderAlpha, runAlpha, runMass, nf, mZ, aMz, mT, muT, mB, muB, mC,
                 muC, mu, R}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      psdelta
:Pattern:       PSDelta[orderAlpha_, runAlpha_, nf_, mZ_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, mu_, R_, lg_]
:Arguments:     {orderAlpha, runAlpha, nf, mZ, aMz, mT, muT, mB, muB,
                 mC, muC, mu, R, lg}
:ArgumentTypes: {Integer, Integer, Integer, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      coefmat
:Pattern:       CoefMat[str_, nf_, s3_]
:Arguments:     {str, nf, s3}
:ArgumentTypes: {String, Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      wtilde
:Pattern:       wTilde[order_, nf_, gamma_, a0_, a1_]
:Arguments:     {order, nf, gamma, a0, a1}
:ArgumentTypes: {Integer, Integer, RealList, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      ktilde
:Pattern:       kTilde[order_, nf_, gamma_, a0_, a1_]
:Arguments:     {order, nf, gamma, a0, a1}
:ArgumentTypes: {Integer, Integer, RealList, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      alphaqcd
:Pattern:       AlphaQCD[str_, method_, order_, run_, nf_, Mz_, aMz_, mT_, muT_,
                mB_, muB_, mC_, muC_, mu_]
:Arguments:     {str, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,
                 muC, mu}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      alphaqed
:Pattern:       AlphaQED[nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_, mu_]
:Arguments:     {nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}
:ArgumentTypes: {Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      alphacomplex
:Pattern:       AlphaComplex[str_, method_, order_, run_, nf_, Mz_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, mu_]
:Arguments:     {str, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,
                muC, Re[mu], Im[mu]}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      msbarmass
:Pattern:       MSbarMass[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,
                mC_, muC_, mu_]
:Arguments:     {order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real, Real,
                Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      polemass
:Pattern:       PoleMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_, muT_,
                mB_, muB_, mC_, muC_, mu_]
:Arguments:     {orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC,
                 mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real,
                Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      msbarmasslow
:Pattern:       MSbarMassLow[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_,
                muB_, mC_, muC_, mu_]
:Arguments:     {order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real, Real,
                Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      msrmass
:Pattern:       MSRMass[type_, method_, orderAlpha_, runAlpha_, order_, run_,
                nf_, Mz_, aMz_, mT_,muT_, mB_, muB_, mC_, muC_, lambda_, mu_,
                R_]
:Arguments:     {type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB,
                 mC, muC, lambda, mu, R}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      rsmass
:Pattern:       RSMass[type_, method_, orderAlpha_, runAlpha_, order_, run_, nf_, Mz_,
                aMz_, mT_, muT_, mB_, muB_, mC_, muC_, lambda_, R_]
:Arguments:     {type, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT,
                 mB, muB, mC, muC, lambda, R}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      msrvfns
:Pattern:       MSRVFNS[up_, type_, method_, orderAlpha_, runAlpha_, order_,
                run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_, lambda_,
                mu1_, mu2_, R_]
:Arguments:     {up, type, method, orderAlpha, runAlpha, order, run, nf, Mz,
                 aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, R}
:ArgumentTypes: {String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      msrtop
:Pattern:       MSRTop[up_, type_, method_, orderAlpha_, runAlpha_, order_,
                run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_, lambda_,
                mu1_, mu2_, mu3_, R_]
:Arguments:     {up, type, method, orderAlpha, runAlpha, order, run, nf, Mz,
                 aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, mu3, R}
:ArgumentTypes: {String, String, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      nrqcd
:Pattern:       NRQCD[n_, l_, j_, s_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu_,
                R_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      nrqcddercharm
:Pattern:       NRQCDDerCharm[n_, l_, j_, s_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu_,
                R_, eps_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu, R, eps}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      nrqcdderalpha
:Pattern:       NRQCDDerAlpha[n_, l_, j_, s_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu_,
                R_, eps_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu, R, eps}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massiter
:Pattern:       MassIter[n_, l_, j_, s_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_,
                lam_, mu_, R_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String, String, String,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      massexpand
:Pattern:       MassExpand[n_, l_, j_, s_, charm_, scheme_, average_, method_, counting_, orderAlpha_,
                runAlpha_, order_, run_, nl_, mZ_, amZ_, mT_, muT_, mB_, muB_,
                mC_, muC_, mass_, lambda1_, lambda2_, lam_, mu_, R_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String, String, String,
                 Integer, Integer, Integer, Integer, Integer, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      findmass
:Pattern:       FindMass[ord_, n_, l_, j_, s_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_, mT_, muT_,
                mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_, lam_, mu_, R_]
:Arguments:     {ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 mass, lambda1, lambda2, lam, mu, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, String, String, String,
                 String, String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      findenergy
:Pattern:       FindEnergy[ord_, n_, l_, j_, s_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_, mT_, muT_,
                mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu_, R_]
:Arguments:     {ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu, R}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, String, String, String,
                 String, String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      masserror
:Pattern:       MassError[ord_, n_, l_, j_, s_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_, mT_, muT_,
                mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, x_]
:Arguments:     {ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR,
                 x}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, String, String, String,
                 String, String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      masslist
:Pattern:       MassList[ord_, n_, l_, j_, s_, iter_, charm_, scheme_, average_,
                method_, counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_,
                lam_, mu0_, mu1_, deltaMu_, R0_, R1_, deltaR_]
:Arguments:     {ord, n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, String, String, String,
                 String, String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      nrqcdlist
:Pattern:       NRQCDList[n_, l_, j_, s_, iter_, charm_, scheme_, average_,
                method_, counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_,
                lam_, mu0_, mu1_, deltaMu_, R0_, R1_, deltaR_]
:Arguments:     {n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      upsilonlist
:Pattern:       UpsilonList[n_, l_, j_, s_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, epsAlpha_, epsCharm_]
:Arguments:     {n, l, j, s, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR,
                 epsAlpha, epsCharm}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      corrmat
:Pattern:       CorrMat[qnlist_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, epsAlpha_, epsCharm_]
:Arguments:     {Flatten[qnlist], Length[qnlist], charm, scheme, average, method,
                 counting, orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT,
                 muT, mB, muB, mC, muC, lambda1, lambda2, lam, mu0, mu1, deltaMu,
                 R0, R1, deltaR, epsAlpha, epsCharm}
:ArgumentTypes: {IntegerList, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      chi2nrqcd
:Pattern:       Chi2NRQCD[qnlist_, datalist_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, n_, nl_, mZ_,
                amZ_, mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_,
                mu_, R_]
:Arguments:     {Flatten[qnlist], Flatten[datalist], Length[qnlist], iter, charm,
                 scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {IntegerList, RealList, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, RealList, RealList}
:ReturnType:     Real
:End:

:Begin:
:Function:      chi2minnrqcd
:Pattern:       Chi2MinNRQCD[qnlist_, datalist_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, n_, nl_, mZ_,
                amZ_, mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_,
                mu_, R_]
:Arguments:     {Flatten[qnlist], Flatten[datalist], Length[qnlist], iter, charm,
                 scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {IntegerList, RealList, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, RealList, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      chi2minalphanrqcd
:Pattern:       Chi2MinAlphaNRQCD[qnlist_, datalist_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, n_, nl_, mZ_,
                amZ_, mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_,
                mu_, R_]
:Arguments:     {Flatten[qnlist], Flatten[datalist], Length[qnlist], iter, charm,
                 scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {IntegerList, RealList, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, RealList, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      chi2minalphambnrqcd
:Pattern:       Chi2MinAlphaMbNRQCD[qnlist_, datalist_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, n_, nl_, mZ_,
                amZ_, mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_,
                mu_, R_]
:Arguments:     {Flatten[qnlist], Flatten[datalist], Length[qnlist], iter, charm,
                 scheme, average, method, counting, orderAlpha, runAlpha, order,
                 run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1,
                 lambda2, lam, mu, R}
:ArgumentTypes: {IntegerList, RealList, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, RealList, RealList}
:ReturnType:     Manual
:End:

:Begin:
:Function:      errmat
:Pattern:       ErrMat[qnlist_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, epsAlpha_, epsCharm_]
:Arguments:     {Flatten[qnlist], Length[qnlist], charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR,
                 epsAlpha, epsCharm}
:ArgumentTypes: {IntegerList, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      errmatrices
:Pattern:       ErrMatrices[qnlist_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_,
                mT_, muT_, mB_, muB_, mC_, muC_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, epsAlpha_, epsCharm_]
:Arguments:     {Flatten[qnlist], Length[qnlist], charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR,
                 epsAlpha, epsCharm}
:ArgumentTypes: {IntegerList, Integer, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      nrqcderror
:Pattern:       NRQCDError[n_, l_, j_, s_, iter_, charm_, scheme_, average_, method_,
                counting_, orderAlpha_, runAlpha_, order_, run_, nl_, mZ_, amZ_, mT_, muT_,
                mB_, muB_, mC_, muC_, mass_, lambda1_, lambda2_, lam_, mu0_,
                mu1_, deltaMu_, R0_, R1_, deltaR_, x_]
:Arguments:     {n, l, j, s, iter, charm, scheme, average, method, counting, orderAlpha,
                 runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,
                 mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR,
                 x}
:ArgumentTypes: {Integer, Integer, Integer, Integer, String, String, String, String,
                 String, String, Integer, Integer, Integer, Integer, Integer,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      optimalr
:Pattern:       OptimalR[type_, n_, method_, orderAlpha_, runAlpha_, order_, run_,
                nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_, lambda_, mu_]
:Arguments:     {type, n, method, orderAlpha, runAlpha, order, run, nf, Mz, aMz,
                 mT, muT, mB, muB, mC, muC, lambda, mu}
:ArgumentTypes: {String, Real, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      optimalrvfns
:Pattern:       OptimalRVFNS[up_, type_, n_, method_, orderAlpha_, runAlpha_,
                order_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,
                lambda_, mu1_, mu2_]
:Arguments:     {up, type, n, method, orderAlpha, runAlpha, order, run, nf, Mz,
                 aMz, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2}
:ArgumentTypes: {String, String, Real, String, Integer, Integer, Integer, Integer,
                 Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      optimalr2
:Pattern:       OptimalR2[n_, orderAlpha_, runAlpha_, order_, run_, nf_, Mz_,
                aMz_, mT_, muT_, mB_, muB_, mC_, muC_, mass_]
:Arguments:     {n, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB,
                 muB, mC, muC, mass}
:ArgumentTypes: {Real, Integer, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      mmfrommsr
:Pattern:       mmfromMSR[type_, orderAlpha_, runAlpha_, order_, run_, nf_, Mz_,
                aMz_, mT_, muT_, mB_, muB_, mC_, muC_, mu_, R_]
:Arguments:     {type, orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT,
                 mB, muB, mC, muC, mu, R}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      jetmass
:Pattern:       JetMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, R_, mu_]
:Arguments:     {orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,
                 muC, muLambda, R, mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      mmfromjetmass
:Pattern:       mmFromJetMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_,
                muT_, mB_, muB_, mC_, muC_, muLambda_, R_, mu_]
:Arguments:     {orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,
                 muC, muLambda, R, mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      deltamsbar
:Pattern:       DeltaMSbar[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,
                 mC_, muC_, mu_]
:Arguments:     {order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}
:ArgumentTypes: {Integer, Integer, Integer, Integer, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:     Manual
:End:

:Begin:
:Function:      rhad
:Pattern:       Rhad[scheme_, orderAlpha_, runAlpha_, order_, nf_, Mz_, aMz_,
                mT_, muT_, mB_, muB_, mC_, muC_, mu_, Q_]
:Arguments:     {scheme, orderAlpha, runAlpha, order, nf, Mz, aMz, mT, muT, mB,
                 muB, mC, muC, mu, Q}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmahad
:Pattern:       SigmaHad[scheme_, current_, orderAlpha_, runAlpha_, order_,
                nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_, muT_, mB_,
                muB_, mC_, muC_, mu_, Q_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ,
                 sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, mu, Q}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmarad
:Pattern:       SigmaRad[scheme_, current_, orderAlpha_, runAlpha_, order_,
                nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_, muT_, mB_,
                muB_, mC_, muC_, eH_, Q_, x_, theta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ,
                 sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x,
                 theta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmaradcum
:Pattern:       SigmaRadCum[scheme_, current_, orderAlpha_, runAlpha_, order_,
                nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_, muT_, mB_,
                muB_, mC_, muC_, eH_, Q_, x0_, x1_, theta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ,
                 sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x0,
                 x1, theta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmaradcone
:Pattern:       SigmaRadCone[scheme_, current_, orderAlpha_, runAlpha_, order_,
                nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_, muT_, mB_,
                muB_, mC_, muC_, eH_, Q_, x_, theta_, deltaTheta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ,
                 sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x,
                 theta, deltaTheta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmaradconecum
:Pattern:       SigmaRadConeCum[scheme_, current_, orderAlpha_, runAlpha_, order_,
                nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_, muT_, mB_,
                muB_, mC_, muC_, eH_, Q_, x0_, x1_, theta_, deltaTheta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, order, nf, Mz, GammaZ,
                 sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH, Q, x0,
                 x1, theta, deltaTheta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Real, Real,
                 Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rhadcoefs
:Pattern:       RhadCoefs[nf_]
:Arguments:     {nf}
:ArgumentTypes: {Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      rhadmass
:Pattern:       RhadMass[scheme_, current_, orderAlpha_, runAlpha_, runMass_,
                order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, mT_, muT_, mB_,
                muB_, mC_, muC_, mu_, Q_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, mu, Q}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamass
:Pattern:       SigmaMass[scheme_, current_, orderAlpha_, runAlpha_, runMass_,
                order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_, mT_,
                muT_, mB_, muB_, mC_, muC_, mu_, Q_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, mu,
                Q}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamassrad
:Pattern:       SigmaMassRad[scheme_, current_, orderAlpha_, runAlpha_,
                runMass_, order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_,
                mT_, muT_, mB_, muB_, mC_, muC_, eH_, Q_, x_, theta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH,
                Q, x, theta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamassradcum
:Pattern:       SigmaMassRadCum[scheme_, current_, orderAlpha_, runAlpha_,
                runMass_, order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_,
                mT_, muT_, mB_, muB_, mC_, muC_, eH_, Q_, x0_, x1_, theta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH,
                Q, x0, x1, theta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamassradcone
:Pattern:       SigmaMassRadCone[scheme_, current_, orderAlpha_, runAlpha_,
                runMass_, order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_,
                mT_, muT_, mB_, muB_, mC_, muC_, eH_, Q_, x_, theta_, deltaTheta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH,
                Q, x, theta, deltaTheta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamassradconecum
:Pattern:       SigmaMassRadConeCum[scheme_, current_, orderAlpha_, runAlpha_,
                runMass_, order_, nf_, Mz_, GammaZ_, sin2ThetaW_, aMz_, aMzQED_,
                mT_, muT_, mB_, muB_, mC_, muC_, eH_, Q_, x0_, x1_, theta_,
                deltaTheta_]
:Arguments:     {scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz,
                GammaZ, sin2ThetaW, aMz, aMzQED, mT, muT, mB, muB, mC, muC, eH,
                Q, x0, x1, theta, deltaTheta}
:ArgumentTypes: {String, String, Integer, Integer, Integer, Integer, Integer,
                Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rqcd
:Pattern:       RQCD[scheme_, runAlpha_, runMass_, ordMass_, order_, ord1S_,
                R1S_, method_, lambda_, gt_, Mz_, aMz_, mT_, mu_, Q_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, aMz, mT, mu, Q}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rexp
:Pattern:       RExp[scheme_, runAlpha_, runMass_, ordMass_, order_, ord1S_,
                R1S_, method_, lambda_, gt_, Mz_, aMz_, mT_, mu_, nu_, Q_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, aMz, mT, mu, nu, Q}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rmatched
:Pattern:       Rmatched[scheme_, runAlpha_, runMass_, ordMass_, order_, ord1S_,
                R1S_, method_, lambda_, gt_, Mz_, aMz_, mT_, mu_, nu_, v1_, v2_,
                Q_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, aMz, mT, mu, nu, v1, v2, Q}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamatched
:Pattern:       SigmaMatched[scheme_, runAlpha_, runMass_, ordMass_, order_, ord1S_,
                R1S_, method_, lambda_, gt_, Mz_, gammaZ_, sinW_, aMz_, aMzQED_,
                mT_, mu_, nu_, v1_, v2_, Q_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, nu, v1, v2, Q}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamatchedrad
:Pattern:       SigmaMatchedRad[scheme_, runAlpha_, runMass_, ordMass_,
                order_, ord1S_, R1S_, method_, lambda_, gt_, Mz_, gammaZ_,
                sinW_, aMz_, aMzQED_, mT_, mu_, nu_, v1_, v2_, Q_, x_, theta_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, nu, v1, v2,
                Q, x, theta}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamatchedradcum
:Pattern:       SigmaMatchedRadCum[scheme_, runAlpha_, runMass_, ordMass_,
                order_, ord1S_, R1S_, method_, lambda_, gt_, Mz_, gammaZ_,
                sinW_, aMz_, aMzQED_, mT_, mu_, nu_, v1_, v2_, Q_, x0_, x1_,
                theta_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, nu, v1, v2,
                Q, x0, x1, theta}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamatchedradcone
:Pattern:       SigmaMatchedRadCone[scheme_, runAlpha_, runMass_, ordMass_,
                order_, ord1S_, R1S_, method_, lambda_, gt_, Mz_, gammaZ_,
                sinW_, aMz_, aMzQED_, mT_, mu_, nu_, v1_, v2_, Q_, x_, theta_,
                deltaTheta_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, nu, v1, v2,
                Q, x, theta, deltaTheta}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      sigmamatchedradconecum
:Pattern:       SigmaMatchedRadConeCum[scheme_, runAlpha_, runMass_, ordMass_,
                order_, ord1S_, R1S_, method_, lambda_, gt_, Mz_, gammaZ_,
                sinW_, aMz_, aMzQED_, mT_, mu_, nu_, v1_, v2_, Q_, x0_, x1_,
                theta_,deltaTheta_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, gammaZ, sinW, aMz, aMzQED, mT, mu, nu, v1, v2,
                Q, x0, x1, theta, deltaTheta}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      rmatchedlist
:Pattern:       RmatchedList[scheme_, runAlpha_, runMass_, ordMass_, order_, ord1S_,
                R1S_, method_, lambda_, gt_, Mz_, aMz_, mT_, mu_, nu_, v1_, v2_,
                Q0_, Q1_, deltaQ_]
:Arguments:     {scheme, runAlpha, runMass, ordMass, order, ord1S, R1S, method,
                lambda, gt, Mz, aMz, mT, mu, nu, v1, v2, Q0, Q1, deltaQ}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Integer, Real,
                String, Real, Real, Real, Real, Real, Real, Real, Real, Real,
                Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      lambdaqcd
:Pattern:       LambdaQCD[scheme_, order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_,
                 mB_, muB_, mC_, muC_, mu_]
:Arguments:     {scheme, order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}
:ArgumentTypes: {String, Integer, Integer, Integer, Integer, Real, Real, Real, Real,
                 Real, Real, Real, Real, Real}
:ReturnType:     Real
:End:

:Begin:
:Function:      kernels
:Pattern:       Kernel[n_, width_, w_, mu_, p_]
:Arguments:     {n, width, w, mu, p}
:ArgumentTypes: {Integer, Real, Real, Real, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      gammaderlist
:Pattern:       GammaDerList[n_, w_]
:Arguments:     {n, w}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      polygamma
:Pattern:       polyGamma[n_, w_]
:Arguments:     {n, w}
:ArgumentTypes: {Integer, Real}
:ReturnType:    Manual
:End:

:Begin:
:Function:      nglkernels
:Pattern:       NGLKernel[n_, n1_, n2_, width_, w_, mu_, p_]
:Arguments:     {n, n1, n2, width, w, mu, p}
:ArgumentTypes: {Integer, Integer, Integer, Real, RealList, Real, RealList}
:ReturnType:    Manual
:End:

:Begin:
:Function:      nglintegral
:Pattern:       NGLIntegral[nf_, pow_, w1_, w2_]
:Arguments:     {nf, pow, w1, w2}
:ArgumentTypes: {Integer, Integer, Real, Real}
:ReturnType:    Real
:End:

:Begin:
:Function:      ngldoubleintegral
:Pattern:       NGLDoubleIntegral[nf_, pow_, w1_, w2_, r_]
:Arguments:     {nf, pow, w1, w2, r}
:ArgumentTypes: {Integer, Integer, Real, Real, Real}
:ReturnType:    Real
:End:

:Evaluate:   realQ = Head[# + 1.] === Real &

:Evaluate:   End[]

:Evaluate:   EndPackage[]

#include "mathlink.h"
#include "ftypes.h"
#include <stdio.h>
#include <unistd.h>
#include <math.h>

extern double f90toppik_(double* energy, double* tm, double* tg,
double* alphas, double* scale, double* cutn, double* cutv, double* c0,
double* c1, double* c2, double* cdeltc, double* cdeltl, double* cfullc,
double* cfulll, double* crm2, double* kincm, double* kinca, int* jknflg,
int* jgcflg, double* kincv, int* jvflg, double* res);

static void ttbar(double energy,double tm, double tg, double alphas,
double scale, double cutn, double cutv, double c0, double c1, double c2,
double cdeltc, double cdeltl, double cfullc, double cfulll, double crm2,
double kincm, double kinca, int jknflg, int jgcflg, double kincv, int jvflg){
double res[2];

 f90toppik_(&energy, &tm, &tg, &alphas, &scale, &cutn, &cutv, &c0, &c1,
 &c2, &cdeltc, &cdeltl, &cfullc, &cfulll, &crm2, &kincm, &kinca,
 &jknflg, &jgcflg, &kincv, &jvflg, res);

 MLPutRealList(stdlink, res, 2);

 MLEndPacket(stdlink);
}

extern double f90toppiklist_(double* energy, double* tm, double* tg,
double* alphas, double* scale, double* cutn, double* cutv, double* c0,
double* c1, double* c2, double* cdeltc, double* cdeltl, double* cfullc,
double* cfulll, double* crm2, double* kincm, double* kinca, int* jknflg,
int* jgcflg, double* kincv, int* jvflg, double* res, double* list);

static void ttbarlist(double energy,double tm, double tg, double alphas,
double scale, double cutn, double cutv, double c0, double c1, double c2,
double cdeltc, double cdeltl, double cfullc, double cfulll, double crm2,
double kincm, double kinca, int jknflg, int jgcflg, double kincv, int jvflg){
double res[2]; double list[5*360];

 f90toppiklist_(&energy, &tm, &tg, &alphas, &scale, &cutn, &cutv, &c0,
 &c1, &c2, &cdeltc, &cdeltl, &cfullc, &cfulll, &crm2, &kincm, &kinca,
 &jknflg, &jgcflg, &kincv, &jvflg, res, list);

   MLPutFunction(stdlink, "List", 2 ); MLPutRealList(stdlink, res, 2);
   MLPutFunction(stdlink, "Partition", 2 ); MLPutRealList(stdlink, list, 5*360);
   MLPutInteger(stdlink, 5);

 MLEndPacket(stdlink);
}

extern double f90vssll_(int* nl, double* ah, double* as, double* result);

static double vssll(int nl, double ah, double as){
  double res;

   f90vssll_(&nl, &ah, &as, &res);

   return res;
}

extern double f90mnllplusnnllnonmixc1_(int* nl, double* ah, double* as, double* au, double* result);

static double mnllplusnnllnonmixc1(int nl, double ah, double as, double au){
  double res;

   f90mnllplusnnllnonmixc1_(&nl, &ah, &as, &au, &res);

   return res;
}

extern double f90mnnllallc1inclsoftmixlog_(int* nl, double* ah, double* as, double* au,
double* nu, double* hh, double* ss, double* result);

static double mnnllallc1inclsoftmixlog(int nl, double ah, double as, double au,
double nu, double hh, double ss){
  double res;

   f90mnnllallc1inclsoftmixlog_(&nl, &ah, &as, &au, &nu, &hh, &ss, &res);

   return res;
}

extern double f90xinnllnonmix_(int* nl, double* ah, double* as, double* au,
double* hh, double* ss, double* result);

static double xinnllnonmix(int nl, double ah, double as, double au,
double hh, double ss){
  double res;

   f90xinnllnonmix_(&nl, &ah, &as, &au, &hh, &ss, &res);

   return res;
}

extern double f90xinnllsoftmixlogc1_(double* ah, double* nu, double* hh, double* result);

static double xinnllsoftmixlogc1(double ah, double nu, double hh){
  double res;

   f90xinnllsoftmixlogc1_(&ah, &nu, &hh, &res);

   return res;
}

extern double f90vk1sll_(int* nl, double* ah, double* as, double* result);

static double vk1sll(int nl, double ah, double as){
  double res;

   f90vk1sll_(&nl, &ah, &as, &res);

   return res;
}

extern double f90vk2sll_(int* nl, double* ah, double* as, double* result);

static double vk2sll(int nl, double ah, double as){
  double res;

   f90vk2sll_(&nl, &ah, &as, &res);

   return res;
}

extern double f90vkeffsll_(int* nl, double* ah, double* as, double* result);

static double vkeffsll(int nl, double ah, double as){
  double res;

   f90vkeffsll_(&nl, &ah, &as, &res);

   return res;
}

extern double f90vrsll_(int* nl, double* ah, double* as, double* result);

static double vrsll(int nl, double ah, double as){
  double res;

   f90vrsll_(&nl, &ah, &as, &res);

   return res;
}

extern double f90v2sll_(int* nl, double* ah, double* au, double* as, double* result);

static double v2sll(int nl, double ah, double au, double as){
  double res;

   f90v2sll_(&nl, &ah, &au, &as, &res);

   return res;
}

extern double f90xinll_(int* nl, double* ah, double* au, double* as, double* result);

static double xinll(int nl, double ah, double au, double as){
  double res;

   f90xinll_(&nl, &ah, &au, &as, &res);

   return res;
}

extern double f90mnllc1_(int* nl, double* ah, double* au, double* as, double* result);

static double mnllc1(int nl, double ah, double au, double as){
  double res;

   f90mnllc1_(&nl, &ah, &au, &as, &res);

   return res;
}

extern double f90xinnllmixusoft_(int* nl, double* ah, double* as, double* result);

static double xinnllmixusoft(int nl, double ah, double as){
  double res;

   f90xinnllmixusoft_(&nl, &ah, &as, &res);

   return res;
}

extern double f90mllc2_(int* nl, double* ah, double* as, double* result);

static double mllc2(int nl, double ah, double as){
  double res;

   f90mllc2_(&nl, &ah, &as, &res);

   return res;
}

extern double f90vceffsnnll_(int* nl, double* ah, double* au, double* as, double* result);

static double vceffsnnll(int nl, double ah, double au, double as){
  double res;

   f90vceffsnnll_(&nl, &ah, &au, &as, &res);

   return res;
}

extern double f90a1pole_(int* nl, int* order, double* En, double* mtpole,
double* gamtop, double* asoft, double* VcsNNLL, double* musoft, double* result);

static double a1pole(int nl, int order, double En, double mtpole,
double gamtop, double asoft, double VcsNNLL, double musoft){
  double res;

   f90a1pole_(&nl, &order, &En, &mtpole, &gamtop, &asoft, &VcsNNLL, &musoft, &res);

   return res;
}

extern double f90rnrqcd_(int* nl, int* order, char const* scheme,
char const* method, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass,
int* ord1S, double* R1S, double* muLam, double* xLam, double* mZ, double* aMz,
double* Q, double* mt, double* gt, double* h, double* nu, double* res);

static double rnrqcd(int nl, int order, char const* scheme, char const* method,
int orderAlpha, int runAlpha, int orderMass, int runMass, int ord1S, double R1S,
double muLam, double xLam, double mZ, double aMz, double Q, double mt,
double gt, double h, double nu){
  double res;

   f90rnrqcd_(&nl, &order, scheme, method, &orderAlpha, &runAlpha, &orderMass,
   &runMass, &ord1S, &R1S, &muLam, &xLam, &mZ, &aMz, &Q, &mt, &gt, &h, &nu,
   &res);

   return res;
}

extern double f90qswitch_(int* nl, int* orderAlpha, int* runAlpha, int* orderMass,
int* runMass, int* ord1S, double* muLam, double* xLam, char const* method, double* mZ,
double* aMz, double* mt, double*gt, double* R, double* res);

static double qswitch(int nl, int orderAlpha, int runAlpha, int orderMass,
int runMass, int ord1S, double muLam, double xLam, char const* method, double mZ,
double aMz, double mt, double gt, double R){
  double res;

   f90qswitch_(&nl, &orderAlpha, &runAlpha, &orderMass, &runMass, &ord1S, &muLam,
   &xLam, method, &mZ, &aMz, &mt, &gt, &R, &res);

   return res;

}

extern double f90delta1s_(int* nl, int* orderAlpha, int* runAlpha, int* orderMass,
int* runMass, double* muLam, double* xLam, char const* method, double* mZ,
double* aMz, double* mt, double* R, double* res);

static void delta1s(int nl, int orderAlpha, int runAlpha, int orderMass,
int runMass, double muLam, double xLam, char const* method, double mZ,
double aMz, double mt, double R){
  double res[5];

   f90delta1s_(&nl, &orderAlpha, &runAlpha, &orderMass, &runMass, &muLam,
   &xLam, method, &mZ, &aMz, &mt, &R, res);

  MLPutRealList(stdlink, res, 5);

}

extern double f90vcsll_(double* as, double* result);

static double vcsll(double as){
  double res;

   f90vcsll_(&as, &res);

   return res;
}

extern double f90qfromv_(double* v, double* m, double* gt, double* result);

static double qfromv(double v, double m, double gt){
  double res;

   f90qfromv_(&v, &m, &gt, &res);

   return res;
}

extern double f90vc_(double* v, double* m, double* gt, double* result);

static void vc(double v, double m, double gt){
  double res[2];

   f90vc_(&v, &m, &gt, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
 }

extern double f90vstar_(double* v, double* m, double* gt, double* result);

static double vstar(double v, double m, double gt){
  double res;

   f90vstar_(&v, &m, &gt, &res);

   return res;
}

extern double f90vrootstar_(double* v, double* m, double* gt, double* result);

static double vrootstar(double v, double m, double gt){
  double res;

   f90vrootstar_(&v, &m, &gt, &res);

   return res;
}

extern double f90switchoff_(double* q, double* m, double* gt, double* v0,
  double* v1, double* result);

static double switchoff(double q, double m, double gt, double v0, double v1){
  double res;

   f90switchoff_(&q, &m, &gt, &v0, &v1, &res);

   return res;
}

extern double f90hypgeo_(double* ar, double* ai, double* br, double* bi,
  double* cr, double* ci, double* zr, double* zi, double* res);

static void hypgeo(double ar, double ai, double br, double bi,
  double cr, double ci, double zr, double zi){
  double res[2];

   f90hypgeo_(&ar, &ai, &br, &bi, &cr, &ci, &zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90cdigamma_(double* zr, double* zi, double* res);

static void cdigamma(double zr, double zi){
  double res[2];

   f90cdigamma_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90ctrigamma_(double* zr, double* zi, double* res);

static void ctrigamma(double zr, double zi){
  double res[2];

   f90ctrigamma_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90ewfactors_(int* nf, double* Q, double* Mz, double* GammaZ, double* sin2ThetaW, double* res);

static void ewfactors(int nf, double Q, double Mz, double GammaZ, double sin2ThetaW){
  double res[2];

   f90ewfactors_(&nf, &Q, &Mz, &GammaZ, &sin2ThetaW, res);

   MLPutRealList(stdlink, res, 2);
   MLEndPacket(stdlink);
}

extern double f90kernels_(int* n, double* width, double* w, double* mu, double* p, double* result);

static void kernels(int n, double width, double w, double mu, double p){
  double res[n+1];

   f90kernels_(&n, &width, &w, &mu, &p, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90gammaderlist_(int* n, double* w, double* result);

static void gammaderlist(int n, double w){
  double res[n+1];

   f90gammaderlist_(&n, &w, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90polygamma_(int* n, double* w, double* result);

static void polygamma(int n, double w){
  double res[n+1];

   f90polygamma_(&n, &w, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90nglkernels_(int* n, int* n1, int*n2, double* width, double* w, double* mu, double* p, double* result);

static void nglkernels(int n, int n1, int n2, double width, double w[], long wlen, double mu, double p[], long plen){
  double res[n];

   f90nglkernels_(&n, &n1, &n2, &width, w, &mu, p, res);

   MLPutRealList(stdlink, res, n);

   MLEndPacket(stdlink);
}

extern double f90nglintegral_(int* nf, int* pow, double* w1, double* w2, double* res);

static double nglintegral(int nf, int pow, double w1, double w2){
  double res;

   f90nglintegral_(&nf, &pow, &w1, &w2, &res);

   return res;
}

extern double f90ngldoubleintegral_(int* nf, int* pow, double* w1, double* w2, double* r, double* res);

static double ngldoubleintegral(int nf, int pow, double w1, double w2, double r){
  double res;

   f90ngldoubleintegral_(&nf, &pow, &w1, &w2, &r, &res);

   return res;
}

extern double f90singularhjmpiece_(char const* hard, char const* gap, char const* space, char const* cum,
 int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf, double* j3,
 double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double singularhjmpiece(char const* hard, char const* gap, char const* space, char const* cum,
 int orderAlpha, int runAlpha, int order, int run, int isoft, int nf, double j3, double s3,
 double s31, double s32, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
 double muJ, double muS, double R, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau){
  double res;

f90singularhjmpiece_(hard, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardoublepiece_(char const* hard, char const* gap, char const* space, char const* cum1,
 char const* cum2, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft,
 int* nf, double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ,
 double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 int* c, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* rho1,
 double* rho2, double* res);

static double singulardoublepiece(char const* hard, char const* gap, char const* space, char const* cum1,
 char const* cum2, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, int c[], long len, double lambda,
 double R0, double mu0, double delta0, double h, double rho1, double rho2){
  double res;

f90singulardoublepiece_(hard, gap, space, cum1, cum2, &orderAlpha, &runAlpha, &order, &run,
 &isoft, &nf, &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda, &Q, &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &rho1,
 &rho2, &res);

return res;

}

extern double f90legendrelist_(int* n, int* k, double* x, double* res);

static void legendrelist(int n, int k, double x){
  double res[n + 1];

   f90legendrelist_(&n, &k, &x, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);

}

extern double f90qlegendrelist_(int* n, double* x, double* res);

static void qlegendrelist(int n, double x){
  double res[n + 1];

   f90qlegendrelist_(&n, &x, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);

}

extern double f90singularhjm_(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singularhjm(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, double c[], long len, int clen,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;

f90singularhjm_(hard, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ,
 &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularhjm1d_(char const* hard, char const* gap,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singularhjm1d(char const* hard, char const* gap,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, double c[], long len,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;

 int clen = len;

f90singularhjm1d_(hard, gap, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularhjm1dpiece_(char const* hard, char const* gap,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* res);

static double singularhjm1dpiece(char const* hard, char const* gap,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, int c[], long len, double lambda,
 double R0, double mu0, double delta0, double h, double tau){
  double res;

f90singularhjm1dpiece_(hard, gap, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardouble_(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum1, char const* cum2, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* isoft, int* nf, double* j3, double* s3, double* s31, double* s32, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* mu, double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0,
 double* h, double* rho1, double* rho2, double* res);

static double singulardouble(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum1, char const* cum2, int orderAlpha, int runAlpha, int order, int run,
 int isoft, int nf, double j3, double s3, double s31, double s32,double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double muH, double muJ, double muS, double R, double mu,
 double c[], long len, int clen, double lambda, double R0, double mu0, double delta0,
 double h, double rho1, double rho2){
  double res;

f90singulardouble_(hard, setup, gap, space, cum1, cum2, &orderAlpha, &runAlpha, &order, &run,
 &isoft, &nf, &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda, &Q, &muH, &muJ, &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h,
 &rho1, &rho2, &res);

return res;

}

extern double f90profiles_(double* Q, double* mu0, double* R0, double* n0, double* n1,
double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH, double* eS,
double* eJ, double* eR, int* ns, double* tau, double* res);

static void profiles(double Q, double mu0, double R0, double n0, double n1, double t2,
double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
double eR, int ns, double tau){
  double result[5];

   f90profiles_(&Q, &mu0, &R0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR,
   &ns, &tau, result);

   MLPutRealList(stdlink, result, 5);

   MLEndPacket(stdlink);
}

extern double f90profilesmass_(double* Q, double* beta, double* mu0, double* delLamb,
double* R0, double* n0, double* delta0, double* n1, double* delta1, double* t2, double* ts,
double* slope, double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
int* ns, char const* def, char const* EShape, double* tau, double* res);

static void profilesmass(double Q, double beta, double mu0, double delLamb, double R0,
double n0, double delta0, double n1, double delta1, double t2, double ts, double slope,
double cnt, double eH, double eS, double eJ, double mass, double muM, int ns,
char const* def, char const* EShape, double tau){
  double result[6];

   f90profilesmass_(&Q, &beta, &mu0, &delLamb, &R0, &n0, &delta0, &n1, &delta1, &t2,
   &ts, &slope, &cnt, &eH, &eS, &eJ, &mass, &muM, &ns, def, EShape, &tau, result);

   MLPutRealList(stdlink, result, 6);

   MLEndPacket(stdlink);
}

extern double f90massintersection_(double* Q, double* beta, double* mu0, double* delLamb,
double* n0, double* delta0, double* n1, double* delta1, double* t2, double* ts,
double* slope, double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
char const* def, char const* EShape, double* res);

static void massintersection(double Q, double beta, double mu0, double delLamb, double n0,
double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
double eH, double eS, double eJ, double mass, double muM, char const* def,
char const* EShape){
  double result[2];

   f90massintersection_(&Q, &beta, &mu0, &delLamb, &n0, &delta0, &n1, &delta1, &t2,
   &ts, &slope, &cnt, &eH, &eS, &eJ, &mass, &muM, def, EShape, result);

   MLPutRealList(stdlink, result, 2);

   MLEndPacket(stdlink);
}

extern double f90gammar_(char const* str, int* nf, double* result);

static void gammar(char const* str, int nf){
  double result[3];

   f90gammar_(str, &nf, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90masslessprof_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* res);

static double masslessprof(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double tau){
  double res;
  int clen = len;

f90masslessprof_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &res);

return res;

}

extern double f90findorigin_(char const* shape, char const* gap, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* mu0, double* Rat0, double* n0, double* n1, double* t2, double* tR,
 double* ts, double* slope, double* cnt, double* eH, double* eS, double* eR,
 double* R0, double* muR0, double* delta0, double* h, double* res);

static double findorigin(char const* shape, char const* gap, int orderAlpha, int runAlpha,
 int order, int run, int nf, double mZ, double aMz, double mT, double muT, double mB,
 double muB, double mC, double muC, double muLambda, double Q, double mu0, double Rat0,
 double n0, double n1, double t2, double tR, double ts, double slope, double cnt,
 double eH, double eS, double eR, double R0, double muR0, double delta0, double h){
  double res;

f90findorigin_(shape, gap, &orderAlpha, &runAlpha, &order, &run, &nf, &mZ, &aMz, &mT,
&muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope,
&cnt, &eH, &eS, &eR, &R0, &muR0, &delta0, &h, &res);

return res;

}

extern double f90masslessprofpiece_(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* res);

static void masslessprofpiece(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90masslessprofpiece_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessprofpiecelist_(char const* terms, char const* hard,
 char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* taulist, int* taulen, double* res);

static void masslessprofpiecelist(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double taulist[], long tlen){
  int taulen = tlen;
  double res[tlen * (clen + 1) * (clen + 2)/2];

f90masslessprofpiecelist_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, tlen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslesspiecebin_(char const* terms, char const* hard,
 char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* taulist, int* taulen, double* res);

static void masslesspiecebin(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double taulist[], long tlen, int taulen){
  double res[taulen * (clen + 1) * (clen + 2)/2];

f90masslesspiecebin_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, taulen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessprofdiffpiece_(char const* terms, char const* hard,
char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* tau, double* tau2, double* res);

static void masslessprofdiffpiece(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double tau, double tau2){
  double res[(clen + 1) * (clen + 2)/2];

f90masslessprofdiffpiece_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessproflist_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* taulist, int* tlen,
 double* result);

static void masslessproflist(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double taulist[], long taulen){
  double result[taulen];
  int clen = len;
  int tlen = taulen;

f90masslessproflist_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &tlen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);
}

extern double f90masslessbinlist_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* taulist, int* tlen,
 double* result);

static void masslessbinlist(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double taulist[], long taulen, int tlen){
  double result[tlen];
  int clen = len;

f90masslessbinlist_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &tlen, result);

   MLPutRealList(stdlink, result, tlen);

   MLEndPacket(stdlink);
}

extern double f90masslessprofdiff_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* tau2,
 double* res);

static double masslessprofdiff(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double tau, double tau2){
  double res;
  int clen = len;

f90masslessprofdiff_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90masslessmoment_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, int* orderAlpha, int* runAlpha,
 int* order, int* run, int* nf, double* j3, double* s3, double* G3, double* mZ,
 double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* mu0, double* Rat0, double* n0, double* n1,
 double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH, double* eS,
 double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* tau, double* tau2, int* pow,
 double* res);

static double masslessmoment(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, int orderAlpha, int runAlpha,
 int order, int run, int nf, double j3, double s3, double G3, double mZ, double aMz,
 double mT, double muT, double mB, double muB, double mC, double muC, double muLambda,
 double Q, double mu0, double Rat0, double n0, double n1, double t2, double tR, double ts,
 double slope, double cnt, double eH, double eS, double eJ, double eR, int ns, double c[],
 long len, double lambda, double R0, double muR0, double delta0, double h, double tau,
 double tau2, int pow){
  double res;
  int clen = len;

f90masslessmoment_(terms, hard, shape, setup, gap, space, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, &pow, &res);

return res;

}

extern double f90singular_(char const* hard, char const* shape, char const* setup, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singular(char const* hard, char const* shape, char const* setup, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double tau){
  double res;
  int clen = len;

f90singular_(hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularlist_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static void singularlist(char const* hard, char const* shape, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, int clen, double lambda,
double R0, double mu0, double delta0, double h, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90singularlist_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);

}


extern double f90singulardiff_(char const* hard, char const* shape, char const* setup, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau1, double* tau2, double* res);

static double singulardiff(char const* hard, char const* shape, char const* setup, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double tau1, double tau2){
  double res;
  int clen = len;

f90singulardiff_(hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau1, &tau2, &res);

return res;

}

extern double f90massiveprof_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* res);

static double massiveprof(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90massiveprof_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &res);

return res;

}

extern double f90massorigin_(char const* shape, char const* Eshape, char const* gap,
 char const* scheme, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass,
 int* order, int* run, int* nf, double* mZ, double* aMz, double* mT, double* muT,
 double* mB, double* muB, double* mC, double* muC, double* muLambda1, double* muLambda2,
 double* Q, double* beta, double* mu0, double* deltaLambda, double* Rat0, double* n0,
 double* delta0, double* n1, double* delta1, double* t2, double* ts, double* slope,
 double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
 double* R0, double* muR0, double* del0, double* h, double* res);

static double massorigin(char const* shape, char const* Eshape, char const* gap,
 char const* scheme, int orderAlpha, int runAlpha, int orderMass, int runMass, int order,
 int run, int nf, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda1, double muLambda2, double Q, double beta,
 double mu0, double deltaLambda, double Rat0, double n0, double delta0, double n1,
 double delta1, double t2, double ts, double slope, double cnt, double eH, double eS,
 double eJ, double mass, double muM, double R0, double muR0, double del0, double h){
  double res;

f90massorigin_(shape, Eshape, gap, scheme, &orderAlpha, &runAlpha, &orderMass, &runMass,
 &order, &run, &nf, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2,
 &Q, &beta, &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt,
 &eH, &eS, &eJ, &mass, &muM, &R0, &muR0, &del0, &h, &res);

return res;

}

extern double f90massiveprofpiece_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* res);

static void massiveprofpiece(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90massiveprofpiece_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofpiecelist_(char const* terms, char const* hard,
 char const* shape, char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, int* taulen, double* res);

static void massiveprofpiecelist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau[], long tlen){
  int taulen = tlen;
  double res[tlen * (clen + 1) * (clen + 2)/2];

f90massiveprofpiecelist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tau, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, tlen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massivepiecebin_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, int* taulen, double* res);

static void massivepiecebin(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau[], long tlen, int taulen){
  double res[taulen * (clen + 1) * (clen + 2)/2];

f90massivepiecebin_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tau, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, taulen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofdiffpiece_(char const* terms, char const* hard,
 char const* shape, char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, double* res);

static void massiveprofdiffpiece(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2){
  double res[(clen + 1) * (clen + 2)/2];

f90massiveprofdiffpiece_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofdiff_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, double* res);

static double massiveprofdiff(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2){
  double res;
  int clen = len;

f90massiveprofdiff_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, &res);

return res;

}

extern double f90massiveproflist_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tauList, int* taulen, double* res);

static void massiveproflist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tauList[], long ctau){
  int clen = len;
  int taulen = ctau;
  double result[taulen];

f90massiveproflist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tauList, &taulen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);

}

extern double f90massivebinlist_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tauList, int* taulen, double* res);

static void massivebinlist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tauList[], long ctau, int taulen){
  int clen = len;
  double result[taulen];

f90massivebinlist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tauList, &taulen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);

}

// extern double f90massivediffprof_(char const* terms, char const* hard, char const* shape,
//  char const* Eshape, char const* setup, char const* gap, char const* space,
//  char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
//  double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
//  int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
//  double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
//  double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
//  double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
//  double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
//  double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
//  double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
//  double* tau, double* tau2, double* res);
//
// static double massivediffprof(char const* terms, char const* hard, char const* shape,
//  char const* Eshape, char const* setup, char const* gap, char const* space,
//  char const* cum, char const* scheme, char const* abs, char const* current, double xi,
//  double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
//  int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
//  double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
//  double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
//  double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
//  double eH, double eS, double eJ, double mass, double muM, int ns, double width,
//  double c[], long len, double lambda, double R0, double muR0, double del0, double h,
//  double gammaZ, double sinW, double tau, double tau2){
//   double res;
//   int clen = len;
//
// f90massivediffprof_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs,
//  current, &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3,
//  &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
//  &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
//  &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
//  &tau, &tau2, &res);
//
// return res;
//
// }

extern double f90massivemoment_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, int* pow, double* res);

static double massivemoment(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2, int pow){
  double res;
  int clen = len;

f90massivemoment_(terms, hard, shape, Eshape, setup, gap, space, scheme, abs,
 current, &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3,
 &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, &pow, &res);

return res;

}

extern double f90singularmass_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* res);

static double singularmass(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
char const* abs, char const* current, double xi, double xiB, int orderAlpha, int runAlpha,
int orderMass, int runMass, int order, int run, int nf, double j3, double s3, double G3,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90singularmass_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW,
&tau, &res);

return res;

}

extern double f90singularmassdiff_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* tau2, double* res);

static double singularmassdiff(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
char const* abs, char const* current, double xi, double xiB, int orderAlpha, int runAlpha,
int orderMass, int runMass, int order, int run, int nf, double j3, double s3, double G3,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;
  int clen = len;

f90singularmassdiff_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW,
&tau, &tau2, &res);

return res;

}

extern double f90massnondist_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf,
 double* G3, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* Rmass, double* muM, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double massnondist(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run, int nf,
 double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda1, double muLambda2, double Q, double muH, double muJ,
 double muS, double R, double Rmass, double muM, double mu, double c[], long len,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;
  int clen = len;

f90massnondist_(hard, shape, Eshape, setup, gap, space, cum, scheme, &orderAlpha,
 &runAlpha, &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB,
 &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu,
 c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90massnondistdiff_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf,
 double* G3, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* Rmass, double* muM, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* tau2, double* res);

static double massnondistdiff(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run, int nf,
 double G3, double mZ, double aMz, double mT, double muT, double mB, double muB, double mC,
 double muC, double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS,
 double R, double Rmass, double muM, double mu, double c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau, double tau2){
  double res;
  int clen = len;

f90massnondistdiff_(hard, shape, Eshape, setup, gap, space, cum, scheme, &orderAlpha,
&runAlpha, &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB,
&muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu,
c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90massnondistpiece_(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int* orderAlpha,
 int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH, double* muJ,
 double* muS, double* R, double* Rmass, double* muM, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* res);

static double massnondistpiece(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int orderAlpha,
 int runAlpha, int orderMass, int runMass, int order, int run, int nf, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
 double Rmass, double muM, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau){
  double res;

f90massnondistpiece_(hard, shape, Eshape, gap, space, cum, scheme, &orderAlpha, &runAlpha,
 &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu, c, &lambda, &R0,
 &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90massnondistdiffpiece_(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int* orderAlpha,
 int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH, double* muJ,
 double* muS, double* R, double* Rmass, double* muM, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* tau2,
 double* res);

static double massnondistdiffpiece(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int orderAlpha,
 int runAlpha, int orderMass, int runMass, int order, int run, int nf, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
 double Rmass, double muM, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau, double tau2){
  double res;

f90massnondistdiffpiece_(hard, shape, Eshape, gap, space, cum, scheme, &orderAlpha, &runAlpha,
&orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
&muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu, c, &lambda,
&R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90singularmasspiece_(char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* res);

static double singularmasspiece(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme, char const* abs,
char const* current, double xi, double xiB, int orderAlpha, int runAlpha, int orderMass,
int runMass, int order, int run, int nf, double j3, double s3, double G3, double mZ,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, int c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90singularmasspiece_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90singularmassdiffpiece_(char const* hard, char const* shape,
 char const* setup, char const* Eshape, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* tau2, double* res);

static double singularmassdiffpiece(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme, char const* abs,
char const* current, double xi, double xiB, int orderAlpha, int runAlpha, int orderMass,
int runMass, int order, int run, int nf, double j3, double s3, double G3, double mZ,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, int c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;

f90singularmassdiffpiece_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2,
&res);

return res;

}

extern double f90nsmass_(char const* shape, char const* setup, char const* gap,
 char const* cum, char const* scheme, char const* abs, char const* current,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muBottom,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* mu,
 double* muM, double* muB, double* muS, double* R, double* Rmass, double* width,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* gammaZ, double* sinW, double* tau, double* res);

static double nsmass(char const* shape, char const* setup, char const* gap,
char const* cum, char const* scheme, char const* abs, char const* current, int orderAlpha,
int runAlpha, int order, int run, int orderMass, int runMass, int nf, double mZ,
double aMz, double mT, double muT, double mB, double muBottom, double mC, double muC,
double muLambda1, double muLambda2, double Q, double mu, double muM, double muB,
double muS, double R, double Rmass, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90nsmass_(shape, setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90nsmassdiff_(char const* shape, char const* setup, char const* gap,
 char const* cum, char const* scheme, char const* abs, char const* current,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muBottom,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* mu,
 double* muM, double* muB, double* muS, double* R, double* Rmass, double* width,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* gammaZ, double* sinW, double* tau, double* tau2, double* res);

static double nsmassdiff(char const* shape, char const* setup, char const* gap,
char const* cum, char const* scheme, char const* abs, char const* current, int orderAlpha,
int runAlpha, int order, int run, int orderMass, int runMass, int nf, double mZ,
double aMz, double mT, double muT, double mB, double muBottom, double mC, double muC,
double muLambda1, double muLambda2, double Q, double mu, double muM, double muB,
double muS, double R, double Rmass, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;
  int clen = len;

f90nsmassdiff_(shape, setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90hjmnsmass_(char const* setup, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* res);

static double hjmnsmass(char const* setup, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, double c[], long len, int clen, double lambda, double R0,
double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90hjmnsmass_(setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90hjmnsmassdiff_(char const* setup, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* tau2, double* res);

static double hjmnsmassdiff(char const* setup, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, double c[], long len, int clen, double lambda, double R0,
double mu0, double delta0, double h, double gammaZ, double sinW, double tau, double tau2){
  double res;

f90hjmnsmassdiff_(setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90nsmasspiece_(char const* shape, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* res);

static double nsmasspiece(char const* shape, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, int c[], long len, double lambda, double R0, double mu0,
double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90nsmasspiece_(shape, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda1, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &lambda, &R0,
&mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90nsmassdiffpiece_(char const* shape, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* tau2, double* res);

static double nsmassdiffpiece(char const* shape, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, int c[], long len, double lambda, double R0, double mu0,
double delta0, double h, double gammaZ, double sinW, double tau, double tau2){
  double res;

f90nsmassdiffpiece_(shape, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda1, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &lambda, &R0,
&mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90singularpiece_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double singularpiece(char const* hard, char const* shape, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int nf, double j3,
 double s3, double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda, double Q, double muH, double muJ, double muS,
 double R, double mu, int c[], long len, double lambda, double R0, double mu0, double delta0,
 double h, double tau){
  double res;

f90singularpiece_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf, &j3,
             &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH,
             &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardiffpiece_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* tau2, double* res);

static double singulardiffpiece(char const* hard, char const* shape, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int nf, double j3,
 double s3, double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda, double Q, double muH, double muJ, double muS,
 double R, double mu, int c[], long len, double lambda, double R0, double mu0, double delta0,
 double h, double tau, double tau2){
  double res;

f90singulardiffpiece_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ,
&muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90anomdim_(char const* str, int* nf, double* G4, double* result);

static void anomdim(char const* str, int nf, double G4){
  double result[5];

   f90anomdim_(str, &nf, &G4, result);

   MLPutRealList(stdlink, result, 5);
   MLEndPacket(stdlink);
}

extern double f90afroms_(char const* str, int* nf, double* G4, double* result);

static void afroms(char const* str, int nf, double G4){
  double result[4];

   f90afroms_(str, &nf, &G4, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90ql_(char const* str, int* nf, int* nl, double* G4, double* result);

static void ql(char const* str, int nf, int order, double G4){
  double result[4];

   f90ql_(str, &nf, &order, &G4, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90anasy_(char const* str, int* nf, int* nl, double* G4, double* result);

static void anasy(char const* str, int nf, int order, double G4){
  double result[4];

   f90anasy_(str, &nf, &order, &G4, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90anlambda_(char const* str, int* nf, double* G4, double* result);

static void anlambda(char const* str, int nf, double G4){
  double result[4];

   f90anlambda_(str, &nf, &G4, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90msbardeltapiece_(int* nl, int* nh, double* result);

static void msbardeltapiece(int nl, int nh){
  double result[20];

   f90msbardeltapiece_(&nl, &nh, result);

   MLPutFunction(stdlink, "Partition", 2);
   MLPutRealList(stdlink, result, 20);
   MLPutInteger(stdlink, 5);
   MLEndPacket(stdlink);
}

extern double f90alphamatchinglog_(char const* str, char const* direction,
  int* nf, double* result);

static void alphamatchinglog(char const* str, char const* direction, int nf){
  double result[25];

   f90alphamatchinglog_(str, direction, &nf, result);

   MLPutFunction(stdlink, "Partition", 2);
   MLPutRealList(stdlink, result, 25);
   MLPutInteger(stdlink, 5);
   MLEndPacket(stdlink);
}

extern double f90alphamatchinginverse_(char const* str, int* nf, double* result);

static void alphamatchinginverse(char const* str, int nf){
  double result[5];

   f90alphamatchinginverse_(str, &nf, result);

  //  MLPutFunction(stdlink, "Partition", 2);
   MLPutRealList(stdlink, result, 5);
  //  MLPutInteger(stdlink, 5);
   MLEndPacket(stdlink);
}

extern double f90ccoef_(int* nf, int* order, int* n, double* result);

static void ccoef(int nf, int order, int n){
  double result[n+1];

   f90ccoef_(&nf, &order, &n, result);

   MLPutRealList(stdlink, result, n+1);
   MLEndPacket(stdlink);
}

extern double f90pscoef_(int* nf, double* lg, double* result);

static void pscoef(int nf, double lg){
  double result[4];

   f90pscoef_(&nf, &lg, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90n12_(char const* str, int* nf, int* order, double* labmda,
  double* err, double* result);

static double n12(char const* str, int nf, int order, double lambda, double err){
  double result;

   f90n12_(str, &nf, &order, &lambda, &err, &result);
   return result;

}

extern double f90n12ratio_(char const* str, int* nf, int* order, double* labmda,
  double* err, double* result);

static double n12ratio(char const* str, int nf, int order, double lambda, double err){
  double result;

   f90n12ratio_(str, &nf, &order, &lambda, &err, &result);
   return result;

}

extern double f90n12residue_(char const* str, int* nf, int* order, double* labmda,
  double* err, double* result);

static double n12residue(char const* str, int nf, int order, double lambda, double err){
  double result;

   f90n12residue_(str, &nf, &order, &lambda, &err, &result);
   return result;

}

extern double f90n12rs_(char const* str, int* nf, int* n, int* order,
  double* labmda, double* err, double* result);

static double n12rs(char const* str, int nf, int n, int order, double lambda,
double err){
  double result;

   f90n12rs_(str, &nf, &n, &order, &lambda, &err, &result);
   return result;

}

extern double f90n12generic_(double* aCoef, int* nf, int* order, double* labmda,
  double* result);

static double n12generic(double aCoef[], long len, int nf, int order, double lambda){
  double result;

   f90n12generic_(aCoef, &nf, &order, &lambda, &result);
   return result;
}

extern double f90scoef_(char const* str, int* nf, double* result);

static void scoef(char const* str, int nf){
  double result[3];

   f90scoef_(str, &nf, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90scoefgamma_(double* gama, int* n, int* nf, double* result);

static void scoefgamma(double gama[], long n, int nf){
  int len = n;
  double result[len];

   f90scoefgamma_(gama, &len, &nf, result);

   MLPutRealList(stdlink, result, len);

   MLEndPacket(stdlink);
}

extern double f90scoeflambda_(char const* str, int* nf, double* lambda, double* result);

static void scoeflambda(char const* str, int nf, double lambda){
  double result[4];

   f90scoeflambda_(str, &nf, &lambda, result);

   MLPutRealList(stdlink, result, 4);

   MLEndPacket(stdlink);
}

extern double f90fomass_(char const* shape, char const* current, double* m, double* Q,
double* Mz, double* gammaZ, double* sin2ThetaW, double* tau, double* result);

static double fomass(char const* shape, char const* current, double m, double Q,
double Mz, double gammaZ, double sin2ThetaW, double tau){
  double result;

  f90fomass_(shape, current, &m, &Q, &Mz, &gammaZ, &sin2ThetaW, &tau, &result);
   return result;
}

extern double f90thrustns1loop_(double* tau, double* result);

static void thrustns1loop(double tau){
  double result[3];

   f90thrustns1loop_(&tau, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90thrustns2loop_(double *er, double* tau, double* result);

static void thrustns2loop(double er, double tau){
  double result[2];

   f90thrustns2loop_(&er, &tau, result);

   MLPutRealList(stdlink, result, 2);

   MLEndPacket(stdlink);
}

extern double f90polylog_(int* nf, double* z, double* result);

static double polylog(int n, double z){
  double res;

   f90polylog_(&n, &z, &res);

   return res;
}

extern double f90dilog_(double* z, double* result);

static double dilog(double z){
  double res;

   f90dilog_(&z, &res);

   return res;
}

extern double f90pfq_(double* a, int* clena, double* b, int* clenb, double* z, double* res);

static double pfq(double a[], long lena, double b[], long lenb, double z){
  int clena = lena; int clenb = lenb;
  double res;

   f90pfq_(a, &clena, b, &clenb, &z, &res);

   return res;
}

extern double f90elliptic3_(double* psi, double* k, double* c, double* result);

static double elliptic3(double psi, double k, double c){
  double res;

   f90elliptic3_(&psi, &k, &c, &res);

   return res;
}

extern double f90nglfunction_(int* nf, double* z, double* result);

static double nglfunction(int n, double z){
  double res;

   f90nglfunction_(&n, &z, &res);

   return res;
}

extern double f90mctop_(char const* str, double* mt, double* Q, int* n, int* k,
  double* x, double* result);

static double mctop(char const* str, double mt, double Q, int n, int k, double x){
  double res;

   f90mctop_(str, &mt, &Q, &n, &k, &x, &res);

   return res;
}


extern double f90breitunstable_(char const* str, double* mt, double* Q,
  double * gamma, int* n, int* k, double* x, double* result);

static double breitunstable(char const* str, double mt, double Q, double gamma,
  int n, int k, double x){
  double res;

   f90breitunstable_(str, &mt, &Q, &gamma, &n, &k, &x, &res);

   return res;
}

extern double f90deltamctop_(char const* str, double* mt, double* Q, double* result);

static double deltamctop(char const* str, double mt, double Q){
  double res;

   f90deltamctop_(str, &mt, &Q, &res);

   return res;
}

extern double f90complexpolylog_(int* nf, double* z1, double* z2, double* result);

static void complexpolylog(int n, double z1, double z2){
  double res[2];

   f90complexpolylog_(&n, &z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90nglsoft_(int* nf, double* z1, double* z2, double* result);

static void nglsoft(int n, double z1, double z2){
  double res[2];

   f90nglsoft_(&n, &z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90deltacharm3_(int* nl, int* nh, double* z, double* result);

static double deltacharm3(int nl, int nh, double z){
  double res;

   f90deltacharm3_(&nl, &nh, &z, &res);
   return res;
}

extern double f90deltacharm3der_(int* nl, int* nh, double* z, double* result);

static double deltacharm3der(int nl, int nh, double z){
  double res;

   f90deltacharm3der_(&nl, &nh, &z, &res);
   return res;
}

extern double f90gammarcharm3_(int* nl, int* nh, double* z, double* result);

static double gammarcharm3(int nl, int nh, double z){
  double res;

   f90gammarcharm3_(&nl, &nh, &z, &res);
   return res;
}

extern double f90deltacharm2_(double* z, double* result);

static double deltacharm2(double z){
  double res;

   f90deltacharm2_(&z, &res);
   return res;
}

extern double f90p2_(double* z, double* result);

static double p2(double z){
  double res;

   f90p2_(&z, &res);
   return res;
}

extern double f90pi0_(double* zr, double* zi, double* result);

static void pi0(double zr, double zi){
  double res[2];

   f90pi0_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi0der_(int* i, double* zr, double* zi, double* result);

static void pi0der(int i, double zr, double zi){
  double res[2];

   f90pi0der_(&i, &zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi1der_(int* i, double* zr, double* zi, double* result);

static void pi1der(int i, double zr, double zi){
  double res[2];

   f90pi1der_(&i, &zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi1_(double* zr, double* zi, double* result);

static void pi1(double zr, double zi){
  double res[2];

   f90pi1_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi3_(double* zr, double* zi, double* result);

static void pi3(double zr, double zi){
  double res[2];

   f90pi3_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi2_(double* zr, double* zi, double* result);

static void pi2(double zr, double zi){
  double res[2];

   f90pi2_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90pi2der_(double* zr, double* zi, double* result);

static void pi2der(double zr, double zi){
  double res[2];

   f90pi2der_(&zr, &zi, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]); MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90p2int_(double* z, double* result);

static double p2int(double z){
  double res;

   f90p2int_(&z, &res);
   return res;
}

extern double f90p2double_(double* z1, double* z2, double* result);

static double p2double(double z1, double z2){
  double res;

   f90p2double_(&z1, &z2, &res);
   return res;
}

extern double f90deltabottomcharm_(double* z1, double* z2, double* result);

static double deltabottomcharm(double z1, double z2){
  double res;

   f90deltabottomcharm_(&z1, &z2, &res);
   return res;
}

extern double f90gammarbottomcharm_(double* z1, double* z2, double* result);

static double gammarbottomcharm(double z1, double z2){
  double res;

   f90gammarbottomcharm_(&z1, &z2, &res);
   return res;
}

extern double f90deltacharmnh_(double* z, double* result);

static double deltacharmnh(double z){
  double res;

   f90deltacharmnh_(&z, &res);
   return res;
}

extern double f90deltacharmglue_(double* z, double* result);

static double deltacharmglue(double z){
  double res;

   f90deltacharmglue_(&z, &res);
   return res;
}

extern double f90deltacharmglueder_(double* z, double* result);

static double deltacharmglueder(double z){
  double res;

   f90deltacharmglueder_(&z, &res);
   return res;
}

extern double f90deltacharmnl_(double* z, double* result);

static double deltacharmnl(double z){
  double res;

   f90deltacharmnl_(&z, &res);
   return res;
}

extern double f90deltacharmnhder_(double* z, double* result);

static double deltacharmnhder(double z){
  double res;

   f90deltacharmnhder_(&z, &res);
   return res;
}

extern double f90deltacharmnlder_(double* z, double* result);

static double deltacharmnlder(double z){
  double res;

   f90deltacharmnlder_(&z, &res);
   return res;
}

extern double f90gammarcharm2_(double* z, double* result);

static double gammarcharm2(double z){
  double res;

   f90gammarcharm2_(&z, &res);
   return res;
}

extern double f90deltacharm2der_(double* z, double* result);

static double deltacharm2der(double z){
  double res;

   f90deltacharm2der_(&z, &res);
   return res;
}

extern double f90upsilondeltacharm_(int* n, int* l, double* alpha, double* mb,
  double* mc, double* result);

static double upsilondeltacharm(int n, int l, double alpha, double mb, double mc){
  double res;

   f90upsilondeltacharm_(&n, &l, &alpha, &mb, &mc, &res);
   return res;
}

extern double f90upsilondeltacharmbin_(int* n, int* l, double* alpha, double* mb,
  double* mc, double* result);

static double upsilondeltacharmbin(int n, int l, double alpha, double mb, double mc){
  double res;

   f90upsilondeltacharmbin_(&n, &l, &alpha, &mb, &mc, &res);
   return res;
}

extern double f90deltacharmexact_(char const* charm, char const* type,
char const* scheme, char const* average, int* n, int* l, int* j, int* s,
int* nl, double* mH, double* mL, double* mu, double* alp, double* result);

static double deltacharmexact(char const* charm, char const* type,
char const* scheme, char const* average, int n, int l, int j, int s, int nl,
double mH, double mL, double mu, double alp){
  double res;

   f90deltacharmexact_(charm, type, scheme, average, &n, &l, &j, &s, &nl, &mH,
   &mL, &mu, &alp, &res);
   return res;
}

extern double f90cli2_(double* z1, double* z2, double* result);

static void cli2(double z1, double z2){
  double res[2];

   f90cli2_(&z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90cli3_(double* z1, double* z2, double* result);

static void cli3(double z1, double z2){
  double res[2];

   f90cli3_(&z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90deltagap_(char const* str, int* orderAlpha, int* runAlpha, int* runMass,
int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double* R, double* result);

static void deltagap(char const* str, int orderAlpha, int runAlpha, int runMass, int nf,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double R){
  double result[4];

   f90deltagap_(str, &orderAlpha, &runAlpha, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &mu, &R, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90psdelta_(int* orderAlpha, int* runAlpha, int* nf, double* mZ,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double* R, double* lg, double* result);

static void psdelta(int orderAlpha, int runAlpha, int nf, double mZ, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu,
double R, double lg){
  double result[4];

   f90psdelta_(&orderAlpha, &runAlpha, &nf, &mZ, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &mu, &R, &lg, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90delta_(char const* str, int* nf, double* mu, double* R, double* result);

static void delta(char const* str, int nf, double mu, double R){
  double result[4];

   f90delta_(str, &nf, &mu, &R, result);

   MLPutRealList(stdlink, result, 4);
   MLEndPacket(stdlink);
}

extern double f90diffdeltagap_(char const* str, char const* scheme, int* order, double*R0,
double* R1, double* mu0, double* mu1, double* muLambda, int* orderAlpha, int* runAlpha,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* result);

static double diffdeltagap(char const* str, char const* scheme, int order, double R0,
double R1, double mu0, double mu1, double muLambda, int orderAlpha, int runAlpha, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC){

  double result;

  f90diffdeltagap_(str, scheme, &order, &R0, &R1, &mu0, &mu1, &muLambda, &orderAlpha,
   &runAlpha, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &result);

   return result;
}

extern double f90diffdeltagapmass_(char const* str, int* order, double* R0, double* R1,
double* mu0, double* mu1, double* muM, double* muLambda1, double* muLambda2,
int* orderAlpha, int* runAlpha, int* runMass, int* nf, double* Mz, double* aMz, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* result);

static double diffdeltagapmass(char const* str, int order, double R0, double R1,
double mu0, double mu1, double muM, double muLambda1, double muLambda2, int orderAlpha,
int runAlpha, int runMass, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC){

  double result;

  f90diffdeltagapmass_(str, &order, &R0, &R1, &mu0, &mu1, &muM, &muLambda1, &muLambda2,
  &orderAlpha, &runAlpha, &runMass, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &result);

   return result;
}

extern double f90coefmat_(char const* str, int* nf, double* s3, double* result);

static void coefmat(char const* str, int nf, double s3){
  double result[15];

   f90coefmat_(str, &nf, &s3, result);

   MLPutFunction(stdlink, "Partition", 2);
   MLPutRealList(stdlink, result, 15);
   MLPutInteger(stdlink, 5);

   MLEndPacket(stdlink);
}

extern double f90model_(double* c, int* clen, double* lambda, int* k, double* l,
  double* result);

static double model(double c[], long clen, double lambda, int k, double l){
  double res;
  int len = clen;

   f90model_(c, &len, &lambda, &k, &l, &res);

   return res;
}

extern double f90modelunstable_(char const* shape, double* mt, double* Q,
double* c, int* clen, double* lambda, int* n, int* k, double* l, double* result);

static double modelunstable(char const* shape, double mt, double Q, double c[], long clen,
double lambda, int n, int k, double l){
  double res;  int len = clen;

   f90modelunstable_(shape, &mt, &Q, c, &len, &lambda, &n, &k, &l, &res);

   return res;
}

extern double f90breitmodelunstable_(char const* shape, double* mt, double* Q,
double *gamma, double* c, int* clen, double* lambda, int* n, int* k, double* l,
double* result);

static double breitmodelunstable(char const* shape, double mt, double Q, double gamma,
  double c[], long clen, double lambda, int n, int k, double l){
  double res;  int len = clen;

   f90breitmodelunstable_(shape, &mt, &Q, &gamma, c, &len, &lambda, &n, &k, &l, &res);

   return res;
}

extern double f90modelunstablediff_(char const* shape, double* mt, double* Q,
double* c, int* clen, double* lambda, int* n, int* k, double* l, double* l2, double* result);

static double modelunstablediff(char const* shape, double mt, double Q, double c[], long clen,
double lambda, int n, int k, double l, double l2){
  double res;  int len = clen;

   f90modelunstablediff_(shape, &mt, &Q, c, &len, &lambda, &n, &k, &l, &l2, &res);

   return res;
}

extern double f90modeldiff_(double* c, int* clen, double* lambda, int* k, double* l,
double* l2, double* result);

static double modeldiff(double c[], long clen, double lambda, int k, double l, double l2){
  double res;
  int len = clen;

   f90modeldiff_(c, &len, &lambda, &k, &l, &l2, &res);

   return res;
}

extern double f90breitmodel_(double* c, int* clen, double* lambda, double* width,
int* k, double* l, double* result);

static double breitmodel(double c[], long clen, double lambda, double width, int k, double l){
  double res;
  int len = clen;

   f90breitmodel_(c, &len, &lambda, &width, &k, &l, &res);

   return res;
}

extern double f90breitmodeldiff_(double* c, int* clen, double* lambda, double* width,
int* k, double* l, double* l2, double* result);

static double breitmodeldiff(double c[], long clen, double lambda, double width,
int k, double l, double l2){
  double res;
  int len = clen;

   f90breitmodeldiff_(c, &len, &lambda, &width, &k, &l, &l2, &res);

   return res;
}

extern double f90taylor_(double* c, int* clen, double* lambda, int* k, double* result);

static double taylor(double c[], long clen, double lambda, int k){
  double res;
  int len = clen;

   f90taylor_(c, &len, &lambda, &k, &res);

   return res;
}

extern double f90momentmodel_(double* c, int* clen, double* lambda, int* k,
double* result);

static double momentmodel(double c[], long clen, double lambda, int k){
  double res;
  int len = clen;

   f90momentmodel_(c, &len, &lambda, &k, &res);

   return res;
}

extern double f90modelpiece_(int* c, double* lambda, int* k, double* l, double* result);

static double modelpiece(int c[], long clen, double lambda, int k, double l){
  double res;

   f90modelpiece_(c, &lambda, &k, &l, &res);

   return res;
}

extern double f90taylorpiece_(int* c, double* lambda, int* k, double* result);

static double taylorpiece(int c[], long clen, double lambda, int k){
  double res;

   f90taylorpiece_(c, &lambda, &k, &res);

   return res;
}

extern double f90wtilde_(int* order, int* nf, double* gamma, double* a0, double* a1,
double* res);

static double wtilde(int order, int nf, double gamma[], long ngamma, double a0, double a1){
  double res;

   f90wtilde_(&order, &nf, gamma, &a0, &a1, &res);

   return res;
}

extern double f90ktilde_(int* order, int* nf, double* gamma, double* a0, double* a1,
double* res);

static double ktilde(int order, int nf, double gamma[], long ngamma, double a0,
double a1){
  double res;

   f90ktilde_(&order, &nf, gamma, &a0, &a1, &res);

   return res;
}

extern double f90alphaqcd_(char const* str, char const* method, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double alphaqcd(char const* str, char const* method, int order, int run, int nf, double Mz, double
aMz, double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90alphaqcd_(str, method, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
   &mu, &res);

  return res;
}

extern double f90alphaqed_(int* nf, double* Mz, double* aMz, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* mu,
double* res);

static double alphaqed(int nf, double Mz, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double mu){

   double res;

   f90alphaqed_(&nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &res);

  return res;
}

extern double f90alphacomplex_(char const* str, char const* method, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* muR, double* muI, double* res);

static void alphacomplex(char const* str, char const* method, int order, int run,
int nf, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double muR, double muI){

   double res[2];

   f90alphacomplex_(str, method, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &muR, &muI, res);

   MLPutFunction(stdlink,"Complex",2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);
   MLEndPacket(stdlink);
}

extern double f90lambdaqcd_(char const* str, int* order, int* runAlpha, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double* res);

static double lambdaqcd(char const* str, int order, int runAlpha, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu){

   double res;

   f90lambdaqcd_(str, &order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
   &muC, &mu, &res);

  return res;
}

extern double f90msbarmass_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double msbarmass(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90msbarmass_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
                 &muC, &mu, &res);

  return res;
}

extern double f90polemass_(int* orderAlpha, int* runAlpha, int*order, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double* res);

static double polemass(int orderAlpha, int runAlpha, int order, int run, int nf, double Mz,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double mu){

   double res;

   f90polemass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
                 &muB, &mC, &muC, &mu, &res);

  return res;
}

extern double f90msbarmasslow_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double msbarmasslow(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90msbarmasslow_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &mu, &res);

  return res;
}

extern double f90msrmass_(char const* type, char const* str, int* orderAlpha,
int* runAlpha, int* order, int* run, int* nf, double* Mz, double* aMz, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* lambda,
double* mu, double* R, double* res);

static double msrmass(char const* type, char const* str, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda,
double mu, double R){

  double res;

  f90msrmass_(type, str, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz,
  &mT, &muT, &mB, &muB, &mC, &muC, &lambda, &mu, &R, &res);

  return res;
}

extern double f90rsmass_(char const* type, char const* method, int* orderAlpha, int* runAlpha,
int* order, int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT,
double* mB, double* muB, double* mC, double* muC, double* lambda, double* R, double* res);

static double rsmass(char const* type, char const* method, int orderAlpha, int runAlpha, int order,
int run, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC, double lambda, double R){

  double res;

  f90rsmass_(type, method, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz,
  &mT, &muT, &mB, &muB, &mC, &muC, &lambda, &R, &res);

  return res;
}

extern double f90msrvfns_(char const* up, char const* type, char const* str,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda, double* mu1, double* mu2, double* R, double* res);

static double msrvfns(char const* up, char const* type, char const* str,
int orderAlpha, int runAlpha, int order, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC,
double lambda, double mu1, double mu2, double R){

  double res;

  f90msrvfns_(up, type, str, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz,
  &aMz, &mT, &muT, &mB, &muB, &mC, &muC,&lambda, &mu1, &mu2, &R, &res);

  return res;
}

extern double f90msrtop_(char const* up, char const* type, char const* str,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda, double* mu1, double* mu2, double* mu3,
double* R, double* res);

static double msrtop(char const* up, char const* type, char const* str,
int orderAlpha, int runAlpha, int order, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC,
double lambda, double mu1, double mu2, double mu3, double R){

  double res;

  f90msrtop_(up, type, str, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz,
  &aMz, &mT, &muT, &mB, &muB, &mC, &muC,&lambda, &mu1, &mu2, &mu3, &R, &res);

  return res;
}

extern double f90findmass_(int* ord, int* n, int* l, int* j, int* s, char const* iter,
char const* charm, char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mass, double* lambda1, double* lambda2, double* lam,
double* mu, double* R, double* res);

static double findmass(int ord, int n, int l, int j, int s, char const* iter,
char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double mass,
double lambda1, double lambda2, double lam, double mu, double R){

  double res;

  f90findmass_(&ord, &n, &l, &j, &s, iter, charm, str, average, method,
  counting, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT,
  &mB, &muB, &mC, &muC, &mass, &lambda1, &lambda2, &lam, &mu, &R, &res);

  return res;

}

extern double f90findenergy_(int* ord, int* n, int* l, int* j, int* s, char const* iter,
char const* charm, char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam,
double* mu, double* R, double* res);

static double findenergy(int ord, int n, int l, int j, int s, char const* iter,
char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC,
double lambda1, double lambda2, double lam, double mu, double R){

  double res;

  f90findenergy_(&ord, &n, &l, &j, &s, iter, charm, str, average, method,
  counting, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT,
  &mB, &muB, &mC, &muC, &lambda1, &lambda2, &lam, &mu, &R, &res);

  return res;

}

extern double f90nrqcderror_(int* n, int* l, int* j, int* s, char const* iter,
char const* charm, char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mass, double* lambda1, double* lambda2, double* lam,
double* mu0, double* mu1, double* deltaMu, double* R0, double* R1,
double* deltaR, double* x, double* res);

static void nrqcderror(int n, int l, int j, int s, char const* iter,
char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double mass,
double lambda1, double lambda2, double lam, double mu0, double mu1,
double deltaMu, double R0, double R1, double deltaR, double x){

  double res[10];

  f90nrqcderror_(&n, &l, &j, &s, iter, charm, str, average, method, counting, &orderAlpha,
  &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &mass, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1, &deltaR,
  &x, res);

 MLPutFunction(stdlink, "Partition", 2 );
 MLPutRealList(stdlink, res, 10);
 MLPutInteger(stdlink, 2);

 MLEndPacket(stdlink);

}

extern double f90masserror_(int* ord, int* n, int* l, int* j, int* s,
char const* iter, char const* charm, char const* str, char const* average,
char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mass, double* lambda1, double* lambda2, double* lam,
double* mu0, double* mu1, double* deltaMu, double* R0, double* R1,
double* deltaR, double* x, double* res);

static void masserror(int ord, int n, int l, int j, int s, char const* iter,
char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double mass,
double lambda1, double lambda2, double lam, double mu0, double mu1,
double deltaMu, double R0, double R1, double deltaR, double x){

  double res[2];

  f90masserror_(&ord, &n, &l, &j, &s, iter, charm, str, average, method, counting, &orderAlpha,
  &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &mass, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1, &deltaR,
  &x, res);

 MLPutRealList(stdlink, res, 2);
 MLEndPacket(stdlink);

}

extern double f90masslist_(int* ord, int* n, int* l, int* j, int* s,
char const* iter, char const* charm, char const* str, char const* average,
char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mass, double* lambda1, double* lambda2, double* lam,
double* mu0, double* mu1, double* deltaMu, double* R0, double* R1,
double* deltaR, double* res);

static void masslist(int ord, int n, int l, int j, int s, char const* iter,
char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double mass,
double lambda1, double lambda2, double lam, double mu0, double mu1,
double deltaMu, double R0, double R1, double deltaR){
  int imax = floor( (mu1 - mu0)/deltaMu ) + 1 ;
  int jmax = floor( (R1 - R0)/deltaR )    + 1 ;

  double res[ 3 * imax * jmax ];

  f90masslist_(&ord, &n, &l, &j, &s, iter, charm, str, average, method, counting, &orderAlpha,
  &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &mass, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1, &deltaR,
  res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, 3 * imax * jmax);
   MLPutInteger(stdlink, 3);
   MLEndPacket(stdlink);

}

extern double f90nrqcdlist_(int* n, int* l, int* j, int* s, char const* iter,
char const* charm, char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mass, double* lambda1, double* lambda2, double* lam,
double* mu0, double* mu1, double* deltaMu, double* R0, double* R1,
double* deltaR, double* res);

static void nrqcdlist(int n, int l, int j, int s, char const* iter, char const* charm,
char const* str, char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double mass,
double lambda1, double lambda2, double lam, double mu0, double mu1,
double deltaMu, double R0, double R1, double deltaR){
  int imax = floor( (mu1 - mu0)/deltaMu ) + 1 ;
  int jmax = floor( (R1 - R0)/deltaR )    + 1 ;

  double res[ 7 * imax * jmax ];

  f90nrqcdlist_(&n, &l, &j, &s, iter, charm, str, average, method, counting, &orderAlpha,
  &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &mass, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1, &deltaR,
  res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, 7 * imax * jmax);
   MLPutInteger(stdlink, 7);
   MLEndPacket(stdlink);

}

extern double f90upsilonlist_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu0,
double* mu1, double* deltaMu, double* R0, double* R1, double* deltaR,
double* epsAlpha, double* epsCharm, double* res);

static void upsilonlist(int n, int l, int j, int s, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda1,
double lambda2, double lam, double mu0, double mu1, double deltaMu, double R0,
double R1, double deltaR, double epsAlpha, double epsCharm){
  int imax = floor( (mu1 - mu0)/deltaMu ) + 1 ;
  int jmax = floor( (R1 - R0)/deltaR )    + 1 ;

  double res[ 15 * imax * jmax ];

  f90upsilonlist_(&n, &l, &j, &s, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1,
  &deltaR, &epsAlpha, &epsCharm, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, 15 * imax * jmax);
   MLPutInteger(stdlink, 5);
   MLPutInteger(stdlink, 3);
   MLEndPacket(stdlink);

}

extern double f90corrmat_(int* qnlist, int* m, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu0,
double* mu1, double* deltaMu, double* R0, double* R1, double* deltaR,
double* epsAlpha, double* epsCharm, double* massList, double* corMat);

static void corrmat(int qnlist[], long len, int m, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda1,
double lambda2, double lam, double mu0, double mu1, double deltaMu, double R0,
double R1, double deltaR, double epsAlpha, double epsCharm){

  double massList[ 25 * m ];
  double corMat[ 5 * m * m ];

  f90corrmat_(qnlist, &m, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1,
  &deltaR, &epsAlpha, &epsCharm, massList, corMat);

   MLPutFunction(stdlink, "List", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, massList, 25 * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, 5);
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, corMat, 5 * m * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, m);
   MLEndPacket(stdlink);

}

extern double f90chi2nrqcd_(int* qnlist, double* datalist, int* m, char const* iter, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* n, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu,
double* R, int* ndim, double* res);

static double chi2nrqcd(int qnlist[], long len1, double datalist[], long len2,
int m, char const* iter, char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
int n, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double lambda1, double lambda2, double lam, double mu[],
long lenMu, double R[], long lenR){

  double res; int ndim = lenMu;

  f90chi2nrqcd_(qnlist, datalist, &m, iter, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &n, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, mu, R, &ndim, &res);

   return res;

}

extern double f90chi2minnrqcd_(int* qnlist, double* datalist, int* m, char const* iter, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* n, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu,
double* R, int* ndim, double* res);

static void chi2minnrqcd(int qnlist[], long len1, double datalist[], long len2,
int m, char const* iter, char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
int n, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double lambda1, double lambda2, double lam, double mu[],
long lenMu, double R[], long lenR){

  double res[3]; int ndim = lenMu;

  f90chi2minnrqcd_(qnlist, datalist, &m, iter, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &n, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, mu, R, &ndim, res);

   MLPutRealList(stdlink, res, 3);

}

extern double f90chi2minalphanrqcd_(int* qnlist, double* datalist, int* m, char const* iter, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* n, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu,
double* R, int* ndim, double* res);

static void chi2minalphanrqcd(int qnlist[], long len1, double datalist[], long len2,
int m, char const* iter, char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
int n, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double lambda1, double lambda2, double lam, double mu[],
long lenMu, double R[], long lenR){

  double res[3]; int ndim = lenMu;

  f90chi2minalphanrqcd_(qnlist, datalist, &m, iter, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &n, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, mu, R, &ndim, res);

   MLPutRealList(stdlink, res, 3);

}

extern double f90chi2minalphambnrqcd_(int* qnlist, double* datalist, int* m, char const* iter, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* n, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu,
double* R, int* ndim, double* res);

static void chi2minalphambnrqcd(int qnlist[], long len1, double datalist[], long len2,
int m, char const* iter, char const* charm, char const* str, char const* average, char const* method,
char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
int n, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double lambda1, double lambda2, double lam, double mu[],
long lenMu, double R[], long lenR){

  double res[6]; int ndim = lenMu;

  f90chi2minalphambnrqcd_(qnlist, datalist, &m, iter, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &n, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, mu, R, &ndim, res);

   MLPutRealList(stdlink, res, 6);

}

extern double f90errmat_(int* qnlist, int* m, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu0,
double* mu1, double* deltaMu, double* R0, double* R1, double* deltaR,
double* epsAlpha, double* epsCharm, double* massList, double* corMat);

static void errmat(int qnlist[], long len, int m, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda1,
double lambda2, double lam, double mu0, double mu1, double deltaMu, double R0,
double R1, double deltaR, double epsAlpha, double epsCharm){

  double massList[ 20 * m ];
  double corMat[ 5 * m * m ];

  f90errmat_(qnlist, &m, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1,
  &deltaR, &epsAlpha, &epsCharm, massList, corMat);

   MLPutFunction(stdlink, "List", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, massList, 20 * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, 4);
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, corMat, 5 * m * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, m);
   MLEndPacket(stdlink);

}

extern double f90errmatrices_(int* qnlist, int* m, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu0,
double* mu1, double* deltaMu, double* R0, double* R1, double* deltaR,
double* epsAlpha, double* epsCharm, double* massList, double* corMat);

static void errmatrices(int qnlist[], long len, int m, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda1,
double lambda2, double lam, double mu0, double mu1, double deltaMu, double R0,
double R1, double deltaR, double epsAlpha, double epsCharm){

  double massList[ 10 * m ];
  double corMat[ 15 * m * m ];

  f90errmatrices_(qnlist, &m, charm, str, average, method, counting,
  &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
  &mC, &muC, &lambda1, &lambda2, &lam, &mu0, &mu1, &deltaMu, &R0, &R1,
  &deltaR, &epsAlpha, &epsCharm, massList, corMat);

   MLPutFunction(stdlink, "List", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, massList, 10 * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, 2);
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, corMat, 15 * m * m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, m);
   MLPutInteger(stdlink, 3);
   MLEndPacket(stdlink);

}

extern double f90nrqcd_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* lambda1, double* lambda2, double* lam, double* mu,
double* R, double* res);

static void nrqcd(int n, int l, int j, int s, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha,
int runAlpha, int order, int run, int nf, double Mz, double aMz, double mT,
double muT, double mB, double muB, double mC, double muC, double lambda1,
double lambda2, double lam, double mu, double R){

  double res[5];

  f90nrqcd_(&n, &l, &j, &s, charm, str, average, method, counting, &orderAlpha,
  &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &lambda1, &lambda2, &lam, &mu, &R, res);

   MLPutRealList(stdlink, res, 5);
   MLEndPacket(stdlink);

}

extern double f90nrqcddercharm_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* lambda1, double* lambda2,
double* lam, double* mu, double* R, double* eps, double* res);

static void nrqcddercharm(int n, int l, int j, int s, char const* charm, char const* str, char const* average,
char const* method, char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double lambda1, double lambda2, double lam, double mu, double R,
double eps){

  double res[5];

  f90nrqcddercharm_(&n, &l, &j, &s, charm, str, average, method, counting, &orderAlpha, &runAlpha, &order,
  &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &lambda1, &lambda2,
  &lam, &mu, &R, &eps, res);

   MLPutRealList(stdlink, res, 5);
   MLEndPacket(stdlink);

}

extern double f90nrqcdderalpha_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* lambda1, double* lambda2,
double* lam, double* mu, double* R, double* eps, double* res);

static void nrqcdderalpha(int n, int l, int j, int s, char const* charm, char const* str, char const* average,
char const* method, char const* counting, int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double lambda1, double lambda2, double lam, double mu, double R,
double eps){

  double res[5];

  f90nrqcdderalpha_(&n, &l, &j, &s, charm, str, average, method, counting, &orderAlpha, &runAlpha, &order,
  &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &lambda1, &lambda2,
  &lam, &mu, &R, &eps, res);

   MLPutRealList(stdlink, res, 5);
   MLEndPacket(stdlink);

}

extern double f90massiter_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* mass, double* lambda1,
double* lambda2, double* lam, double* mu, double* R, double* res);

static void massiter(int n, int l, int j, int s, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha, int runAlpha, int order,
int run, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC, double mass, double lambda1, double lambda2,
double lam, double mu, double R){

  double res[5];

  f90massiter_(&n, &l, &j, &s, charm, str, average, method, counting, &orderAlpha, &runAlpha,
  &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &mass,
  &lambda1, &lambda2, &lam, &mu, &R, res);

   MLPutRealList(stdlink, res, 5);
   MLEndPacket(stdlink);

}

extern double f90massexpand_(int* n, int* l, int* j, int* s, char const* charm,
char const* str, char const* average, char const* method, char const* counting,
int* orderAlpha, int* runAlpha, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* mass, double* lambda1,
double* lambda2, double* lam, double* mu, double* R, double* res);

static void massexpand(int n, int l, int j, int s, char const* charm, char const* str,
char const* average, char const* method, char const* counting, int orderAlpha, int runAlpha, int order,
int run, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC, double mass, double lambda1, double lambda2,
double lam, double mu, double R){

  double res[5];

  f90massexpand_(&n, &l, &j, &s, charm, str, average, method, counting, &orderAlpha, &runAlpha,
  &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &mass,
  &lambda1, &lambda2, &lam, &mu, &R, res);

   MLPutRealList(stdlink, res, 5);
   MLEndPacket(stdlink);

}

extern double f90optimalr_(char const* type, double* n, char const* str,
int* orderAlpha, int* runAlpha, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* lambda, double* mu, double* res);

static double optimalr(char const* type, double n, char const* str,
int orderAlpha, int runAlpha, int order, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC,
double muC, double lambda, double mu){

  double res;

  f90optimalr_(type, &n, str, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz,
   &aMz, &mT, &muT, &mB, &muB, &mC, &muC,&lambda, &mu, &res);

  return res;
}

extern double f90optimalrvfns_(char const* up, char const* type, double* n,
char const* str, int* orderAlpha, int* runAlpha, int* order, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* lambda, double* mu1, double* mu2, double* res);

static double optimalrvfns(char const* up, char const* type, double n,
char const* str, int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double lambda, double mu1, double mu2){

  double res;

  f90optimalrvfns_(up, type, &n, str, &orderAlpha, &runAlpha, &order, &run,
  &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,&lambda, &mu1, &mu2, &res);

  return res;
}

extern double f90optimalr2_(double* n, int* orderAlpha, int* runAlpha,
int* order, int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT,
double* mB, double* muB, double* mC, double* muC, double* mass, double* res);

static double optimalr2(double n, int orderAlpha, int runAlpha, int order,
int run, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC, double mass){

  double res;

  f90optimalr2_(&n, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT,
  &muT, &mB, &muB, &mC, &muC, &mass, &res);

  return res;
}

extern double f90mmfrommsr_(char const* type, int* orderAlpha, int* runAlpha,
int* order, int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT,
double* mB, double* muB, double* mC, double* muC, double* mu, double* R,
double* res);

static double mmfrommsr(char const* type, int orderAlpha, int runAlpha, int order,
int run, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC, double mu, double R){

   double res;

   f90mmfrommsr_(type, &orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz,
   &mT, &muT, &mB, &muB, &mC, &muC, &mu, &R, &res);

  return res;
}

extern double f90jetmass_(int* orderAlpha, int* runAlpha, int* order, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* muLambda, double* R, double* mu, double* res);

static double jetmass(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double muLambda, double R, double mu){

   double res;

   f90jetmass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
               &muB, &mC, &muC, &muLambda, &R, &mu, &res);

  return res;
}

extern double f90mmfromjetmass_(int* orderAlpha, int* runAlpha, int* order, int* run,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* muLambda, double* R, double* mu, double* res);

static double mmfromjetmass(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double muLambda, double R, double mu){

   double res;

   f90mmfromjetmass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
               &muB, &mC, &muC, &muLambda, &R, &mu, &res);

  return res;
}

extern double f90deltamsbar_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static void deltamsbar(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res[4];

   f90deltamsbar_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
   &muC, &mu, res);

   MLPutRealList(stdlink, res, 4);
   MLEndPacket(stdlink);
}

extern double f90rhad_(char const* str, int* orderAlpha, int* runAlpha, int* order,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double * Q, double* res);

static double rhad(char const* str, int orderAlpha, int runAlpha, int order, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double mu, double Q){

   double res;

   f90rhad_(str, &orderAlpha, &runAlpha, &order, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
   &mC, &muC, &mu, &Q, &res);

  return res;
}

extern double f90sigmahad_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double * Q, double* res);

static double sigmahad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double mu, double Q){

   double res;

   f90sigmahad_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
   &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

  return res;
}

extern double f90sigmarad_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x, double* theta,
double* res);

static double sigmarad(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x, double theta){

   double res;

   f90sigmarad_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
   &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
   &Q, &x, &theta, &res);

  return res;
}

extern double f90sigmaradcum_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x0, double* x1, double* theta,
double* res);

static double sigmaradcum(char const* str, char const* curr, int orderAlpha,
int runAlpha, int order, int nf, double Mz, double gammaZ, double thetaW,
double aMz, double aMzQED, double mT, double muT, double mB, double muB,
double mC, double muC, double eH, double Q, double x0, double x1, double theta){

   double res;

   f90sigmaradcum_(str, curr, &orderAlpha, &runAlpha, &order, &nf, &Mz,
   &gammaZ, &thetaW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &eH,
   &Q, &x0, &x1, &theta, &res);

  return res;
}

extern double f90sigmaradcone_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ, double* thetaW,
double* aMz, double* aMzQED, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* eH, double * Q, double* x, double* theta,
double* deltaTheta, double* res);

static double sigmaradcone(char const* str, char const* curr, int orderAlpha,
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

extern double f90sigmaradconecum_(char const* str, char const* curr,
int* orderAlpha, int* runAlpha, int* order, int* nf, double* Mz, double* gammaZ,
double* thetaW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* theta,
double* x0, double* x1, double* deltaTheta, double* res);

static double sigmaradconecum(char const* str, char const* curr, int orderAlpha,
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

extern double f90rhadcoefs_(int * nf, double* res);

static void rhadcoefs(int nf){

   double res[4];

   f90rhadcoefs_(&nf, res);

   MLPutRealList(stdlink, res, 4);
   MLEndPacket(stdlink);
 }

extern double f90rhadmass_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double * Q, double* res);

static double rhadmass(char const* str, char const* curr, int orderAlpha, int runAlpha,
int runMass, int order, int nf, double Mz, double gammaZ, double sinW, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu, double Q){

   double res;

   f90rhadmass_(str, curr, &orderAlpha, &runAlpha, &runMass, &order, &nf, &Mz,
   &gammaZ, &sinW, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

  return res;
}

extern double f90sigmamass_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* mu, double * Q, double* res);

static double sigmamass(char const* str, char const* curr, int orderAlpha,
int runAlpha, int runMass, int order, int nf, double Mz, double gammaZ,
double sinW, double aMz, double aMzQED, double mT, double muT, double mB,
double muB, double mC, double muC, double mu, double Q){

  double res;

  f90sigmamass_(str, curr, &orderAlpha, &runAlpha, &runMass, &order, &nf, &Mz,
  &gammaZ, &sinW, &aMz, &aMzQED, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

 return res;
}

extern double f90sigmamassrad_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x,
double* theta, double* res);

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

extern double f90sigmamassradcum_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x0,
double* x1, double* theta, double* res);

static double sigmamassradcum(char const* str, char const* curr, int orderAlpha,
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

extern double f90sigmamassradcone_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* aMzQED, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* eH, double * Q, double* x,
double* theta, double* deltaTheta, double* res);

static double sigmamassradcone(char const* str, char const* curr, int orderAlpha,
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

extern double f90sigmamassradconecum_(char const* str, char const* curr,
int* orderAlpha, int* runAlpha, int* runMass, int* order, int* nf, double* Mz,
double* gammaZ, double* sinW, double* aMz, double* aMzQED, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* eH,
double * Q, double* x0, double* x1, double* theta, double* deltaTheta,
double* res);

static double sigmamassradconecum(char const* str, char const* curr, int orderAlpha,
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

extern double f90rqcd_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz,  double* aMz, double* mT, double* mu,
double * Q, double* res);

static double rqcd(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double aMz, double mT, double mu, double Q){

  double res;

  f90rqcd_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S, method,
  &lambda, &gt, &Mz, &aMz, &mT, &mu, &Q, &res);

 return res;
}

extern double f90rexp_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz,  double* aMz, double* mT, double* mu,
double* nu, double * Q, double* res);

static double rexp(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double aMz, double mT, double mu, double nu, double Q){

  double res;

  f90rexp_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S, method,
  &lambda, &gt, &Mz, &aMz, &mT, &mu, &nu, &Q, &res);

 return res;
}

extern double f90rmatched_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* aMz, double* mT, double* mu,
double* nu, double* v1, double* v2, double* Q, double* res);

static double rmatched(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double aMz, double mT, double mu, double nu, double v1, double v2,
double Q){

  double res;

  f90rmatched_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S, method,
  &lambda, &gt, &Mz, &aMz, &mT, &mu, &nu, &v1, &v2, &Q, &res);

 return res;
}

extern double f90sigmamatched_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* res);

static double sigmamatched(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double gammaZ, double sinW, double aMz, double aMzQED, double mT,
double mu, double nu, double v1, double v2, double Q){

  double res;

  f90sigmamatched_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S,
  method, &lambda, &gt, &Mz, &gammaZ, &sinW, &aMz, &aMzQED, &mT, &mu, &nu,
  &v1, &v2, &Q, &res);

 return res;
}

extern double f90sigmamatchedrad_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x, double* theta, double* res);

static double sigmamatchedrad(char const* str, int runAlpha, int runMass,
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

extern double f90sigmamatchedradcum_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x0, double* x1, double* theta, double* res);

static double sigmamatchedradcum(char const* str, int runAlpha, int runMass,
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

extern double f90sigmamatchedradcone_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x, double* theta, double* deltatheta, double* res);

static double sigmamatchedradcone(char const* str, int runAlpha, int runMass,
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

extern double f90sigmamatchedradconecum_(char const* str, int* runAlpha,
int* runMass, int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz, double* gammaZ, double* sinW,
double* aMz, double* aMzQED, double* mT, double* mu, double* nu, double* v1,
double* v2, double* Q, double* x0, double* x1, double* theta, double* deltatheta,
double* res);

static double sigmamatchedradconecum(char const* str, int runAlpha, int runMass,
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

extern double f90rmatchedlist_(char const* str, int* runAlpha, int* runMass,
int* ordMass, int* order, int* ord1S, double* R1S, char const* method,
double* lambda, double* gt, double* Mz,  double* aMz, double* mT, double* mu,
double* nu, double* v1, double* v2, double* Q0, double* Q1, double* deltaQ,
double* res);

static void rmatchedlist(char const* str, int runAlpha, int runMass, int ordMass,
int order, int ord1S, double R1S, char const* method, double lambda, double gt,
double Mz, double aMz, double mT, double mu, double nu, double v1, double v2,
double Q0, double Q1, double deltaQ){
  int imax = floor( (Q1 - Q0)/deltaQ ) + 1 ;
  double res[2 * imax];

  f90rmatchedlist_(str, &runAlpha, &runMass, &ordMass, &order, &ord1S, &R1S, method,
  &lambda, &gt, &Mz, &aMz, &mT, &mu, &nu, &v1, &v2, &Q0, &Q1, &deltaQ, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, 2 * imax);
   MLPutInteger(stdlink, 2);
   MLEndPacket(stdlink);

}

extern double f90hyper2f1_(double *a, double *b, double *c, double *x, double *res);

static double hyper2f1( double a, double b, double c, double x){
  double res;

   f90hyper2f1_(&a, &b, &c, &x, &res);

   return res;
}

extern double f90intecorre_(double *b, double *x0, double *x1, double *res);

static double intecorre( double b, double x0, double x1){
  double res;

   f90intecorre_(&b, &x0, &x1, &res);

   return res;
}

extern double f90hyperf32exact_(double *w, double *x, double *res);

static double hyperf32exact(double w, double x){
  double res;

   f90hyperf32exact_(&w, &x, &res);

   return res;
}

int main(int argc, char *argv[]){
    return MLMain(argc, argv);
}
