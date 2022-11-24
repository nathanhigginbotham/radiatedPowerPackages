// PatchAntenna.cxx

#include <cassert>

#include "Antennas/PatchAntenna.h"

#include "TVector3.h"
#include "TMath.h"

rad::PatchAntenna::PatchAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antYAx,
                                double width, double height, double f0,
                                double permittivity, double delay)
{
  antYAx = antYAx.Unit();
  antXAx = antXAx.Unit();

  // Make sure axes are perpendicular
  assert(antXAx.Dot(antYAx) == 0);

  antennaPosition = antPos;
  antennaXAxis = antXAx;
  antennaYAxis = antYAx;
  antennaZAxis = antXAx.Cross(antYAx);

  H = height;
  W = width;

  relativePerm = permittivity;
  centralFreq = f0;
  timeDelay = delay;

  // First calculate the effective perimittivity
  double epEff{((permittivity + 1) / 2) + ((permittivity - 1) / 2) * pow(1 + 12 * H / W, -0.5)};
  // Then calculate the effective length
  LEff = TMath::C() / (2 * f0 * sqrt(epEff));
  // From this, we can calculate the actual length of the patch
  double deltaL{0.412 * H * (epEff + 0.3) * (W / H + 0.264) / ((epEff - 0.258) * (W / H + 0.8))};
  L = LEff - 2 * deltaL;

  SetBandwidth();
}

TVector3 rad::PatchAntenna::GetETheta(const TVector3 electronPosition)
{
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  const double k = 2 * TMath::Pi() * centralFreq / TMath::C();

  double premult = TMath::Sin(k * W * TMath::Sin(thetaAng) * TMath::Sin(phiAng) / 2.0) * TMath::Cos(k * L * TMath::Sin(thetaAng) * TMath::Cos(phiAng) / 2) * TMath::Cos(phiAng) / (k * W * TMath::Sin(thetaAng) * TMath::Sin(phiAng) / 2);
  thetaHat *= premult;
  return thetaHat;
}

TVector3 rad::PatchAntenna::GetEPhi(const TVector3 electronPosition)
{
  TVector3 phiHat = GetPhiHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  const double k = 2 * TMath::Pi() * centralFreq / TMath::C();

  double premult = -1 * TMath::Sin(k * W * TMath::Sin(thetaAng) * TMath::Sin(phiAng) / 2.0) * TMath::Cos(k * L * TMath::Sin(thetaAng) * TMath::Cos(phiAng) / 2) * TMath::Cos(thetaAng) * TMath::Sin(phiAng) / (k * W * TMath::Sin(thetaAng) * TMath::Sin(phiAng) / 2);
  phiHat *= premult;
  return phiHat;
}
