// PatchAntenna.cxx

#include <cassert>

#include "Antennas/PatchAntenna.h"

#include "boost/math/special_functions/sinc.hpp"

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

  PRad = GetPatternIntegral();
}

TVector3 rad::PatchAntenna::GetETheta(const TVector3 electronPosition)
{
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  return thetaHat * GetETheta(thetaAng, phiAng);
}

TVector3 rad::PatchAntenna::GetEPhi(const TVector3 electronPosition)
{
  TVector3 phiHat = GetPhiHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  return phiHat * GetEPhi(thetaAng, phiAng);
}

double rad::PatchAntenna::GetETheta(double theta, double phi)
{
  double k{2 * TMath::Pi() * centralFreq / TMath::C()};
  double F1{boost::math::sinc_pi(k * H * sin(theta) * cos(theta) / 2.0) *
            boost::math::sinc_pi(k * W * sin(theta) * sin(phi) / 2.0)};
  double F2{2.0 * cos(k * L * sin(theta) * cos(phi) / 2.0)};
  double ETheta{cos(phi) * F1 * F2};
  return ETheta;
}

double rad::PatchAntenna::GetEPhi(double theta, double phi)
{
  double k{2 * TMath::Pi() * centralFreq / TMath::C()};
  double F1{boost::math::sinc_pi(k * H * sin(theta) * cos(theta) / 2.0) *
            boost::math::sinc_pi(k * W * sin(theta) * sin(phi) / 2.0)};
  double F2{2.0 * cos(k * L * sin(theta) * cos(phi) / 2.0)};
  double EPhi{cos(theta) * sin(phi) * F1 * F2};
  return EPhi;
}

double rad::PatchAntenna::GetAEff(TVector3 ePos)
{
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() * (GetETheta(theta, phi) * GetETheta(theta, phi)) +
              GetEPhi(theta, phi) * GetEPhi(theta, phi) / PRad};
  return pow(TMath::C() / centralFreq, 2) * gain / (4 * TMath::Pi());
}

