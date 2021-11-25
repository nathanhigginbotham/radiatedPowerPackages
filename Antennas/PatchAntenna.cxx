// PatchAntenna.cxx

#include <cassert>

#include "Antennas/PatchAntenna.h"

#include "TVector3.h"
#include "TMath.h"

rad::PatchAntenna::PatchAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antZAx,
				const double length, const double width, const double permittivity) {
  antXAx = antXAx.Unit();
  antZAx = antZAx.Unit();

  // Make sure axes are perpendicular
  assert(antXAx.Dot(antZAx) == 0);

  antennaPosition = antPos;
  antennaXAxis = antXAx;
  antennaZAxis = antZAx;
  antennaYAxis = antZAx.Cross(antXAx);
  
  L = length;
  W = width;
  relativePerm = permittivity;
  centralFreq = TMath::C() / (2.0 * L * TMath::Sqrt(relativePerm));
  
  SetBandwidth();
}

TVector3 rad::PatchAntenna::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  const double k = 2 * TMath::Pi() * centralFreq / TMath::C();

  double premult = TMath::Sin( k*W*TMath::Sin(thetaAng)*TMath::Sin(phiAng)/2.0 ) * TMath::Cos( k*L*TMath::Sin(thetaAng)*TMath::Cos(phiAng)/2 ) * TMath::Cos(phiAng) / (k*W*TMath::Sin(thetaAng)*TMath::Sin(phiAng)/2);
  thetaHat *= premult;
  return thetaHat;
}

TVector3 rad::PatchAntenna::GetEPhi(const TVector3 electronPosition) {
  TVector3 phiHat = GetPhiHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  const double k = 2 * TMath::Pi() * centralFreq / TMath::C();

  double premult = -1 * TMath::Sin( k*W*TMath::Sin(thetaAng)*TMath::Sin(phiAng)/2.0 ) * TMath::Cos( k*L*TMath::Sin(thetaAng)*TMath::Cos(phiAng)/2 ) * TMath::Cos(thetaAng)*TMath::Sin(phiAng) / (k*W*TMath::Sin(thetaAng)*TMath::Sin(phiAng)/2);
  phiHat *= premult;
  return phiHat;
}

double rad::PatchAntenna::GetHEff() {
  double LEff = TMath::C() / (2*GetCentralFrequency());
  return LEff;
}
