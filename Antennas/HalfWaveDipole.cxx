// HalfWaveDipole.cxx

#include "Antennas/HalfWaveDipole.h"

#include "TVector3.h"

// Calculate the radiation pattern in the theta hat direction
TVector3 rad::HalfWaveDipole::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  thetaHat *= TMath::Cos(TMath::Pi() * TMath::Cos(thetaAng) / 2) / TMath::Sin(thetaAng);
  return thetaHat;
}

TVector3 rad::HalfWaveDipole::GetEPhi(const TVector3 electronPosition) {
  // No radiation in the phi direction for a hertzian dipole
  return TVector3(0, 0, 0);
}
