// HertzianDipole.cxx

#include "Antennas/HertzianDipole.h"

#include "TVector3.h"

// Calculate the radiation pattern in the theta hat direction
TVector3 rad::HertzianDipole::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  thetaHat *= TMath::Sin(thetaAng);
  return thetaHat;
}

TVector3 rad::HertzianDipole::GetEPhi(const TVector3 electronPosition) {
  // No radiation in the phi direction for a hertzian dipole
  return TVector3(0, 0, 0);
}


