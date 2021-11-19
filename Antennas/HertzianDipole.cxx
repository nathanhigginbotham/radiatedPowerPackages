// HertzianDipole.cxx

#include <cassert>

#include "Antennas/HertzianDipole.h"

#include "TVector3.h"

rad::HertzianDipole::HertzianDipole(TVector3 antPos, TVector3 antXAx, TVector3 antZAx,
				    const double freq) {
  // Make sure all axes are unit vectors initially
  antXAx = antXAx.Unit();
  antZAx = antZAx.Unit();

  // Make sure that axes are perpendicular to one another
  assert(antXAx.Dot(antZAx) == 0);
  
  antennaPosition = antPos;
  antennaXAxis = antXAx;
  antennaYAxis = antZAx.Cross(antXAx);
  antennaZAxis = antZAx;

  centralFreq = freq;
}

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


