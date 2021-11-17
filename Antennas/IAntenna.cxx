// IAntenna.cxx

#include <cassert>

#include "Antennas/IAntenna.h"

#include "TVector3.h"
#include "TMath.h"

rad::IAntenna::IAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antYAx, TVector3 antZAx) {
  // Make sure all axes are unit vectors initially
  antXAx = antXAx.Unit();
  antYAx = antYAx.Unit();
  antZAx = antZAx.Unit();

  // Make sure that axes are all perpendicular to one another
  assert(antXAx.Dot(antYAx) == 0 && antYAx.Dot(antZAx) == 0 && antXAx.Dot(antZAx) == 0);
  
  antennaPosition = antPos;
  antennaXAxis = antXAx;
  antennaYAxis = antYAx;
  antennaZAxis = antZAx;
}

TVector3 rad::IAntenna::GetThetaHat(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  const double theta = TMath::ACos(antennaZAxis.Dot(rHat));
  const double phi = TMath::ATan2(antennaYAxis.Dot(rHat), antennaXAxis.Dot(rHat));
  TVector3 vec(TMath::Cos(theta)*TMath::Cos(phi), TMath::Cos(theta)*TMath::Sin(phi), -1*TMath::Sin(theta));
  vec.RotateUz(antennaZAxis);
  return vec;
}

TVector3 rad::IAntenna::GetPhiHat(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  const double phi = TMath::ATan2(antennaYAxis.Dot(rHat), antennaXAxis.Dot(rHat));
  TVector3 vec(-1*TMath::Sin(phi), TMath::Cos(phi), 0);
  vec.RotateUz(antennaZAxis);
  return vec;
}
