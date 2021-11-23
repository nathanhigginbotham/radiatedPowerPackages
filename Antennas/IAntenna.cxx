// IAntenna.cxx

#include <cassert>

#include "Antennas/IAntenna.h"

#include "TVector3.h"
#include "TMath.h"

double rad::IAntenna::GetCentralWavelength() {
  double freq = GetCentralFrequency();
  double lambda = TMath::C() / freq;
  return lambda;
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

double rad::IAntenna::GetTheta(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double theta = TMath::ACos(antennaZAxis.Dot(rHat));
  return theta;
}

double rad::IAntenna::GetPhi(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double phi = TMath::ATan2(antennaYAxis.Dot(rHat), antennaXAxis.Dot(rHat));
  return phi;
}

void rad::IAntenna::SetBandwidth(const double lowerLimit, const double upperLimit) {
  assert(lowerLimit < upperLimit);
  lowerBandwidth = lowerLimit;
  upperBandwidth = upperLimit;
}

