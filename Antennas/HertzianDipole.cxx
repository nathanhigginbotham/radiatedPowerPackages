// HertzianDipole.cxx

#include <cassert>

#include "Antennas/HertzianDipole.h"

#include "TMath.h"
#include "TVector3.h"

rad::HertzianDipole::HertzianDipole(TVector3 antPos, TVector3 antXAx, TVector3 antZAx,
				    double freq, double delay) {
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
  timeDelay = delay;
  
  SetBandwidth();
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

double rad::HertzianDipole::GetETheta(double theta, double phi)
{
  return sin(theta);
}

double rad::HertzianDipole::GetEPhi(double theta, double phi)
{
  return 0;
}

double rad::HertzianDipole::GetHEff() {
  double heff = GetCentralWavelength() * TMath::Sqrt( 3 / (8*TMath::Pi()) );
  return heff;
}

double rad::HertzianDipole::GetAEff(TVector3 ePos)
{
  // Calculate angle between dipole direction and electron
  double thetaAng{GetTheta(ePos)};
  double lambda{TMath::C() / centralFreq};
  return (3 / (8 * TMath::Pi())) * pow(lambda * sin(thetaAng), 2);
}