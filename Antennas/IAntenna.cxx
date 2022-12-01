// IAntenna.cxx

#include <cassert>
#include <iostream>

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

double rad::IAntenna::GetPatternIntegral()
{
  // Now calculate the surface integral of the radiation pattern
  // Need this for proper calculation of the gain
  const int nPntsTheta{200};
  const int nPntsPhi{200};
  // Hypothetical width for integration
  const double binArea{TMath::Pi() * 2 * TMath::Pi() / double(nPntsPhi * nPntsTheta)}; 
  const double binWidthTheta{TMath::Pi() / double(nPntsTheta)};
  const double binWidthPhi{2 * TMath::Pi() / double(nPntsPhi)};
  double PRad{0};
  for (int ith{0}; ith < nPntsTheta; ith++)
  {
    double theta{double(ith) * binWidthTheta + binWidthTheta / 2};
    for (int iph{0}; iph < nPntsPhi; iph++)
    {
      double phi{double(iph) * binWidthPhi + binWidthPhi / 2};
      double uSin{(GetETheta(theta, phi) * GetETheta(theta, phi) + 
                  GetEPhi(theta, phi) * GetEPhi(theta, phi)) * sin(theta)};
      PRad += uSin * binArea;
    }
  }
  std::cout<<"PRad = "<<PRad<<std::endl;
  return PRad;
}

