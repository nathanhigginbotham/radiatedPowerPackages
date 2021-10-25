// BasicFunctions.cxx
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TMath.h"

// Electric field at the field point, calculated from Lienard-Wiechert potentials
TVector3 rad::CalcEField(const TVector3 fieldPoint, const TVector3 ePosition,
			 const TVector3 eVelocity, const TVector3 eAcceleration)
{  
  const TVector3 beta = eVelocity * (1.0 / TMath::C());
  const TVector3 betaDot = eAcceleration * (1.0 / TMath::C());
  double premult = TMath::Qe() / (4.0 * EPSILON0 * TMath::Pi());
  const double r = (fieldPoint - ePosition).Mag();
  const TVector3 rHat = (fieldPoint - ePosition).Unit();
  TVector3 term1 = (rHat - beta)*(1 - beta.Dot(beta)) * (1.0/( pow(1-beta.Dot(rHat) , 3) * r*r));
  TVector3 term2 = rHat.Cross( (rHat - beta).Cross(betaDot) ) * (1.0/( TMath::C()*r*pow(1-beta.Dot(rHat), 3)));
  TVector3 field = (term1 + term2) * premult;
  return field;
}

// Magnetic field at the field point, calculated from Lienard-Wiechert potentials
TVector3 rad::CalcBField(const TVector3 fieldPoint, const TVector3 ePosition,
			 const TVector3 eVelocity, const TVector3 eAcceleration)
{
  double premult = -1.0 * MU0 * TMath::Qe() / (4.0 * TMath::Pi());
  const TVector3 beta = eVelocity * (1.0 / TMath::C());
  const TVector3 betaDot = eAcceleration * (1.0 / TMath::C());
  const double r = (fieldPoint - ePosition).Mag();
  const TVector3 rHat = (fieldPoint - ePosition).Unit();
  TVector3 term1 = TMath::C() * rHat.Cross(beta) * (1 - beta.Dot(beta)) * (1.0 / ( r*r * pow(1.0 - beta.Dot(rHat), 3) ));
  TVector3 term2 = rHat.Cross(betaDot + rHat.Cross(beta.Cross(betaDot))) * (1.0 / (r * pow(1.0 - beta.Dot(rHat), 3) ));
  TVector3 field = premult * (term1 + term2);
  return field;
}

// Calculate Poynting vector
TVector3 rad::CalcPoyntingVec(const TVector3 fieldPoint, const TVector3 ePosition,
			      const TVector3 eVelocity, const TVector3 eAcceleration)
{
  TVector3 EField = CalcEField(fieldPoint, ePosition, eVelocity, eAcceleration);
  TVector3 BField = CalcBField(fieldPoint, ePosition, eVelocity, eAcceleration);
  TVector3 vec = EField.Cross(BField) * (1.0 / MU0);
  return vec;
}

TVector3 rad::CalcPoyntingVec(const TVector3 EField, const TVector3 BField)
{
  TVector3 vec = EField.Cross(BField) * (1.0 / MU0);
  return vec;
}

void rad::setGraphAttr(TGraph *gr)
{
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}
