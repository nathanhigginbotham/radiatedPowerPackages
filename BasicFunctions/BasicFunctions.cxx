// BasicFunctions.cxx
#include "BasicFunctions.h"

#include "TVector3.h"
#include "TMath.h"

TVector3 CalcEField(const TVector3 sourcePosition, const TVector3 ePosition,
		    const TVector3 eVelocity, const TVector3 eAcceleration)
{
  const double epsilon0 = 8.8541878128e-12;
  
  const TVector3 beta = eVelocity * (1.0 / TMath::C());
  const TVector3 betaDot = eAcceleration * (1.0 / TMath::C());
  double premult = TMath::Qe() / (4.0 * epsilon0 * TMath::Pi());
  const double r = (sourcePosition - ePosition).Mag();
  const TVector3 rHat = (sourcePosition - ePosition).Unit();
  TVector3 term1 = (rHat - beta)*(1 - beta.Dot(beta)) * (1.0/( pow(1-beta.Dot(rHat) , 3) * r*r));
  TVector3 term2 = rHat.Cross( (rHat - beta).Cross(betaDot) ) * (1.0/( TMath::C()*r*pow(1-beta.Dot(rHat), 3)));
  TVector3 field = (term1 + term2) * premult;
  return field;
}
