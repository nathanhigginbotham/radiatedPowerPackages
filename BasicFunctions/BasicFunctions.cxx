// BasicFunctions.cxx
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TMath.h"

#include "FFTtools.h"
#include "FFTWComplex.h"

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

double rad::CalcAeHertzianDipole(const double wavelength, const TVector3 dipoleDir,
				 const TVector3 ePosition, const TVector3 position)
{
  double Ae = 3 * pow(wavelength, 2) / (8*TMath::Pi());
  const double psi = ((ePosition - position).Unit()).Angle(dipoleDir);
  Ae *= pow(TMath::Sin(psi), 2);
  return Ae;
}

double rad::CalcRetardedTime(const TVector3 fieldPoint, const TVector3 ePosition, const double labTime)
{
  double time = labTime - ((ePosition - fieldPoint).Mag() / TMath::C());
  return time;
}

// Very similar to the FFTtools implementation but without the scaling of the x axis the MHz
TGraph* rad::MakePowerSpectrumNorm(const TGraph* grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = FFTtools::doFFT(length, oldY);

  int newLength = (length/2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1/(deltaT*length);

  double tempF = 0;
  for(int i=0;i<newLength;i++) {
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry                                  
    power*=deltaT/(length); //For time-integral squared amplitude                               
    power/=deltaF;//Just to normalise bin-widths                                                   
    //Ends up the same as dt^2, need to integrate the power (multiply by df)                      
    //to get a meaningful number out.                                                            
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;
}

double rad::IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin=-1, Int_t lastBin=-1)
{
  double integral = FFTtools::integratePower(grFFT, firstBin, lastBin);
  // Multiply by frequency bin width
  double deltaF = grFFT->GetPointX(1) - grFFT->GetPointX(0);
  integral /= deltaF;
  return integral;
}
