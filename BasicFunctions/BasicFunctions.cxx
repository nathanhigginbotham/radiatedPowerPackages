// BasicFunctions.cxx

#include <cmath>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

#include "FFTtools.h"
#include "FFTWComplex.h"

// Electric field at the field point, calculated from Lienard-Wiechert potentials
ROOT::Math::XYZVector rad::CalcEField(const ROOT::Math::XYZPoint fieldPoint,
				      const ROOT::Math::XYZPoint ePosition,
				      const ROOT::Math::XYZVector eVelocity,
				      const ROOT::Math::XYZVector eAcceleration)
{  
  const ROOT::Math::XYZVector beta = eVelocity * (1.0 / TMath::C());
  const ROOT::Math::XYZVector betaDot = eAcceleration * (1.0 / TMath::C());
  double premult = TMath::Qe() / (4.0 * EPSILON0 * TMath::Pi());
  const double r = TMath::Sqrt((fieldPoint - ePosition).Mag2());
  const ROOT::Math::XYZVector rHat = (fieldPoint - ePosition).Unit();
  ROOT::Math::XYZVector term1 = (rHat - beta)*(1 - beta.Dot(beta)) * (1.0/( pow(1-beta.Dot(rHat) , 3) * r*r));
  ROOT::Math::XYZVector term2 = rHat.Cross( (rHat - beta).Cross(betaDot) ) * (1.0/( TMath::C()*r*pow(1-beta.Dot(rHat), 3)));
  ROOT::Math::XYZVector field = (term1 + term2) * premult;
  return field;
}

// Magnetic field at the field point, calculated from Lienard-Wiechert potentials
ROOT::Math::XYZVector rad::CalcBField(const ROOT::Math::XYZPoint fieldPoint,
				      const ROOT::Math::XYZPoint ePosition,
				      const ROOT::Math::XYZVector eVelocity,
				      const ROOT::Math::XYZVector eAcceleration)
{
  double premult = -1.0 * MU0 * TMath::Qe() / (4.0 * TMath::Pi());
  const ROOT::Math::XYZVector beta = eVelocity * (1.0 / TMath::C());
  const ROOT::Math::XYZVector betaDot = eAcceleration * (1.0 / TMath::C());
  const double r = TMath::Sqrt((fieldPoint - ePosition).Mag2());
  const ROOT::Math::XYZVector rHat = (fieldPoint - ePosition).Unit();
  ROOT::Math::XYZVector term1 = TMath::C() * rHat.Cross(beta) * (1 - beta.Dot(beta)) * (1.0 / ( r*r * pow(1.0 - beta.Dot(rHat), 3) ));
  ROOT::Math::XYZVector term2 = rHat.Cross(betaDot + rHat.Cross(beta.Cross(betaDot))) * (1.0 / (r * pow(1.0 - beta.Dot(rHat), 3) ));
  ROOT::Math::XYZVector field = premult * (term1 + term2);
  return field;
}

// Calculate Poynting vector
ROOT::Math::XYZVector rad::CalcPoyntingVec(const ROOT::Math::XYZPoint fieldPoint,
					   const ROOT::Math::XYZPoint ePosition,
					   const ROOT::Math::XYZVector eVelocity,
					   const ROOT::Math::XYZVector eAcceleration)
{
  ROOT::Math::XYZVector EField = CalcEField(fieldPoint, ePosition, eVelocity, eAcceleration);
  ROOT::Math::XYZVector BField = CalcBField(fieldPoint, ePosition, eVelocity, eAcceleration);
  ROOT::Math::XYZVector vec = EField.Cross(BField) * (1.0 / MU0);
  return vec;
}

ROOT::Math::XYZVector rad::CalcPoyntingVec(const ROOT::Math::XYZVector EField, const ROOT::Math::XYZVector BField)
{
  ROOT::Math::XYZVector vec = EField.Cross(BField) * (1.0 / MU0);
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

void rad::SetHistAttr(TH1 *h)
{
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH2 *h)
{
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

double rad::CalcAeHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
				 const ROOT::Math::XYZPoint ePosition,
				 const ROOT::Math::XYZPoint antennaPoint)			 
{
  double Ae = 3 * pow(wavelength, 2) / (8*TMath::Pi());
  const double psi = TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Ae *= pow(TMath::Sin(psi), 2);
  return Ae;
}

double rad::CalcAlHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
				 const ROOT::Math::XYZPoint ePosition,
				 const ROOT::Math::XYZPoint antennaPoint)
{
  double Al = wavelength * TMath::Sqrt(3 / (8*TMath::Pi()));
  const double psi = TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Al *= TMath::Sin(psi);
  return Al;
}

double rad::CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const double labTime)
{
  double time = labTime - TMath::Sqrt(((ePosition - fieldPoint).Mag2()) / TMath::C());
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
  setGraphAttr(grPower);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;
}

double rad::IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin, Int_t lastBin)
{
  double integral = FFTtools::integratePower(grFFT, firstBin, lastBin);
  // Multiply by frequency bin width
  double deltaF = grFFT->GetPointX(1) - grFFT->GetPointX(0);
  integral *= deltaF;
  return integral;
}

// Re-implementation of the filter from FFTtools but without and converstion factors
TGraph* rad::BandPassFilter(const TGraph* grWave, const double minFreq, const double maxFreq) {
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT = FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;
    double deltaF=1/(deltaT*length); //Hz

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      if(tempF<minFreq || tempF>maxFreq) {
	theFFT[i].re=0;
	theFFT[i].im=0;
      }      
      tempF+=deltaF;
    }

    double *filteredVals = FFTtools::doInvFFT(length,theFFT);

    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;
}

double rad::RayleighPDF(const double x, const double sigma) {
  double sigmaSq = sigma*sigma;
  double prob = (x / sigmaSq) * exp(-1*x*x / (2*sigmaSq));
  return prob;
}

double rad::RayleighCDF(const double x, const double sigma) {
  double f = 1.0 - exp( -1*x*x / (2*sigma*sigma) );
  return f;
}
