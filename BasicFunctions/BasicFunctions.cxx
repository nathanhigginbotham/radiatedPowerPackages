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
#include "TF1.h"
#include "TRandom3.h"

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

TGraph* rad::MakePowerSpectrumPeriodogram(const TGraph* grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = FFTtools::doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length/2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1/(deltaT*length);

  double tempF = 0;
  for(int i=0;i<newLength;i++) {
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry
    double scale = lengthDub*lengthDub;
    power /= scale;
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;    
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
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

// Re-implementation of the filter from FFTtools but without the conversion factors
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

double rad::RayleighPDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighPDF(xx, par[0]);
  return retVal;
}

double rad::RayleighCDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighCDF(xx, par[0]);
  return retVal;
}

void rad::AddWhiteNoiseFrequencyDomainPowerNorm(TGraph* grIn, const double Teff, const int seed) {
  gRandom->SetSeed(seed);
  const double sampleRate = 2*grIn->GetPointX(grIn->GetN()-1);
  const double deltaT = 1.0 / sampleRate;
  const double deltaF = grIn->GetPointX(1) - grIn->GetPointX(0);
  const double sigma = TMath::Sqrt( TMath::K() * Teff * sampleRate );
  TF1* f1 = new TF1("f1", RayleighPDFFunc, 0, 4*sigma, 1);
  f1->SetParameter(0, sigma);

  // Now loop through the bins of the graph and add the noise
  for (int i = 0; i < grIn->GetN(); i++) {
    double noise = pow(f1->GetRandom()*sqrt(0.5), 2)*(1/deltaF)*(deltaT);
    grIn->SetPointY(i, grIn->GetPointY(i) + noise);
  }
  
  delete f1;
}

TVector3 rad::calculate_omega(const TVector3 BField, const double charge, const double energy, const double mass) {
  double gamma_m0 = mass + energy * TMath::Qe() / pow(TMath::C(), 2);
  return (charge * BField * (1.0 / gamma_m0));
}

TGraph* rad::DownmixInPhase(TGraph* grInput, const double freq) {
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(i, time, grInput->GetPointY(i)*TMath::Cos(2*TMath::Pi()*freq*time));
  }
  return grOut;
}

TGraph* rad::DownmixQuadrature(TGraph* grInput, const double freq) {
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(i, time, grInput->GetPointY(i)*TMath::Sin(2*TMath::Pi()*freq*time));
  }
  return grOut;
}

void rad::ScaleGraph(TGraph* grInput, const double scale) {
  for (int i = 0; i < grInput->GetN(); i++) {
    grInput->SetPointY(i, grInput->GetPointY(i) * scale);
  }
}

double rad::CalcCyclotronFreq(const double KE, const double B) {
  double freq = TMath::Qe()*B / (ME + (KE*TMath::Qe()/pow(TMath::C(), 2)));
  return freq / (2*TMath::Pi());
}
