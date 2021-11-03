// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "TVector3.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "TH1.h"
#include "TH2.h"

namespace rad
{

  TVector3 CalcEField(const TVector3 fieldPoint, const TVector3 ePosition,
		      const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcBField(const TVector3 fieldPoint, const TVector3 ePosition,
		      const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcPoyntingVec(const TVector3 fieldPoint, const TVector3 ePosition,
			   const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcPoyntingVec(const TVector3 EField, const TVector3 BField);

  void setGraphAttr(TGraph *gr);
  
  void SetHistAttr(TH1 *h);
  
  void SetHistAttr(TH2 *h);

  double CalcAeHertzianDipole(const double wavelength, const TVector3 dipoleDir,
			      const TVector3 ePosition, const TVector3 antennaPoint);
  
  double CalcAlHertzianDipole(const double wavelength, const TVector3 dipoleDir,
			      const TVector3 ePosition, const TVector3 antennaPoint);

  double CalcRetardedTime(const TVector3 fieldPoint, const TVector3 ePosition, const double labTime);

  // Produces power spectrum with the desired normalisation
  TGraph* MakePowerSpectrumNorm(const TGraph* grWave);

  // Integrate the power spectrum
  double IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin=-1, Int_t lastBin=-1);
  
}
 
#endif
