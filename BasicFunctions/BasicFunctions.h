// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "BasicFunctions/Constants.h"

#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "TH1.h"
#include "TH2.h"

#include <vector>

namespace rad
{

  ROOT::Math::XYZVector CalcEField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);

  ROOT::Math::XYZVector CalcBField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);
 
  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);

  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZVector EField, const ROOT::Math::XYZVector BField);

  void setGraphAttr(TGraph *gr);
  
  void SetHistAttr(TH1 *h);
  
  void SetHistAttr(TH2 *h);

  double CalcAeHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
			      const ROOT::Math::XYZPoint ePosition,
			      const ROOT::Math::XYZPoint antennaPoint);
  
  double CalcAlHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
			      const ROOT::Math::XYZPoint ePosition,
			      const ROOT::Math::XYZPoint antennaPoint);
  
  double CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const double labTime);

  // Produces power spectrum with the desired normalisation
  TGraph* MakePowerSpectrumNorm(const TGraph* grWave);

  // Produdces power spectrum normalised to the mean squared amplitude of the time domain
  TGraph* MakePowerSpectrumPeriodogram(const TGraph* grWave);  
  
  // Integrate the power spectrum
  double IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin=-1, Int_t lastBin=-1);

  // Implements a simple band pass filter
  TGraph* BandPassFilter(const TGraph* grWave, const double minFreq, const double maxFreq);

  // PDF of the Rayleigh distribution
  double RayleighPDF(const double x, const double sigma);

  // CDF of the Rayleigh distribution
  double RayleighCDF(const double x, const double sigma);

  double RayleighPDFFunc(double *x, double *par);

  double RayleighCDFFunc(double *x, double *par);
  
  void AddWhiteNoiseFrequencyDomainPowerNorm(TGraph* grIn, const double Teff, const int seed=0);

  /// \param BField A 3-vector describing the magnetic field at a point in space
  /// \param charge Particle charge. Default is electron charge
  /// \param energy Particle kinetic energy in electronvolts
  /// \param mass Particle mass in kilograms. Default is electron rest mass
  /// \Returns the cyclotron frequency
  TVector3 calculate_omega(const TVector3 BField, const double charge=-TMath::Qe(), const double energy=0.0, const double mass=ME);

  /// Downmixes a time series of data with the appropriate frequency
  /// \param grInput The input time series data
  /// \param freq The frequency in Hertz with which to downmix
  /// \Returns The downmixed time series
  TGraph* DownmixInPhase(TGraph* grInput, const double freq);

  /// Downmixes a time series of data with the appropriate frequency
  /// \param grInput The input time series data
  /// \param freq The frequency in Hertz with which to downmix
  /// \Returns The downmixed time series
  TGraph* DownmixQuadrature(TGraph* grInput, const double freq);
  
  /// Scales all points of a TGraph by a constant
  /// \param grInput The input graph
  /// \param scale factor
  /// \Returns grInput * scale
  void ScaleGraph(TGraph* grInput, const double scale);

  /// Calculates relativistic electron cyclotron frequency
  /// \param KE The electron kinetic energy in eV
  /// \param B Magnetic field strength in Tesla
  double CalcCyclotronFreq(const double KE, const double B=1.0);

  /// Sums the points in a number of TGraphs
  /// \param grInput A vector of the graphs to be summed
  /// \Returns A graph containing the summed points
  TGraph* SumGraphs(std::vector<TGraph*> grInput);

  /// Downsamples a time series graph at a given sample rate using linear interpolation
  /// \param grInput The input graph to be sampled
  /// \param sRate The sample rate
  /// \Returns The downsampled graph
  TGraph* SampleWaveform(TGraph* grInput, const double sRate);

  /// Converts an input TGraph to a histogram
  /// \param grInput Input graph to be converted
  /// \Returns The converted histogram
  TH1D* GraphToHistogram(TGraph* grInput);

  /// Dowmixes, filters and downsamples a time series graph
  /// \param grInput The inputted time domain voltage graph
  /// \param downmixFreq The frequency at which to downmix, in Hertz
  /// \param sampleRate The frequency at which to sample, in Hertz
  /// \Returns The downmixed, filtered and sampled time domain graph
  TGraph* SignalProcessGraph(TGraph* grInput, const double downmixFreq, const double sampleRate);
}
 
#endif
