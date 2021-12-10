// LockInAmplifier.cxx

#include "SignalProcessing/LockInAmplifier.h"
#include "SignalProcessing/LocalOscillator.h"
#include "BasicFunctions/BasicFunctions.h"

#include "TGraph.h"

#include "FFTtools.h"

#include <cmath>

rad::LockInAmplifier::LockInAmplifier(LocalOscillator referenceFreq) {
  refFreq = referenceFreq;
}

TGraph* rad::LockInAmplifier::GetAmplitudeTimeDomain(Signal sig) {
  TGraph* grV = sig.GetVITimeDomain();
  // Use the existing downmixing functions to get the in phase and quadrature components
  TGraph* VI = sig.DownmixInPhase(grV, refFreq);
  TGraph* VQ = sig.DownmixQuadrature(grV, refFreq);
  delete grV;
  TGraph* VIFilter = BandPassFilter(VI, 0, 375e6);
  TGraph* VQFilter = BandPassFilter(VQ, 0, 375e6);
  delete VI;
  delete VQ;
  TGraph* grR = new TGraph();
  setGraphAttr(grR);
  grR->GetYaxis()->SetTitle("R");
  grR->GetXaxis()->SetTitle("Time [s]");
  
  for (int i = 0; i < VIFilter->GetN(); i++) {
    grR->SetPoint(grR->GetN(), VIFilter->GetPointX(i), sqrt(pow(VIFilter->GetPointY(i), 2) + pow(VQFilter->GetPointY(i), 2)));
  }
  
  delete VIFilter;
  delete VQFilter;  
  return grR;
}

TGraph* rad::LockInAmplifier::GetPhaseTimeDomain(Signal sig) {
  TGraph* grV = sig.GetVITimeDomain();
  // Use the existing downmixing functions to get the in phase and quadrature components
  TGraph* VI = sig.DownmixInPhase(grV, refFreq);
  TGraph* VQ = sig.DownmixQuadrature(grV, refFreq);
  delete grV;
  TGraph* VIFilter = BandPassFilter(VI, 0, 375e6);
  TGraph* VQFilter = BandPassFilter(VQ, 0, 375e6);
  delete VI;
  delete VQ;
  TGraph* grTheta = new TGraph();
  setGraphAttr(grTheta);
  grTheta->GetYaxis()->SetTitle("R");
  grTheta->GetXaxis()->SetTitle("Time [s]");
  
  for (int i = 0; i < VIFilter->GetN(); i++) {
    grTheta->SetPoint(grTheta->GetN(), VIFilter->GetPointX(i), atan2(VQFilter->GetPointY(i), VIFilter->GetPointY(i)));
  }
  
  delete VIFilter;
  delete VQFilter;  
  return grTheta;  
}
