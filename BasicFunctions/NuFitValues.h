/*
  NuFitValues.h
  
  Values from NuFit 5.1 analysis of oscillation data
  JHEP 09 (2020) 178 [arXiv:2007.14792]
*/

#ifndef NUFIT_VALUES_H
#define NUFIT_VALUES_H

#include "TMath.h"

namespace rad
{
  constexpr double kNuFitDmsq21 = 7.42e-5;
  
  // Assuming normal hierarchy
  constexpr double kNuFitTh12NH = 33.44 * TMath::Pi() / 180.0;
  constexpr double kNuFitTh23NH = 49.2 * TMath::Pi() / 180.0;
  constexpr double kNuFitTh13NH = 8.57 * TMath::Pi() / 180.0;
  constexpr double kNuFitdCPNH  = 194 * TMath::Pi() / 180.0;
  constexpr double kNuFitDmsq32NH = 2.515e-3 - kNuFitDmsq21;
      
  // Assuming inverted hierarchy
  constexpr double kNuFitTh12IH = 33.45 * TMath::Pi() / 180.0;
  constexpr double kNuFitTh23IH = 49.5 * TMath::Pi() / 180.0;
  constexpr double kNuFitTh13IH = 8.60 * TMath::Pi() / 180.0;
  constexpr double kNuFitdCPIH  = 287 * TMath::Pi() / 180.0;
  constexpr double kNuFitDmsq32IH = -2.498e-3;
}

#endif
