/*
  HalfWaveDipole.h
  Hypothetical antenna with a length of half a wavelength
*/

#ifndef HALF_WAVE_DIPOLE_H
#define HALF_WAVE_DIPOLE_H

#include "Antennas/IAntenna.h"

#include "TVector3.h"

namespace rad {

  class HalfWaveDipole : public IAntenna {
    
  public:
    TVector3 GetETheta(const TVector3 electronPosition);
    TVector3 GetEPhi(const TVector3 electronPosition);
  };
}

#endif
