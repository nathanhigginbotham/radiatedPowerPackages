/*
  HertzianDipole.h
  Simplest case of a Hertzian dipole
*/

#ifndef HERTZIAN_DIPOLE_H
#define HERTZIAN_DIPOLE_H

#include "Antennas/IAntenna.h"

#include "TVector3.h"

namespace rad {

  class HertzianDipole : public IAntenna {
    
  public:
    TVector3 GetETheta(const TVector3 electronPosition);
    TVector3 GetEPhi(const TVector3 electronPosition);
  };
}

#endif
