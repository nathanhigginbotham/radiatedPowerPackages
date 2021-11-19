/*
  PatchAntenna.h
  Implementation of a microstrip patch antenna of variable width and length
*/

#ifndef PATCH_ANTENNA_H
#define PATCH_ANTENNA_H

#include "Antennas/IAntenna.h"
#include "TVector3.h"

namespace rad {

  class PatchAntenna : public IAntenna {

  private:
    double L;
    double W;
    double relativePerm;
    
  public:
    PatchAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antZAx,
		 const double length, const double width, const double permittivity=1.0);
    
    TVector3 GetETheta(const TVector3 electronPosition);
    TVector3 GetEPhi(const TVector3 electronPosition);

    double GetHEff();
  };
}



#endif
