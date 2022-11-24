/*
  PatchAntenna.h
  Implementation of a microstrip patch antenna of variable width and length
*/

#ifndef PATCH_ANTENNA_H
#define PATCH_ANTENNA_H

#include "Antennas/IAntenna.h"
#include "TVector3.h"

namespace rad
{
  class PatchAntenna : public IAntenna
  {

  private:
    double L;            // Length of the patch (in metres)
    double W;            // Width of the patch antenna (in metres)
    double H;            // Height of the substrate (in metres)
    double relativePerm; // Relative permettivity
    double LEff;         // Effective length (in metres)

  public:
    PatchAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antYAx,
                 double width, double height, double f0,
                 double permittivity = 1.0, double delay = 0.0);

    /// Calculates the electric field of the antenna in the theta direction
    /// \param electronPosition Electron position (in metres)
    /// \return The electric field (in V/m)
    TVector3 GetETheta(const TVector3 electronPosition);

    /// Calculates the electric field of the antenna in the phi direction
    /// \param electronPosition Electron position (in metres)
    /// \return The electric field (in V/m)
    TVector3 GetEPhi(const TVector3 electronPosition);

    double GetHEff() { return LEff; }

    double GetL() { return L; }

    double GetH() {return H; }

    double GetW() { return W; }
  };
}

#endif
