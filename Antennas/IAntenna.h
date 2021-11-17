/*
  IAntenna.h
  Contains the abstract base class for a generalised antenna
*/

#ifndef IANTENNA_H
#define IANTENNA_H

#include "TVector3.h"

namespace rad {

  class IAntenna {

  private:
    TVector3 antennaPosition;
    TVector3 antennaXAxis;
    TVector3 antennaYAxis;
    TVector3 antennaZAxis;
    
  public:
    ~IAntenna(){}

    /// \param antPos The position of the antenna
    /// \param antXAx Antenna defined X axis
    /// \param antYAx Antenna defined Y axis
    /// \param antZAx Antenna defined Z axis
    IAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antYAx, TVector3 antZAx);
    
    // Get the radiation patterns
    // Overidden in each antenna derived class
    virtual TVector3 GetETheta();
    virtual TVector3 GetEPhi();

  protected:
    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the theta direction (relative to antenna axis) in global coords
    TVector3 GetThetaHat(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the phi direction (relative to antenna axis) in global coords
    TVector3 GetPhiHat(const TVector3 electronPosition);
    
  };

}

#endif
