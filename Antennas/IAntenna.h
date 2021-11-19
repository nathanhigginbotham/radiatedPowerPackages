/*
  IAntenna.h
  Contains the abstract base class for a generalised antenna
*/

#ifndef IANTENNA_H
#define IANTENNA_H

#include "TVector3.h"

namespace rad {

  class IAntenna {
    
  public:
    ~IAntenna(){}
    
    // Get the radiation patterns
    // Overidden in each antenna derived class
    virtual TVector3 GetETheta();
    virtual TVector3 GetEPhi();
    
  protected:
    TVector3 antennaPosition;
    TVector3 antennaXAxis;
    TVector3 antennaYAxis;
    TVector3 antennaZAxis;
    
    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the theta direction (relative to antenna axis) in global coords
    TVector3 GetThetaHat(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the phi direction (relative to antenna axis) in global coords
    TVector3 GetPhiHat(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The polar angle theta w.r.t. the antenna axis
    double GetTheta(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The azimuthal angle phi w.r.t. the antenna axis
    double GetPhi(const TVector3 electronPosition);
    
  };

}

#endif
