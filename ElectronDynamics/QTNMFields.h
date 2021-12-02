/*
  QTNMFields.h

  Derived QTNM field classes

  Transcription of T. Goffrey's original Python code 
*/

#ifndef QTNM_FIELDS_H
#define QTNM_FIELDS_H

#include "ElectronDynamics/BaseField.h"
#include "BasicFunctions/Constants.h"

#include "TVector3.h"

namespace rad
{
  /// Class describing the magnetic field of an infinitesimally thin current coil
  /// in the x-y plane
  class CoilField : public BaseField {
  private:
    double coilRadius;
    double coilCurrent;
    double coilZ;
    double coilMu;

    double central_field();

    double on_axis_field(const double z);
    
  public:
    /// \param radius The radius of the coil (in metres)
    /// \param current The current passing through the coil (in Amps)
    /// \param z The z position of the coil (in metres)
    /// \param Magnetic permeability
    CoilField(const double radius=0.005, const double current=40, const double z=0.0, const double mu=MU0);

    /// \param vec Position vector of charge
    /// \returns The magnetic field vector at the point
    TVector3 evaluate_field_at_point(const TVector3 vec);
  };

  /// Class describing the bathtub field
  /// Made up of two superposed current loops and a background magnetic field
  class BathtubField : public BaseField {
  private:
    CoilField coil1;
    CoilField coil2;
    TVector3 btBkg;
    
  public:
    BathtubField(const double radius, const double current, const double Z1, const double Z2, TVector3 background);

    TVector3 evaluate_field_at_point(const TVector3 vec);
  };
}

#endif
