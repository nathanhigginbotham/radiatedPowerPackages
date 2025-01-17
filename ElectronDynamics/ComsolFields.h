/*
    ComsolFields.h

    Class used for realistic fields from COMSOL simulations
    Values come from M. Fleck's simulations

    S. Jones 26-10-2022
*/

#ifndef COMSOL_FIELDS_H
#define COMSOL_FIELDS_H

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/QTNMFields.h"

#include "TGraph2D.h"
#include "TVector3.h"

#include <string>

namespace rad
{
    class ComsolField : public BaseField
    {
    private:
        TGraph2D *fieldValues = 0;
        double scaleFactor;

    public:
        /// Parametrised constructor
        /// \param fieldFile CSV file path
        /// \param centralField The desired central magnetic field
        ComsolField(std::string fieldFile, double centralField = 0.0);

        /// Destructor
        ~ComsolField();

        /// Calculates the magnetic field at a point in space
        /// \param vec The position vector (units of metres)
        /// \return The magnetic field vector (units of tesla)
        TVector3 evaluate_field_at_point(const TVector3 vec);
    };

    class ComsolHarmonicField : public BaseField
    {
    private:
        CoilField coil;
        TGraph2D *fieldValues = 0;
        double scaleFactor;

    public:
        /// Parametrised constructor
        /// \param fieldFile CSV file path
        /// \param centralField The desired central magnetic field
        ComsolHarmonicField(double radius, double current,
                            std::string fieldFile, double centralField = 0.0);

        /// Destructor
        ~ComsolHarmonicField();

        /// Calculates the magnetic field at a point in space
        /// \param vec The position vector (units of metres)
        /// \return The magnetic field vector (units of tesla)
        TVector3 evaluate_field_at_point(const TVector3 vec);
    };
}

#endif