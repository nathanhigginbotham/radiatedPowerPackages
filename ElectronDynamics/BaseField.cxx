// BaseField.cxx

#include "ElectronDynamics/BaseField.h"

double rad::BaseField::evaluate_field_magnitude(const TVector3 vec) {
  TVector3 BField = evaluate_field_at_point(vec);
  return (BField.Mag());
}
