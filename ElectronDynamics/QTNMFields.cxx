// QTNMFields.cxx

#include "ElectronDynamics/QTNMFields.h"
#include "BasicFunctions/Constants.h"

#include <cmath>

#include "TVector3.h"
#include "TMath.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

rad::CoilField::CoilField(const double radius, const double current, const double z, const double mu) {
  coilRadius = radius;
  coilCurrent = current;
  coilZ = z;
  coilMu = mu;
}

double rad::CoilField::central_field() {
  return (coilCurrent * coilMu / coilRadius / 2.0);
}

double rad::CoilField::on_axis_field(const double z) {
  double field = coilMu * coilCurrent * pow(coilRadius, 2) / 2.0 / pow(coilRadius*coilRadius + pow(z - coilZ, 2), 1.5);
  return field;
}

TVector3 rad::CoilField::evaluate_field_at_point(const TVector3 vec) {
  double rad = sqrt(vec.X()*vec.X() + vec.Y()*vec.Y()); 

  // Is on axis
  if (rad / coilRadius < 1e-10) {
    TVector3 BField(0, 0, on_axis_field(vec.Z()));
    return BField;
  }

  // z relative to coil position
  double z_rel = vec.Z() - coilZ;

  double b_central = central_field();
  double rad_norm = rad / coilRadius;
  double z_norm = z_rel / coilRadius;
  double alpha = pow(1.0 + rad_norm, 2) + z_norm*z_norm;
  double root_alpha_pi = sqrt(alpha) * TMath::Pi();
  double beta = 4 * rad_norm / alpha;
  double int_k = boost::math::ellint_1(beta);
  double int_e = boost::math::ellint_2(beta);
  double gamma = alpha - 4 * rad_norm;

  double b_r = b_central * (int_e * ((1.0 + rad_norm*rad_norm + z_norm*z_norm) / gamma) - int_k) / root_alpha_pi *(z_rel / rad);
  double b_z = b_central * (int_e * ((1.0 - rad_norm*rad_norm - z_norm*z_norm) / gamma) + int_k) / root_alpha_pi;

  TVector3 BField(b_r*vec.X()/rad, b_r*vec.Y()/rad, b_z);
  return BField;
}
