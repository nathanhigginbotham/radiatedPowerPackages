/*
  CreateElectronTrajectories.cxx

  Skeleton code showing how to create electron trajectories for trapping geometries 
*/

// radiatedPowerPackages includes
#include "ElectronDynamics/TrajectoryGen.h"
#include "ElectronDynamics/QTNMFields.h"

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

// ROOT includes
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

using namespace rad;

int main(int argc, char *argv[])
{
  /// Create fields
  const double trapCoilRadius{0.03};       // Trapping coil radius in metres
  const double backgroundFieldStrength{1}; // Desired background field strength in tesla
  const double desiredTrapDepth{5e-3};     // 5 mT trap depth
  // Trap coils are simulated as current loops
  // Need to calculate the required current for the loop
  const double trapCoilCurrent{2 * desiredTrapDepth * trapCoilRadius / MU0 };

  // Create the harmonic field
  HarmonicField *harmonicField = new HarmonicField(trapCoilRadius, trapCoilCurrent, backgroundFieldStrength);
  // Create the bathtub field
  const double trapLength{0.3}; // trap length in metres
  TVector3 bkgFieldVector(0, 0, backgroundFieldStrength); // Vector of the background field (assume it's constant everywhere)
  BathtubField *bathtubField = new BathtubField(trapCoilRadius, trapCoilCurrent, -trapLength / 2, trapLength / 2, bkgFieldVector);

  /// Initial electron kinematics
  const double pitchAngleDeg{88};          // Electron pitch angle in degrees
  const double pitchAngleRad{pitchAngleDeg * TMath::Pi() / 180}; // The above in radians
  const double electronKE{18.5e3};                               // Electron kinetic energy in eV
  const double electronSpeed{GetSpeedFromKE(electronKE, ME)};    // Electron speed in m/s
  const double tau = 2 * R_E / (3 * TMath::C()); // Electron energy loss, set to 0 if you don't want electron to lose energy
  TVector3 V0(electronSpeed * sin(pitchAngleRad), 0, electronSpeed * cos(pitchAngleRad)); // Initial electron velocity vector
  
  // Calculate the electron gyroradius
  const double gyroradius{ GetGyroradius(V0, bathtubField->evaluate_field_at_point(TVector3(0, 0, 0)), ME) };
  TVector3 X0(0, -gyroradius, 0); // Initial electron position vector (offset means the gyroradius should be at the trap centre

  // Some electron simulation parameters
  const double simStepSize{1e-12}; // Simulation step size in seconds
  const double simTime{1e-6}; // Total simulation time in seconds
  
  TString harmonicFile{"harmonicTrackFile.root"}; // File path for harmonic file
  // Generate the trajectory for the harmonic trap
  ElectronTrajectoryGen harmTraj(harmonicFile, harmonicField, X0, V0, simStepSize, simTime, 0.0, tau);
  harmTraj.GenerateTraj();

  TString bathtubFile{"bathtubTrackFile.root"}; // File path for bathtub file  
  ElectronTrajectoryGen bathTraj(bathtubFile, bathtubField, X0, V0, simStepSize, simTime, 0.0, tau);
  bathTraj.GenerateTraj();
  
  return 0;
}
