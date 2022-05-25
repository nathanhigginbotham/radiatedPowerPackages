# radiatedPowerPackages
Designed to offer easy calculation of EM fields and the associated voltage signals for a moving electron.

## Requirements
* ROOT 
* C++17 
* libRootFftwWrapper (available here https://github.com/nichol77/libRootFftwWrapper)
* Boost libraries

## Build instructions
The package is designed to be built with CMake and the build instructions are as follows.

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

## Creating electron trajectories
All the code required to create magnetic fields and electron trajectories is contained within the ```ElectronDynamics``` folder. 
It is possible to use the contained classes to generate your own electron trajectories (this is demonstrated in the ```writeTrajectory``` executable).
However, a helper class (```ElectronTrajectoryGen``` contained within ```ElectronDynamics/TrajectoryGen.h```) exists which will create ROOT files containing the trajectories in a format that can be read by other package components.

The required arguments are:
* The output file path
* A pointer to the magnetic field map (more on this later)
* A ROOT TVector3 of the initial electron position
* A ROOT TVector3 of the initial electron velocity
* The simulation step size to use in seconds (for reference an endpoint electron in a 1T field will have a cyclotron period of about $3.7 \times 10^{-11}$ seconds so 20 steps per orbit would be $1.85 \times 10^{-12}$. 
* The total time to simulate in seconds

The optional arguments are:
* The time (in seconds) at which to start the simulation
* The energy loss (which by default is zero). For realistic energy loss this should be set to $2 r_{e} / (3 c)$

Upon successfuly creation of the ```ElectronTrajectoryGen``` object, call the ```GenerateTraj``` member function to produce the file.

### Solver
The solver used to propagate the particles is the Boris solver which is energy conserving (except for any specified energy losses). 
An explanation of the solver (named as the Boris C solver) is given in [Ref. 1][1].

### Magnetic fields
Classes representing several types of magnetic fields and trapping configurations are found in ```ElectronDynamics/QTNMFields.h```.
The most familiar ones are ```BathtubField``` and ```HarmonicField```.

## Using the signal processing

### Antennas
Different types of antenna are implemented as derived classes of the abstract base class ```IAntenna```. 
Typically the one I used for my analyses is the ```HalfWaveDipole```.
Implementations do exist for other antennas but I provide no guarantee of accuracy for them (yet).
The ```HalfWaveDipole``` requires a spatial position, the directions of two (perpendicular) cartesian axes and a central frequency. 
Additionally a time delay (in seconds) may be specified for the antenna.

### FieldPoint
If you want to view the EM fields at a given point in space the best class to use is ```FieldPoint``` located in ```FieldClasses/FieldClasses.h``` which takes an its inputs a ROOT file containing an electron trajectory and a pointer to an antenna.
Calling the function ```GenerateFields``` is required to generate the fields between two specified times.

### InducedVoltage
The ```InducedVoltage``` voltage class provides a lighter weight implementation of just the voltage induced on the specific antenna (or array of antennas). 
The required inputs are an electron trajectory file and either a pointer to an antenna or (in the case of simulating an array of antennas wired together) and vector of antennas.
Additionally, a boolean can be specified over whether or not to include the signal propagation time when generating the signal.

If one just wants to view this signal without any downmixing, downsampling or noise added then simply call the function ```GenerateVoltage```. 
A ```TGraph``` of the time series signal can then be produced using ```GetVoltageGraph```.

### Signal
Full signal processing is done using the ```Signal``` class which takes as an input an instance of ```InducedVoltage```, an ```LocalOscillator``` used to define the down-mixing and the sampling rate (in Hertz).
Additionally a vector of ```GaussianNoise``` terms can be supplied as well as a maximum acquisition time.

[1]: <https://aip.scitation.org/doi/pdf/10.1063/1.5051077>
