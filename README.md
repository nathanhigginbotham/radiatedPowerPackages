# radiatedPowerPackages
Designed to offer easy calculation of electric fields and radiated powers by a moving electron

## Requirements
* ROOT 
* C++17 
* libRootFftwWrapper (available here https://github.com/nichol77/libRootFftwWrapper)
* Boost libraries

## Build instructions
The package is designed to be built with CMake

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

## Creating electron trajectories
The ```writeTrajectory``` executable provides the method to produce the electron trajectories. 
The outputted file is a ROOT file containing a TTree. 
Within this TTree there are branches corresponding to the time and electron dynamics. 
With this executable it is (currently) possible to configure:
* The output file
* The simulation time
* The simulation time step size
* The electron pitch angle (in degrees)
* The input file (for example, if you want to generate a trajectory using the final point of a previous file as the initial condition)

At present the trajectories are simulated in a 50cm long bathtub trap with a trapping field of 4.9mT. 
The background field is a uniform 1T in the Z direction
All electrons start at the middle of the trap. 
Options to configure the start points as well as the trap dimensions are to be added soon.

## Using the signal processing
Once an electron trajectory has been produced it is possible to generate the voltage induced on a given antenna using the ```InducedVoltage``` class which requires an input ROOT file and an antenna.
Full signal processing is done using the ```Signal``` class which takes as an input an instance of ```InducedVoltage```.
