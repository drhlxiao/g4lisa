
# g4padre
A Geant4 Monte Carlo Simulation Package for PADRE (The Solar Polarization and Directivity X-Ray Experiment), a groundbreaking NASA CubeSat mission.

PADRE website: https://padre.ssl.berkeley.edu/

It models a basic detector setup consisting of:

* A single **CdTe detector**
* **X-ray windows**
* An **aluminum window**

This setup is intended for x-ray detection studies or as a clean starting point for detector development and testing in simulation.

This simulation package is based on the Monte Carlo Simulation package for STIX (Solar X-ray Imager/telescope onboard Solar Orbiter), 
     which was developed by Hualin.

## Project Structure

* Based on the `g4stix` package (simplified)
* Minimal geometry focused on essential components only

## Author

**Hualin Xiao**
[hualin.xiao@fhnw.ch](mailto:hualin.xiao@fhnw.ch)

## History

* **Jun 11, 2025**: Initial version extracted and simplified from `g4stix`

## Requirements

* Geant4 (tested with version ≥10.7)
* CMake
* ROOT



## Workflow

* compile source code:
  - cd <SOURCE_CODE>
  - cmake .
  - make

* validation with visualization
  - ./g4main --gui -m vis.mac
* run simulations
  - ./g4main  -m response.mac -o response.root
* create response matrix from the simulation outputs
  - cd analysis
  - python process_root.py



