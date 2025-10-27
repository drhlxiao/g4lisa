


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



