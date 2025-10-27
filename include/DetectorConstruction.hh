

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

// STL //
#include <string>
#include "DetectorMessenger.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4AssemblyVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume *Construct();

  void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green,
                    G4double blue, G4double alpha, G4bool wireFrame,
                    G4bool solid);
  void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green,
                    G4double blue, G4double alpha);
  void SetVisColors();


private:

  G4String fWorldFile;
  G4LogicalVolume *worldLogical;
  G4Material *Tungsten,*Vacuum, *blackHole;
  bool checkOverlaps;
  G4RotationMatrix rotMatrix;

  G4VPhysicalVolume *worldPhysical;
  DetectorMessenger *detMsg;
};

#endif
