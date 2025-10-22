//
/// \file /include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class DetectorMessenger : public G4UImessenger {
public:
  DetectorMessenger(DetectorConstruction *);
  ~DetectorMessenger();

  void SetNewValue(G4UIcommand *, G4String);

private:
  DetectorConstruction *fDetector;
  // G4UIdirectory *fDetectorDir;
  G4UIcmdWithABool *fSetAttenStatusCmd;
  G4UIcmdWithABool *fSetGridStatusCmd;
  G4UIcmdWithABool *fSetDetectorStatusCmd, *fSetDetectorSimpleCmd;
  G4UIcmdWithAString *fSetGdmlCmd;
  G4UIcmdWithAnInteger *fDetectorSelectionCmd,*fSetGridMaskCmd;
  G4UIcmdWithADoubleAndUnit *fSetGridThicknessCmd;
  // G4UIcmdWithAString *fSetCADTypeCommand;
};

#endif
