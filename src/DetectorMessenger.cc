/***************************************************************
 * Class to define detector messeger.
 * Author  : Hualin Xiao
 * Date    : Jan, 2015
 * Version : 1.10
 *
 ***************************************************************/

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction *theDet)
    : fDetector(theDet) {


}

DetectorMessenger::~DetectorMessenger() {
}

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
}
