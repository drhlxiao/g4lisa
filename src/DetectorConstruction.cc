/// detector description
// Author: Hualin Xiao(hualin.xiao@fhnw.ch)
// Date: Fri Jun 10, 2025
// stripTungstenEquivalent
#include "DetectorConstruction.hh"

#include <TFile.h>
#include <TTree.h>

#include <G4AssemblyVolume.hh>
#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Cons.hh>
#include <G4Element.hh>
#include <G4EllipticalTube.hh>
#include <G4GDMLParser.hh>
#include <G4GenericTrap.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4Polycone.hh>
#include <G4Polyhedra.hh>
#include <G4RotationMatrix.hh>
#include <G4RunManager.hh>
#include <G4SolidStore.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4Tubs.hh>
#include <G4TwoVector.hh>
#include <G4UImanager.hh>
#include <G4UnionSolid.hh>
#include <G4VisAttributes.hh>
#include <vector>

#include "AnalysisManager.hh"
#include "G4ExtrudedSolid.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
const G4double pi = CLHEP::pi;

DetectorConstruction::DetectorConstruction() {
	fWorldFile = "";
	checkOverlaps = true;

	detMsg = new DetectorMessenger(this);

}



DetectorConstruction::~DetectorConstruction() { delete detMsg; }

G4VPhysicalVolume *DetectorConstruction::Construct() {
	// AnalysisManager->SetAttenuatorStatus(attenuatorIn);

	G4NistManager *nist = G4NistManager::Instance();
	//Alum = nist->FindOrBuildMaterial("G4_Al");
	Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	Tungsten = nist->FindOrBuildMaterial("G4_W");

	//Platinum = nist->FindOrBuildMaterial("G4_Pt");
	//	Copper = nist->FindOrBuildMaterial("G4_Cu");

	//Siliver = nist->FindOrBuildMaterial("G4_Ag");
	//Gold = nist->FindOrBuildMaterial("G4_Au");
	//Nickle = nist->FindOrBuildMaterial("G4_Ni");

	G4Element *elAl = nist->FindOrBuildElement("Al");
	G4Element *elH = nist->FindOrBuildElement("H");
	G4Element *elC = nist->FindOrBuildElement("C");
	G4Element *elTi = nist->FindOrBuildElement("Ti");
	G4Element *elCd = nist->FindOrBuildElement("Cd");
	G4Element *elTe = nist->FindOrBuildElement("Te");
	G4Element *elAu = nist->FindOrBuildElement("Au");
	G4Element *elAg = nist->FindOrBuildElement("Ag");

	G4Element *elCr = nist->FindOrBuildElement("Cr");
	G4Element *elZn = nist->FindOrBuildElement("Zn");
	G4Element *elMn = nist->FindOrBuildElement("Mn");
	G4Element *elCu = nist->FindOrBuildElement("Cu");
	G4Element *elFe = nist->FindOrBuildElement("Fe");
	G4Element *elO = nist->FindOrBuildElement("O");
	G4Element *elSi = nist->FindOrBuildElement("Si");
	G4Element *elMg = nist->FindOrBuildElement("Mg");
	G4Element *elNi = nist->FindOrBuildElement("Ni");
	G4Element *elN = nist->FindOrBuildElement("N");
	G4Element *elW = nist->FindOrBuildElement("W");

    blackHole=new G4Material("blackHole",1e100*g/cm3, 1);
    blackHole->AddElement(elH,1);

	// construct world
	G4Box *worldSolid = new G4Box("worldSolid", 50 * cm, 50 * cm, 50 * cm);
	worldLogical =
		new G4LogicalVolume(worldSolid, Vacuum, "worldLogical", 0, 0, 0);
	worldPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), worldLogical,
			"worldPhysical", 0, false, 0);

	G4double plateHalfWidth =110*mm/2;
	G4double plateHalfDepth=15*mm;
	G4double pitchHalfWidth=115*um/2;

	G4double detectorHalfDepth=10*mm;

	G4Box *tungstenWindow= new G4Box("TungstenWindow", plateHalfWidth, plateHalfWidth, plateHalfDepth); //
	G4LogicalVolume *TungstenWindowLog=new G4LogicalVolume(tungstenWindow, Tungsten, "tungstenWindow", 0, 0, 0);

	new G4PVPlacement(0, G4ThreeVector(0,0,0), TungstenWindowLog, "TungstenPlate", worldLogical,
				false, 0, true);


	G4Box *tungstenPitch= new G4Box("TungstenPitch", pitchHalfWidth, plateHalfWidth, plateHalfDepth); //
	G4LogicalVolume *tungstenPitchLog=new G4LogicalVolume(tungstenPitch, Vacuum, "tungstenPitch", 0, 0, 0);

	G4double px = -plateHalfWidth +2*pitchHalfWidth;
	int i=0;

	while(px<plateHalfWidth-2*pitchHalfWidth)
	{
		new G4PVPlacement(0, G4ThreeVector(px, 0, 0), tungstenPitchLog,
				"tungstenPitch", TungstenWindowLog, false, i, true);
		px+= pitchHalfWidth *4;
		i++;
	}
	G4cout<<i<<" pitch created!"<<G4endl;

	G4Box *detectorBox= new G4Box("detectorBox", plateHalfWidth, plateHalfWidth, detectorHalfDepth); //
																									 //
	G4double detectorZ=41.5*cm;

	G4LogicalVolume *detectorLog=new G4LogicalVolume(detectorBox, blackHole, "detectorBox", 0, 0, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,detectorZ), detectorLog, "detector", worldLogical,
				false, 0, true);

	SetVisColors();
	worldLogical->SetVisAttributes(G4VisAttributes(false));
	G4cout << "World construction completed" << G4endl;
	return worldPhysical;
}



void DetectorConstruction::SetVisAttrib(G4LogicalVolume *log, G4double red,
		G4double green, G4double blue,
		G4double alpha, G4bool wireFrame,
		G4bool solid) {
	G4VisAttributes *visAttrib =
		new G4VisAttributes(G4Colour(red, green, blue, alpha));
	visAttrib->SetForceWireframe(wireFrame);
	visAttrib->SetForceSolid(solid);
	log->SetVisAttributes(visAttrib);
}
void DetectorConstruction::SetVisAttrib(G4LogicalVolume *log, G4double red,
		G4double green, G4double blue,
		G4double alpha) {
	SetVisAttrib(log, red, green, blue, alpha, true, true);
}


void DetectorConstruction::SetVisColors() {
	G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
	std::vector<G4LogicalVolume *>::const_iterator lvciter;
	for (lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++) {
		G4String volumeName = (*lvciter)->GetName();
		G4double red = G4UniformRand() * 0.7 + 0.15;
		G4double green = G4UniformRand() * 0.7 + 0.15;
		G4double blue = G4UniformRand() * 0.7 + 0.15;
		G4double alpha = 0.2;

		if (volumeName == "CdTeAnodeLog" || volumeName == "CdTeCathodeLog") {
			red = 1;
			green = 0;
			blue = 0;
			alpha = 0.05;
			SetVisAttrib(*lvciter, red, green, blue, alpha, true, false);
		} else if (volumeName.contains("frontGrid") ||
				volumeName.contains("rearGrid")) {
			(*lvciter)->SetVisAttributes(G4VisAttributes(false));
		} else if (volumeName.contains("grid")) {
			red = 0.28;
			green = 0.28;
			blue = 0.28;
			alpha = 0.1;
			SetVisAttrib(*lvciter, red, green, blue, alpha, true, true);
		} else if(volumeName=="CdTeLog"||volumeName=="calisteWorld"){
			red = 0;
			green = 0;
			blue = 1;
			alpha = 0.05;
			SetVisAttrib(*lvciter, red, green, blue, alpha, true, false);

		}else{
			SetVisAttrib(*lvciter, red, green, blue, alpha);
		}
		// randomize color

		double mass = (*lvciter)->GetMass() / g;
		G4cout << "~~~ The MASS of " << volumeName << " is " << mass << " g. ~~~"
			<< G4endl;
	}
}
