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
const G4ThreeVector singleDetectorPosition(0, 0, -0.8 * mm);
const G4double bigW = 2.15 * mm;
const G4double bigH = 4.55 * mm;
const G4double bigSW = 1.1 * mm;
const G4double bigSH = 0.455 * mm;
const G4double smallW = 1.05 * mm;
const G4double smallH = 0.86 * mm;
const G4double pixel0CenterX = -3.3 * mm;
const G4double pixel0CenterY = 2.3 * mm;
const G4double pixel4CenterX = -3.3 * mm;
const G4double pixel4CenterY = -2.3 * mm;
const G4double deltaW = 2.2 * mm;
const G4double pixel8CenterX = -3.85 * mm;
const G4double pixel8CenterY = 0 * mm;

const G4double cdteThickness = 1 * mm;
const G4double anodeThickness = (15 + 15 + 100) * 1e-9 * m;  // 130 nm
const G4double cathodeThickness = 15 * 1e-9 * m;             // 15 nm

const G4double calisteWidth = 11 * mm;
const G4double calisteLength = 12 * mm;
const G4double calisteBaseThickness = 14.4 * mm;
const G4double CdTeCalisteOffsetY =
0.8 * mm;  //  (top margion) 0.2 + 10 + 1.8 (Bottom margin)  = 12

const G4double TungstenGridDefaultThickness= 0.4*mm;

const G4double bondingLandZoneOnModuleLength = 0.4 * mm;
const G4double bondingLandZoneOnModuleWidth = 2 * mm;
const G4double bondingLandZoneOnModuleThickness = 0.25 * mm;
const G4double pinsThickness = 1.2 * mm;  //
G4double padCdTeBondingWidth = 0.4 * mm;
G4double padCdTeBondingLength = 1 * mm;
G4double padCdTeBondingThickness = 0.250 * mm;
G4double silverThickness = 0.2 * mm;
G4double platingThickness = (2 + 20 + 2.5 + 2) * 1e-3 * mm;
G4double padTotalThickness =
silverThickness + platingThickness;  // pads underneath the  CdTe
									 //
									 //


G4RotationMatrix *noRotation = new G4RotationMatrix(0., 0., 0.);
G4double CdTeTotalThickness = anodeThickness + cdteThickness + cathodeThickness;

const G4double calisteTotalHight =
calisteBaseThickness + cdteThickness + anodeThickness + cathodeThickness +
padTotalThickness + pinsThickness + padCdTeBondingThickness;
// 16.6

DetectorConstruction::DetectorConstruction() {
	// G4endl;
	fWorldFile = "";
	//tungstenGridThickness=TungstenGridDefaultThickness;
	checkOverlaps = true;

	detMsg = new DetectorMessenger(this);

}


G4LogicalVolume *DetectorConstruction::ConstructCdTe() {
	// cdte detector with pixels
	// the orgin is the top of surface

	G4Box *CdTeBox =
		new G4Box("CdTeBox", 5 * mm, 5 * mm, CdTeTotalThickness / 2.);
	G4LogicalVolume *CdTeLog =
		new G4LogicalVolume(CdTeBox, Vacuum, "CdTeLog", 0, 0, 0);
	// it doesn't matter what materials it is

	// thickness negligible  compared to the thickness  uncertainty
	checkOverlaps = true;
	////
	///// calisteResin base with different coating layers
	// caliste base material, is unknown
	// can be considered filled resin

	/////construct CdTe and pixels

	G4Box *CdTeSDBox = new G4Box("CdTeDetSD", 5 * mm, 5 * mm, cdteThickness / 2.);
	G4Box *CdTeAnodeBox =
		new G4Box("CdTeAnode", 5 * mm, 5 * mm, anodeThickness / 2.);
	G4Box *CdTeCathodeBox =
		new G4Box("CdTeCathode", 5 * mm, 5 * mm, cathodeThickness / 2.);

	// Top and Small pixels
	std::vector<G4TwoVector> vertexCoordsTop;
	vertexCoordsTop.push_back(G4TwoVector(-bigW / 2, bigH / 2));
	vertexCoordsTop.push_back(G4TwoVector(bigW / 2, bigH / 2));
	vertexCoordsTop.push_back(G4TwoVector(bigW / 2, -bigH / 2));
	vertexCoordsTop.push_back(G4TwoVector(-bigW / 2 + bigSW, -bigH / 2));
	vertexCoordsTop.push_back(G4TwoVector(-bigW / 2 + bigSW, -bigH / 2 + bigSH));
	vertexCoordsTop.push_back(G4TwoVector(-bigW / 2, -bigH / 2 + bigSH));

	G4ExtrudedSolid *bigPixelTopGeo =
		new G4ExtrudedSolid("topBigPixel", vertexCoordsTop, cdteThickness / 2.,
				G4TwoVector(), 1, G4TwoVector(), 1);
	//top row
	std::vector<G4TwoVector> vertexCoordsBott;
	vertexCoordsBott.push_back(G4TwoVector(-bigW / 2, bigH / 2 - bigSH));
	vertexCoordsBott.push_back(G4TwoVector(-bigW / 2 + bigSW, bigH / 2 - bigSH));
	vertexCoordsBott.push_back(G4TwoVector(-bigW / 2 + bigSW, bigH / 2));
	vertexCoordsBott.push_back(G4TwoVector(bigW / 2, bigH / 2));
	vertexCoordsBott.push_back(G4TwoVector(bigW / 2, -bigH / 2));
	vertexCoordsBott.push_back(G4TwoVector(-bigW / 2, -bigH / 2));

	G4ExtrudedSolid *bigPixelBottomGeo =
		new G4ExtrudedSolid("bottBigPixel", vertexCoordsBott, cdteThickness / 2.,
				G4TwoVector(), 1, G4TwoVector(), 1);
	//bottom row
	G4Box *smallPixelGeo =
		new G4Box("smallPixel", smallW / 2, smallH / 2, cdteThickness / 2.);


	G4Material *cdtesdmat=CdTe;
	//set to vaccum
	G4cout<<">>>CdTe Coating included!"<<G4endl;

	G4LogicalVolume *CdTeDetSDLog =
		new G4LogicalVolume(CdTeSDBox, cdtesdmat, "CdTeDetSDLog", 0, 0, 0);
	//detector, used to place senstive volumes






	G4LogicalVolume *bigPixelTopLog =
		new G4LogicalVolume(bigPixelTopGeo, CdTe, "bigPixelTopLog", 0, 0, 0);
	G4LogicalVolume *bigPixelBottomLog = new G4LogicalVolume(
			bigPixelBottomGeo, CdTe, "bigPixelBottomLog", 0, 0, 0);
	G4LogicalVolume *smallPixelLog =
		new G4LogicalVolume(smallPixelGeo, CdTe, "smallPixelLog", 0, 0, 0);

	G4int copyNb;
	G4double detZ = 0;
	// place electrode
	G4String name = "pixel";
	for (int i = 0; i < 4; i++) {
		// it becomes pixel4 after rotation
		copyNb = i;
		G4cout<<"Adding pixels:"<<copyNb<<G4endl;
		G4ThreeVector posBigPixelTop(pixel0CenterX + deltaW * i, pixel0CenterY,
				detZ);
		new G4PVPlacement(0, posBigPixelTop, bigPixelTopLog, name, CdTeDetSDLog,
				false, copyNb, checkOverlaps);

		// it becomes pixel4 after rotation
		G4cout<<"Adding pixels:"<<copyNb<<G4endl;
		copyNb = i + 4;
		G4ThreeVector posBigPixelBottom(pixel4CenterX + deltaW * i, pixel4CenterY,
				detZ);
		new G4PVPlacement(0, posBigPixelBottom, bigPixelBottomLog, name,
				CdTeDetSDLog, false, copyNb, checkOverlaps);

		// small pixels
		copyNb = i + 8;
		G4ThreeVector posSmallPixel(pixel8CenterX + deltaW * i, pixel8CenterY,
				detZ);

		G4cout<<"Adding pixels:"<<copyNb<<G4endl;

		new G4PVPlacement(0, posSmallPixel, smallPixelLog, name, CdTeDetSDLog,
				false, copyNb, checkOverlaps);
	}


	G4ThreeVector TmDet(
			0, 0, CdTeTotalThickness / 2. - cathodeThickness - cdteThickness / 2);
	new G4PVPlacement(0, TmDet, CdTeDetSDLog, "CdTeDetSD", CdTeLog, false, 0,
			checkOverlaps);

	G4cout<<">>>Constructing Full cdte , constructing coating, anode, cathod"<<G4endl;


	G4ThreeVector TmCdTeCathod(0, 0,
			CdTeTotalThickness / 2. - cathodeThickness / 2);
	G4ThreeVector TmAnode(0, 0,
			CdTeTotalThickness / 2. - cathodeThickness -
			cdteThickness - anodeThickness / 2);

	//useful to record incoming photons, kill tracks in CdTe
	G4cout<<"Adding CdTe gold layer..."<<G4endl;
	G4LogicalVolume *CdTeAnodeLog = new G4LogicalVolume(
			CdTeAnodeBox, goldLayerMaterial, "CdTeAnodeLog", 0, 0, 0);

	G4cout<<"Adding CdTe container..."<<G4endl;
	G4LogicalVolume *CdTeCathodeLog =
		new G4LogicalVolume(CdTeCathodeBox, Platinum, "CdTeCathodeLog", 0, 0, 0);

	G4cout<<"Adding CdTe cathod..."<<G4endl;
	new G4PVPlacement(0, TmCdTeCathod, CdTeCathodeLog, "CdTeCathodPhys", CdTeLog,
			false, 0, checkOverlaps);
	G4cout<<"Adding CdTe anode..."<<G4endl;
	new G4PVPlacement(0, TmAnode, CdTeAnodeLog, "CdTeAnodePhys", CdTeLog, false,
			0, checkOverlaps);

	return CdTeLog;
	// done
}
G4LogicalVolume *DetectorConstruction::ConstructCalisteBase() {
	G4cout<<"Adding CdTe base..."<<G4endl;

	G4Box *calisteBaseOuter =
		new G4Box("CdTeModuleBaseOuter", calisteWidth / 2, calisteLength / 2,
				calisteBaseThickness / 2);

	G4LogicalVolume *calisteBaseOuterLog =
		new G4LogicalVolume(calisteBaseOuter, Gold, "calisteBaseOuter", 0, 0, 0);

	// Au+Ni+Cu+2u
	////////Plating from inside to outside:
	//+ All surfaces
	//+ Ni 2u + Cu 20u + Ni 2.5u + Au 2 u
	//
	G4double coatingThickness[4] = {2e-3 * mm, 2.5e-3 * mm, 20e-3 * mm,
		2e-3 * mm};
	// we need to reverse the order
	G4Material *coatingMateris[4] = {Nickle, Copper, Nickle, Resin};
	// outest layer is gold
	//  different layer of the coating material is  considered as mixture

	G4Box *calisteBaseInner[4];
	G4String calisteBaseInnerLayerName[4] = {"nickle", "copper", "nickle",
		"resin"};

	G4LogicalVolume *calisteBaseInnerLog[4];
	G4double totalCoatingThickness = 0;

	G4LogicalVolume *coatingMotherLog;

	for (int i = 0; i < 4; i++) {
		totalCoatingThickness += coatingThickness[i];
		calisteBaseInner[i] = new G4Box(
				G4String("CdTeModuleBaseInner_") + calisteBaseInnerLayerName[i],
				calisteWidth / 2 - totalCoatingThickness,
				calisteLength / 2 - totalCoatingThickness, calisteBaseThickness / 2);
		calisteBaseInnerLog[i] = new G4LogicalVolume(
				calisteBaseInner[i], coatingMateris[i],
				G4String("log_CalisteBaseInner_") + calisteBaseInnerLayerName[i], 0, 0,
				0);
		if (i == 0) {
			coatingMotherLog = calisteBaseOuterLog;
		} else {
			coatingMotherLog = calisteBaseInnerLog[i - 1];
		}

		new G4PVPlacement(0, G4ThreeVector(0, 0, 0), calisteBaseInnerLog[i],
				calisteBaseInnerLayerName[i] + "_phys", coatingMotherLog,
				false, 0, checkOverlaps);
	}
	return calisteBaseOuterLog;
}
G4AssemblyVolume *DetectorConstruction::ConstructPads() {
	//////construct pads
	//pads are placed behine cdTe
	//////
	G4cout<<"Adding electrode pads..."<<G4endl;
	G4AssemblyVolume *padAssembly = new G4AssemblyVolume();
	G4double padZ = 0 * mm;

	G4ThreeVector TmStackInPadMother(
			0, 0, -padTotalThickness / 2 + platingThickness / 2);

	// round pad
	G4Tubs *roundPadGeo = new G4Tubs("roundPadGeo", 0, 0.3 * mm,
			padTotalThickness / 2., 0, 360. * deg);
	G4Tubs *roundPadStackGeo = new G4Tubs("roundPadStackGeo", 0, 0.3 * mm,
			platingThickness / 2., 0, 360. * deg);

	G4LogicalVolume *roundPadLog =
		new G4LogicalVolume(roundPadGeo, SilverEpoxy, "roundPadLog", 0, 0, 0);
	G4LogicalVolume *roundPadStackLog = new G4LogicalVolume(
			roundPadStackGeo, padStackMaterial, "roundPadStackLog", 0, 0, 0);
	new G4PVPlacement(0, TmStackInPadMother, roundPadStackLog,
			"roundPadStackPhys", roundPadLog, false, 0, checkOverlaps);

	// elliptical pads
	G4EllipticalTube *ellipticalPadGeo = new G4EllipticalTube(
			"ellipticalGeo", 0.3, 0.4 * mm, padTotalThickness / 2);
	G4EllipticalTube *ellipticalPadStackGeo = new G4EllipticalTube(
			"ellipticalStackGeo", 0.3, 0.4 * mm, platingThickness / 2);
	G4LogicalVolume *ellipticalPadLog = new G4LogicalVolume(
			ellipticalPadGeo, SilverEpoxy, "ellipticalPadLog", 0, 0, 0);
	G4LogicalVolume *ellipticalPadStackLog =
		new G4LogicalVolume(ellipticalPadStackGeo, padStackMaterial,
				"ellipticalPadStackLog", 0, 0, 0);
	new G4PVPlacement(0, TmStackInPadMother, ellipticalPadStackLog,
			"ellipticalPadStackLog", ellipticalPadLog, false, 0,
			checkOverlaps);

	G4Box *rectPadGeo =
		new G4Box("rectPadGeo", 0.3 * mm, 4.45 / 2 * mm, padTotalThickness / 2.);
	G4LogicalVolume *rectPadLog =
		new G4LogicalVolume(rectPadGeo, SilverEpoxy, "rectPadLog", 0, 0, 0);

	G4Box *rectPadStackGeo = new G4Box("rectPadStckGeo", 0.3 * mm, 4.45 / 2 * mm,
			platingThickness / 2.);
	G4LogicalVolume *rectPadStackLog = new G4LogicalVolume(
			rectPadStackGeo, padStackMaterial, "rectPadStackLog", 0, 0, 0);
	new G4PVPlacement(0, TmStackInPadMother, rectPadStackLog, "rectPadStackPhys",
			rectPadLog, false, 0, checkOverlaps);

	G4double roundPadX[] = {
		-3.744, -2.832, -1.541, -0.646, 0.714,  1.489,  2.849,  3.744,  -3.727,
		-2.832, -1.541, -0.628, 0.680,  1.558,  2.832,  3.744,  -3.744, -2.849,
		-1.523, -0.646, 0.646,  1.541,  2.832,  3.727,  -3.865, 2.746,  -3.710,
		-2.832, -1.558, -0.628, 1.541,  2.832,  3.744,  3.727,  2.815,  1.523,
		0.611,  -0.663, -1.523, -2.832, -3.710, 3.744,  2.849,  1.523,  0.628,
		-0.663, -1.523, -2.832, -3.744, 0.559,  -1.627, 0.628};
	G4double roundPadY[] = {
		-4.678, -4.678, -4.695, -4.678, -4.678, -4.695, -4.644, -4.661, -3.768,
		-3.768, -3.785, -3.768, -3.768, -3.751, -3.768, -3.768, -2.858, -2.876,
		-2.910, -2.893, -2.893, -2.876, -2.876, -2.893, -0.747, -0.798, 1.296,
		1.313,  1.313,  1.313,  1.348,  1.348,  1.365,  2.258,  2.240,  2.240,
		2.258,  2.223,  2.240,  2.240,  2.240,  3.150,  3.150,  3.150,  3.150,
		3.133,  3.116,  3.133,  3.133,  -0.764, -0.764, 1.365};
	G4double ellipitalPadX[] = {-2.746, -2.763, -0.525, -0.542,
		1.644,  1.644,  3.830,  3.847};
	G4double ellipitalPadY[] = {-1.433, -0.060, -1.451, -0.077,
		-1.451, -0.060, -1.433, -0.077};
	G4double rectPadX[] = {4.71, -4.75};
	G4double rectPadY[] = {-0.755, -0.755};
	// positions extracted from the graph and verified using pyplot
	//	cdteThickness - anodeThickness - padTotalThickness/2.;

	// placing pads
	//
	for (int i = 0; i < 52; i++) {
		G4ThreeVector Tm(roundPadX[i] * mm, -roundPadY[i] * mm, padZ);
		padAssembly->AddPlacedVolume(roundPadLog, Tm, noRotation);
	}
	for (int i = 0; i < 8; i++) {
		G4ThreeVector Tm(ellipitalPadX[i] * mm, -ellipitalPadY[i] * mm, padZ);
		padAssembly->AddPlacedVolume(ellipticalPadLog, Tm, noRotation);
	}
	for (int i = 0; i < 2; i++) {
		G4ThreeVector Tm(rectPadX[i] * mm, -rectPadY[i] * mm, padZ);
		padAssembly->AddPlacedVolume(rectPadLog, Tm, noRotation);
	}
	// y reversed
	//  see email from olivier, sent on Feb. 09 2023: CAD model Caliste-SO
	//
	// Place the gold bonding rod for HV
	G4Box *bondingLandZoneOnModuleGeo = new G4Box(
			"bondingLandZoneOnModuleGeo", bondingLandZoneOnModuleWidth / 2,
			bondingLandZoneOnModuleLength / 2, bondingLandZoneOnModuleThickness / 2);

	G4LogicalVolume *bondingLandZoneOnModuleLog = new G4LogicalVolume(
			bondingLandZoneOnModuleGeo, Gold, "bondingLandZoneOnModuleLog", 0, 0, 0);
	G4Box *padCdTeBondingGeo =
		new G4Box("padCdTeBondingGeo", padCdTeBondingWidth / 2,
				padCdTeBondingLength / 2, padCdTeBondingThickness / 2);
	G4LogicalVolume *padCdTeBondingLog = new G4LogicalVolume(
			padCdTeBondingGeo, Gold, "padCdTeBondingLog", 0, 0, 0);
	G4ThreeVector TmHVPad(
			calisteWidth / 2 - bondingLandZoneOnModuleWidth / 2,
			-calisteLength / 2. + bondingLandZoneOnModuleLength / 2.,
			padZ + (bondingLandZoneOnModuleThickness - padTotalThickness) / 2.);

	padAssembly->AddPlacedVolume(bondingLandZoneOnModuleLog, TmHVPad, noRotation);
	// placed at the corner corner
	G4ThreeVector TmCdTePadInAssembly(4.8 * mm, -3.5 * mm,
			padZ + padTotalThickness / 2. +
			CdTeTotalThickness +
			padCdTeBondingThickness / 2);
	padAssembly->AddPlacedVolume(padCdTeBondingLog, TmCdTePadInAssembly,
			noRotation);
	// center is at the pad center
	return padAssembly;
}

G4LogicalVolume *DetectorConstruction::ConstructCaliste() {
	G4Box *calisteWorld = new G4Box("CdTeModuleWorld", calisteWidth / 2,
			calisteLength / 2, calisteTotalHight / 2);
	G4LogicalVolume *calisteLog;

	calisteLog= new G4LogicalVolume(calisteWorld, Vacuum, "calisteWorld", 0, 0, 0);
	//construct cdte with coating
	G4cout<<"Constructing CdTe detector..."<<G4endl;
	G4LogicalVolume *CdTeLog = ConstructCdTe();


	G4ThreeVector TmCdTe(
			0, CdTeCalisteOffsetY,
			calisteTotalHight / 2 - padCdTeBondingThickness - CdTeTotalThickness / 2);
	// shifted by 0.8 mm
	new G4PVPlacement(0, TmCdTe, CdTeLog, "CdTeLog", calisteLog, false, 0,
			checkOverlaps);

	//if(isSimpleCaliste){
	//	G4cout<<">>>>Simple CdTe detector constructed, only 8 big pixels, no coating pads, or base"<<G4endl;
	//	return calisteLog;
	//}
	//don't contruct more for simple detector

	//place pads
	G4ThreeVector TmHVBondingPad(
			0, 0,
			calisteTotalHight / 2 - (padCdTeBondingThickness + CdTeTotalThickness) -
			padTotalThickness / 2);

	G4cout<<"Constructing CdTe pads..."<<G4endl;
	G4AssemblyVolume *padAssembly = ConstructPads();
	padAssembly->MakeImprint(calisteLog, TmHVBondingPad, noRotation);



	//create base
	G4ThreeVector TmCalisteBase(0, 0,
			calisteTotalHight / 2 - padCdTeBondingThickness -
			CdTeTotalThickness - padTotalThickness -
			calisteBaseThickness / 2);
	G4cout<<"Constructing CdTe base..."<<G4endl;
	G4LogicalVolume *calisteBaseLog = ConstructCalisteBase();
	new G4PVPlacement(0, TmCalisteBase, calisteBaseLog, "calisteBasePhys",
			calisteLog, false, 0, checkOverlaps);




	return calisteLog;
}

DetectorConstruction::~DetectorConstruction() { delete detMsg; }

G4VPhysicalVolume *DetectorConstruction::Construct() {
	// AnalysisManager->SetAttenuatorStatus(attenuatorIn);

	G4NistManager *nist = G4NistManager::Instance();
	Alum = nist->FindOrBuildMaterial("G4_Al");
	Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	Tungsten = nist->FindOrBuildMaterial("G4_W");

	Platinum = nist->FindOrBuildMaterial("G4_Pt");
	Copper = nist->FindOrBuildMaterial("G4_Cu");

	Siliver = nist->FindOrBuildMaterial("G4_Ag");
	Gold = nist->FindOrBuildMaterial("G4_Au");
	Nickle = nist->FindOrBuildMaterial("G4_Ni");

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

	G4double density = 5.85 * g / cm3;  // STIX-DS-0017-PSI
	G4int nelements, natoms;
	G4double fractionmass;
	CdTe = new G4Material("CdTe", density, nelements = 2);
	CdTe->AddElement(elCd, natoms = 1);
	CdTe->AddElement(elTe, natoms = 1);

	goldLayerMaterial =
		new G4Material("goldLayerMaterial", density, nelements = 3);
	goldLayerMaterial->AddElement(elAu, 0.94700687);  // mass fraction
	goldLayerMaterial->AddElement(elTi, 0.033120707);
	goldLayerMaterial->AddElement(elAl, 0.019872424);

	// Same as the body: Ni 2u + Cu 20u + Ni 2.5u + Au 2 u

	density = 9.0 * g / cm3;
	padStackMaterial = new G4Material("padStackMaterial", density, nelements = 3);
	padStackMaterial->AddElement(elCu, 0.695);  // mass fraction
	padStackMaterial->AddElement(elAu, 0.15);
	padStackMaterial->AddElement(elNi, 0.155);
	// it is considered as a mixture, maybe it is OK for X rays

	SiO2 = new G4Material("SiO2", density = 2.2 * g / cm3, nelements = 2);
	SiO2->AddElement(elSi, 1);
	SiO2->AddElement(elO, 2);

	density=19.3 * g/cm3;
	//gridDensity=density * 0.4 * mm/tungstenGridThickness;
	//used for study of grid shadow, the thickness to reduced to compare the thickness


	blackHole= new G4Material("blackHole", density = 1e100* g / cm3, nelements = 1);
	blackHole->AddElement(elH, 1);
	//

	Alum7075 = new G4Material("Alum7075", density = 2.8 * g / cm3, nelements = 7);
	Alum7075->AddElement(elAl, 0.876);
	Alum7075->AddElement(elZn, 0.056);
	Alum7075->AddElement(elCr, 0.023);
	Alum7075->AddElement(elMg, 0.025);
	Alum7075->AddElement(elCu, 0.016);
	Alum7075->AddElement(elFe, 0.0025);
	Alum7075->AddElement(elMn, 0.0015);

	density = 1.43* g / cm3;
	Kapton= new G4Material("Kapton", density, nelements = 4);
	Kapton->AddElement(elH, 0.026362);
	Kapton->AddElement(elC, 0.691133);
	Kapton->AddElement(elN, 0.07327);
	Kapton->AddElement(elO, 0.209235);

	density = 1.2 * g / cm3;

	Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, nelements = 2);
	Epoxy->AddElement(elH, natoms = 2);
	Epoxy->AddElement(elC, natoms = 2);



	density = 1.86 * g / cm3;
	FR4 = new G4Material("FR4", density, nelements = 2);
	FR4->AddMaterial(SiO2, 0.528);
	FR4->AddMaterial(Epoxy, 0.472);

	//x-ray window with solar black 

    // Define elements
    G4Element* elCa = nist->FindOrBuildElement("Ca");
    G4Element* elP = nist->FindOrBuildElement("P");
    G4Element* elBe = nist->FindOrBuildElement("Be");
    //G4Element* elFe = nist->FindOrBuildElement("Fe");
    //G4Element* elSi = nist->FindOrBuildElement("Si");

    // Define isotopes if necessary
    // G4Isotope* isoX = new G4Isotope("name", Z, A, mass);

    // Define the material
    berylliumMat= new G4Material("beryllium_with_solar_black", 1.84339 * g/cm3, 9); //STIX X-ray window material
    berylliumMat->AddElement(elC, 0.00877056277056277);
    berylliumMat->AddElement(elCa, 0.0021774891774891773);
    berylliumMat->AddElement(elP, 0.0008441558441558443);
    berylliumMat->AddElement(elBe, 0.9807359307359307);
    berylliumMat->AddElement(elO, 0.005525974025974026);
    berylliumMat->AddElement(elAl, 0.0004978354978354978);
    berylliumMat->AddElement(elFe, 0.0007467532467532468);
    berylliumMat->AddElement(elMg, 0.0003982683982683983);
    berylliumMat->AddElement(elSi, 0.0002987012987012987);

	// 2/ Molding resin
	// Use the following mixture
	//+ density of the mixture 1.77
	
G4Material* aluminium6061Mat = new G4Material("aluminium6061", 2.70 * g/cm3, 8);
    aluminium6061Mat->AddElement(elAl, 0.9738); // Adjusted to 97.38%
    aluminium6061Mat->AddElement(elMg, 0.0100); // 1.00%
    aluminium6061Mat->AddElement(elSi, 0.0060); // 0.60%
    aluminium6061Mat->AddElement(elFe, 0.0035); // 0.35%
    aluminium6061Mat->AddElement(elCu, 0.00275); // 0.275%
    aluminium6061Mat->AddElement(elCr, 0.00195); //
    aluminium6061Mat->AddElement(elZn, 0.00125); //
    aluminium6061Mat->AddElement(elTi, 0.00075); //
	//geometry: 

	// 4.1 mm  4.2 mm Al ===> 22 mm X-ray window

	density = 1.77 * g / cm3;
	Resin = new G4Material("Resin", density, nelements = 2);
	Resin->AddMaterial(SiO2, fractionmass = 0.27);
	Resin->AddMaterial(Epoxy, fractionmass = 0.73);

	SilverEpoxy =
		new G4Material("SilverEpoxy", density = 2.566 * g / cm3, nelements = 2);
	// not 3.3 g/cm3, should be 2.566 measured by Olivier
	SilverEpoxy->AddMaterial(Epoxy, fractionmass = 0.7);
	SilverEpoxy->AddMaterial(Siliver, fractionmass = 0.3);

	LeadPadMat =
		new G4Material("LeadPadMat", density = 11 * g / cm3, nelements = 2);
	LeadPadMat->AddMaterial(Nickle, fractionmass = 0.587);
	LeadPadMat->AddMaterial(Gold, fractionmass = 0.413);

	// construct world
	G4Box *worldSolid = new G4Box("worldSolid", 50 * cm, 50 * cm, 50 * cm);
	worldLogical =
		new G4LogicalVolume(worldSolid, Vacuum, "worldLogical", 0, 0, 0);
	worldPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), worldLogical,
			"worldPhysical", 0, false, 0);

	/*
	 * The following code is for PARDRE 
	//construct  Beryllium  window
	G4Box *berylliumWindow= new G4Box("berylliumWindow", 10 * mm, 10 * mm, 0.5 * mm); //X-ray window dimensions:  20 x 20 x 1 mm^3
	berylliumWindowLog =new G4LogicalVolume(berylliumWindow, berylliumMat, "berylliumWindowLog", 0, 0, 0);
	//construct Aluminium  Window
	//
	G4Box *alumWindow= new G4Box("alumWindow", 10 * mm, 10 * mm, 0.2 * mm); //X-ray window dimensions:  20 x 20 x 1 mm^3
	alumWindowLog =new G4LogicalVolume(alumWindow, aluminium6061Mat, "alumWindowLog", 0, 0, 0);
	new G4PVPlacement(0, G4ThreeVector(0, 0, 4.6*mm), berylliumWindowLog, 
			"berylliumWindowPhys", worldLogical,  false, 0, true);
	//distance to the CdTe surface  = 4.1 + 0.5
	new G4PVPlacement(0, G4ThreeVector(0, 0, 20.2*mm), alumWindowLog,
			"alumWindowPhys", worldLogical, false, 0, true);
		 //distance to the CdTe surface  = 20 + 0.2
	G4RotationMatrix rotY;
	//rotY.rotateY(-90 * deg);
	rotY.rotateY(0 * deg);
	G4RotationMatrix rotX;
	//rotX.rotateX(90 * deg);
	rotX.rotateX(0 * deg);
	rotMatrix = rotX * rotY;

	G4LogicalVolume *CalisteLog = ConstructCaliste();
	G4cout << "Placing detector Caliste..." << G4endl;
	G4ThreeVector pos(0,0, -calisteTotalHight/2); // the origin is at the surface of CadTe
	new G4PVPlacement(G4Transform3D(rotMatrix, pos), CalisteLog,"Caliste",
			worldLogical, false, 0, true);
	*/
	G4double plateHalfWidth =50*mm;
	G4double plateHalfDepth=15*mm;
	G4double pitchHalfWidth=0.05*mm;

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
