/***************************************************************
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/
#include "AnalysisManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4UnitsTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNamed.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"

bool DEBUG = false;
const int MAX_NUM_TREE_TO_FILL = 1000000;

AnalysisManager *AnalysisManager::fManager = 0;
AnalysisManager *AnalysisManager::GetInstance()
{
	if (!fManager)
	{
		fManager = new AnalysisManager();
	}
	return fManager;
}
AnalysisManager::AnalysisManager()
{

	// isDummyDetector=false;
	numKilled = 0;
	isNewEvent = true;
	numEventIn = 0;
	numEventOut = 0;
	killTracksEnteringGrids = false;
	killTracksEnteringDetectors = false;
	numInpTreeFilled = 0;
	numSourceTreeFilled = 0;

	numPhysTreeFilled = 0;
}
void AnalysisManager::CopyMacrosToROOT(TFile *f, TString &macfilename)
{
	G4cout << "Copying macros from file " << macfilename << " to root file"
		   << G4endl;
	if (macfilename == "" || !f)
		return;
	std::ifstream infile(macfilename.Data());
	if (!infile.good())
	{
		G4cout << "can not open the macro file, existing..." << G4endl;
		return;
	}
	G4String macros = Form("\nCommand:%s\n", commandLine.data());
	macros += Form("\nMacro filename: %s\n ", macroFilename.Data());
	std::string line;
	G4cout << "Macros read from the macro file:" << G4endl;
	while (std::getline(infile, line))
	{
		macros += line + "\n";
		G4cout << line << G4endl;
	}
	TNamed cmd;
	cmd.SetTitle(macros);
	f->cd();
	cmd.Write("metadata");
	infile.close();
}

void AnalysisManager::InitRun(const G4Run *run)
{
	rootFile = new TFile(outputFilename.Data(), "recreate");
	
	inpTree = new TTree("inp", "inp");
	inpTree->Branch("pos", inpPos, "pos[3]/D");
	inpTree->Branch("E0", &gunEnergy, "E0/D");
	inpTree->Branch("eventID", &eventID, "eventID/I");
	inpTree->Branch("itrack", &itrack, "itrack/I");
	inpTree->Branch("boundary", &boundary, "boundary/I");
	inpTree->Branch("pixelID", &pixelID, "pixelID/I");
	inpTree->Branch("detectorID", &detectorID, "detectorID/I");
	inpTree->Branch("v", inpVec, "v[3]/D");
	inpTree->Branch("theta", &inpTheta, "theta/D");
	inpTree->Branch("energy", &inpEnergy, "energy/D");
	inpTree->Branch("pdg", &inpPDG, "pdg/I");
	inpTree->Branch("parent", &parentID, "parent/I");

	primTree = new TTree("source", "source");
	primTree->Branch("pos", gunPosition, "pos[3]/D");
	primTree->Branch("eventID", &eventID, "eventID/I");
	primTree->Branch("vec", gunDirection, "vec[3]/D");
	primTree->Branch("E0", &gunEnergy, "E0/D");

	h2xy = new TH2F("h2xy", "Locations of hits; X (mm); Y(mm)", 1800, -90, 90,
					1800, -90, 90);
}

//////////////////////////////////////////////////////////////////////////

void AnalysisManager::ProcessRun(const G4Run *run)
{
	CloseROOT();
	G4cout << "Events entered detectors:" << numEventIn << G4endl;
	G4cout << "Events escaped from detectors:" << numEventOut << G4endl;
	G4cout << "Events captured:" << numEventIn - numEventOut << G4endl;
}

/// EventAction
void AnalysisManager::InitEvent(const G4Event *event)
{


	itrack = 0;
	totalNumSteps = 0;
	isNewEvent = true;
}
void AnalysisManager::ProcessEvent(const G4Event *event)
{
	eventID = event->GetEventID();
	G4bool effectiveEvent = false;
	

	if (numSourceTreeFilled < 100000)
	{
		primTree->Fill();
		numSourceTreeFilled++;
	}
	// if (effectiveEvent)  evtTree->Fill();
	isNewEvent = true;
}
// Stepping Action
void AnalysisManager::UpdateParticleGunInfo()
{
	// write particle information to the tree
	if (!isNewEvent)
		return;
	// don't update

	G4RunManager *runManager = G4RunManager::GetRunManager();
	PrimaryGeneratorAction *primaryAction =
		(PrimaryGeneratorAction *)runManager->GetUserPrimaryGeneratorAction();
	G4ThreeVector position, direction;
	primaryAction->GetGPS(position, direction, gunEnergy);

	gunPosition[0] = position.getX() / mm;
	gunPosition[1] = position.getY() / mm;
	gunPosition[2] = position.getZ() / mm;
	gunDirection[0] = direction.getX();
	gunDirection[1] = direction.getY();
	gunDirection[2] = direction.getZ();
	isNewEvent = false;
}

////////////////////////////////////////////////////////////////////

inline void AnalysisManager::AddEnergy(G4int detId, G4double edep)
{
}

AnalysisManager::~AnalysisManager()
{
	if (fManager)
		delete fManager;
	// if (rootFile) delete rootFile;
	fManager = 0;
}

void AnalysisManager::CloseROOT()
{
	rootFile->cd();
	inpTree->Write();
	primTree->Write();
	CopyMacrosToROOT(rootFile, macroFilename);

	rootFile->Close();
}

void AnalysisManager::ProcessStep(const G4Step *aStep)
{

	UpdateParticleGunInfo();

	G4double px, py, pz;
	G4double edep;
	G4String volName;
	const G4Track *track = aStep->GetTrack();

	if (track->GetVolume())
		volName = track->GetVolume()->GetName();

	if (volName == "detector")
	{ 
		// don't store too many tracks
		FillDetectorIncidentParticle(aStep);
		aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		numKilled++;
	}
	// dummy detector, material as black hole
	// used to study grid effect
}
void AnalysisManager::FillDetectorIncidentParticle(const G4Step *aStep)
{

	G4StepPoint *preStep = aStep->GetPreStepPoint();
	G4StepPoint *postStep = aStep->GetPostStepPoint();
	inpEnergy = preStep->GetKineticEnergy() / keV;
	const G4Track *track = aStep->GetTrack();
	G4ThreeVector prePos = postStep->GetPosition();

	G4double px, py, pz;
	// G4cout<<"filling2 "<<G4endl;
	px = prePos.x() / mm;
	py = prePos.y() / mm;
	pz = prePos.z() / mm;
	h2xy->Fill(py , pz );

	// if (numInpTreeFilled < MAX_NUM_TREE_TO_FILL) {
	G4ThreeVector inpV = track->GetMomentumDirection();
	inpVec[0] = inpV.x();
	inpVec[1] = inpV.y();
	inpVec[2] = inpV.z();
	inpTheta = asin(sqrt(inpVec[1] * inpVec[1] + inpVec[2] * inpVec[2])) * 180 /
			   3.1415926;

	inpPos[0] = px;
	inpPos[1] = py ;
	inpPos[2] = pz ;

	inpTree->Fill();
	numInpTreeFilled++;
	//}
}
