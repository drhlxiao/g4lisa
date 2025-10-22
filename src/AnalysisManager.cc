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
// number of photons to fill to the tracking tree

const G4double highVoltage =
300;                      // CdTe HV is 300 during the nominal operations
const G4double ENOISE = 0.43; // keV
							  // 0.52 is from the best fit
							  // from gussian fit, it should be 1.2/2.35=0.5

const G4double FANO_FACTOR = 0.15; // from best fit

G4double NEAR_SURFACE_L =
5.28e-3; // see Oliver's paper, in units of mm, take the mean value
G4double NEAR_SURFACE_R0 =
0.116; // see Oliver's paper, in units of mm, mean value are take
	   // best FIT L-=5.28e-3,R0=0.1 , set enoise=0.52, fanao=0.15

	   // CdTe fano factor is 0.15 according to
	   // https://www.researchgate.net/figure/Fano-factor-for-different-semiconductor-at-room-temperature_tbl5_343053397
	   // On Jun 27, the enoise and fano factor were found to be .24 and 0.74 through
	   // fitting of Ba133 exp sim spectrum, see:
	   // ~/FHNW/STIX/SolarFlareAnalysis/fitG4CalibrationSpectrum

const G4double threshold = 4;
/* At 30 keV, the energy resolution is 2%  * 30 keV= 0.6 keV
 * using fano factor, one could know the intrinsic resolution of CdTe is
 * rho=(pairs * 0.15)/pairs =0.0047 absolute resolution is rho*30 = 0.14 keV ADC
 * resolution is 0.5 ADC channel,this is equivalent  to 0.5 /2.3 =0.2 keV for
 * calibration spectrum electronics noise = sqrt(0.6*0.6 -0.14*0.14 - 0.2*0.2)
 */

const G4double CdTe_SURFACE_X =
12.7741; // surface x-coordinates of CdTe detectors, using geant4 tracks to
		 // find the position, 2023-06-26, not it can be 13.774
const G4double PY_ORIGIN = 103.1;
const G4double PZ_ORIGIN = 127.5;

const G4double PAIR_CREATION_ENERGY = 4.46e-3; // in units of keV

double energyRanges[] = {0,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13,
	14, 15, 16, 18, 20, 22, 25, 28,  32,  36,  40,
	45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 250};

void normalizedEnergySpectrum(TH1F *h) {
	// energy bin widths are different, we need to  divide the counts by energy
	// bin width
	int nbins = h->GetXaxis()->GetNbins();
	if (nbins != 32) {
		G4cout << "Can not normalize histogram" << G4endl;
		return;
	}
	for (int i = 0; i < 32; i++) {
		double binW = energyRanges[i + 1] - energyRanges[i];
		if (binW > 0) {
			h->SetBinContent(i + 1, h->GetBinContent(i + 1) / binW);
		}
	}
}

int getScienceBin(Double_t energy) {
	if (energy < 4) {
		return 0;
	} else if (energy >= 150) {
		return 31;
	} else {
		for (int i = 0; i < 31; i++) {
			if (energy >= energyRanges[i] && energy < energyRanges[i + 1])
				return i;
		}
	}
	return 31;
}

const double histMaxEnergy = 150;
int histNbins = (int)(histMaxEnergy / 0.1);

AnalysisManager *AnalysisManager::fManager = 0;
AnalysisManager *AnalysisManager::GetInstance() {
	if (!fManager) {
		fManager = new AnalysisManager();
	}
	return fManager;
}
AnalysisManager::AnalysisManager() {

	//isDummyDetector=false;
	numKilled = 0;
	isNewEvent = true;
	numEventIn = 0;
	numEventOut = 0;
	killTracksEnteringGrids = false;
	killTracksEnteringDetectors = false;
	numInpTreeFilled = 0;
	numSourceTreeFilled= 0;

	numPhysTreeFilled = 0;
}
void AnalysisManager::CopyMacrosToROOT(TFile *f, TString &macfilename) {
	G4cout << "Copying macros from file " << macfilename << " to root file"
		<< G4endl;
	if (macfilename == "" || !f)
		return;
	std::ifstream infile(macfilename.Data());
	if (!infile.good()) {
		G4cout << "can not open the macro file, existing..." << G4endl;
		return;
	}
	G4String macros = Form("\nCommand:%s\n", commandLine.data());
	macros += Form("\nMacro filename: %s\n ", macroFilename.Data());
	std::string line;
	G4cout << "Macros read from the macro file:" << G4endl;
	while (std::getline(infile, line)) {
		macros += line + "\n";
		G4cout << line << G4endl;
	}
	macros += "---------------  Realistic simulation parameters---------------";
	macros += Form("\nNear surface R0: %f\n ", NEAR_SURFACE_R0);
	macros += Form("\nNear surface L: %f\n ", NEAR_SURFACE_L);
	macros += Form("\nFanno factor : %f\n ", FANO_FACTOR);
	macros += Form("\nENOIS (eV): %f\n ", ENOISE);
	TNamed cmd;
	cmd.SetTitle(macros);
	f->cd();
	cmd.Write("metadata");
	infile.close();
}

void AnalysisManager::InitRun(const G4Run *run) {
	rootFile = new TFile(outputFilename.Data(), "recreate");

	evtTree = new TTree("events", "events");
	evtTree->Branch("edep", edepSum, Form("edep[%d]/D", NUM_CHANNELS));
	evtTree->Branch("sci", sci, Form("sci[%d]/D", NUM_CHANNELS));
	//	evtTree->Branch("colEff", colEff, Form("colEff[%d]/D", NUM_CHANNELS));
	evtTree->Branch("collected", collectedEnergySum,
			Form("collected[%d]/D", NUM_CHANNELS));
	evtTree->Branch("charge", edepWithoutNoise,
			Form("charge[%d]/D", NUM_CHANNELS));
	evtTree->Branch("charge2", collectedEdepSumRealistic,
			Form("charge2[%d]/D", NUM_CHANNELS));
	evtTree->Branch("eventID", &eventID, "eventID/I");
	evtTree->Branch("E0", &gunEnergy, "E0/D");
	evtTree->Branch("gunPos", gunPosition, Form("gunPos[%d]/D", 3));
	evtTree->Branch("gunVec", gunDirection, Form("gunVec[%d]/D", 3));
	evtTree->Branch("numTracks", &itrack, Form("numTracks/I"));
	evtTree->Branch("nHits", nHits, Form("nHits[32]/I"));
	evtTree->Branch("hitx", hitx, Form("hitx[%d]/D", MAX_TRACKS));
	evtTree->Branch("hity", hity, Form("hity[%d]/D", MAX_TRACKS));
	evtTree->Branch("hitz", hitz, Form("hitz[%d]/D", MAX_TRACKS));
	evtTree->Branch("parent", parent, Form("parent[%d]/I", MAX_TRACKS));
	evtTree->Branch("pdg", pdg, Form("pdg[%d]/I", MAX_TRACKS));
	evtTree->Branch("pixel", hitPixelID, Form("pixel[%d]/I", MAX_TRACKS));
	evtTree->Branch("energy", energy, Form("energy[%d]/D", MAX_TRACKS));
	evtTree->Branch("time", time, Form("time[%d]/D", MAX_TRACKS));
	evtTree->Branch("totalNumSteps", &totalNumSteps, "totalNumSteps/I");

	if (DEBUG) {
		evtTree->Branch("R0", &NEAR_SURFACE_R0, "R0/D");
		evtTree->Branch("L", &NEAR_SURFACE_L, "L/D");
	}
	physTree = new TTree("phys", "phys");
	physTree->Branch("type", &processType, "type/I");
	physTree->Branch("subType", &processSubtype, "subType/I");
	physTree->Branch("E0", &gunEnergy, "gunEnergy/D");
	physTree->Branch("pdg", &inpPDG, "pdg/I");
	physTree->Branch("parent", &parentID, "parentID/I");

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

	c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	/*	for (int i = 0; i < NUM_CHANNELS; i++) {
		hd[i] = new TH1F(Form("hd%d", i),
		Form("Spectrum of pixel %d energy depositions; Energy "
		"deposition (keV); Counts",
		i),
		200, 0, 500);
		}
		*/
	TString channelName = "";
	for (int i = 0; i < 34; i++) {
		if (i == 8) {
			channelName = "CFL";
		} else if (i == 9) {
			channelName = "BKG";
		} else if (i == 32) {
			channelName = "Detector summed";
		} else if (i == 33) {
			channelName = "Big pixels except CFL and BKG";
		} else {
			channelName = Form("D%d", i);
		}

	   hpat[i] = new TH1F(Form("hpat_%d", i),
	   "Detector count pattern; Pixel #; counts;", 12, 0, 12);

		hRealSci[i] = new TH1F(Form("hRealSci%d", i),
				Form("Recorded energy spectrum - Rebinned to SCI "
					"channels (%s); Energy (keV)",
					channelName.Data()),
				32, energyRanges);
		hEdepSci[i] = new TH1F(Form("hEdepSci%d", i),
				Form("Deposited energy spectrum - Rebinned to SCI "
					"channels (%s); Energy (keV)",
					channelName.Data()),
				32, energyRanges);

		hReal[i] = new TH1F(Form("hReal%d", i),
				Form("Recorded energy spectrum  (%s); Energy (keV)",
					channelName.Data()),
				histNbins, 0, histMaxEnergy);
		hEdep[i] = new TH1F(Form("hEdep%d", i),
				Form("Deposited energy spectrum (%s); Energy (keV)",
					channelName.Data()),
				histNbins, 0, histMaxEnergy);

		hRealSciSingleHit[i] =
			new TH1F(Form("hRealSciSingleHit%d", i),
					Form("Recorded energy spectrum - Rebinned to SCI channels "
						"(%s); Energy (keV)",
						channelName.Data()),
					32, energyRanges);
		hEdepSciSingleHit[i] =
			new TH1F(Form("hEdepSciSingleHit%d", i),
					Form("Deposited energy spectrum - Rebinned to SCI channels "
						"(%s); Energy (keV)",
						channelName.Data()),
					32, energyRanges);
		hRealSingleHit[i] =
			new TH1F(Form("hRealSingleHit%d", i),
					Form("Recorded energy spectrum  (%s); Energy (keV)",
						channelName.Data()),
					histNbins, 0, histMaxEnergy);
		hEdepSingleHit[i] =
			new TH1F(Form("hEdepSingleHit%d", i),
					Form("Deposited energy spectrum (%s); Energy (keV)",
						channelName.Data()),
					histNbins, 0, histMaxEnergy);
	}
	h2xy = new TH2F("h2xy", "Locations of hits; X (mm); Y(mm)", 1800, -90, 90,
			1800, -90, 90);
	hz = new TH1F("h1depth", "Energy deposition depth; Depth (mm); Counts;", 100,
			0, 1);
	hEdepSum = new TH1F("hEdepSum",
			"Detector summed energy spectrum; Energy (keV); Counts;",
			200, 0, 100);
	hdc = new TH1F("h1DetTotCnts", "Detector total counts; Detector ID; counts;",
			32, 0, 32);
	hpc = new TH1F("h1PixelTotCnts", "Pixel total counts; Pixel ID; counts;", 384,
			0, 384);
	hcol = new TH1F("h1ChargeColEff",
			Form("Distribution of Charge collection efficiency (Fanno: "
				"%f, ENOISE: %f) ; Efficiency; Counts ;",
				FANO_FACTOR, ENOISE),
			200, 0, 1);
	hNS =
		new TH1F("hNearSurfaceFactor",
				Form("CF of surface effect (L:%f; R0:%f); Efficiency; Counts ;",
					NEAR_SURFACE_L, NEAR_SURFACE_R0),
				200, 0, 1);

	hEdepSum->SetCanExtend(TH1::kXaxis);

	// for ROOT version >6.0
}

//////////////////////////////////////////////////////////////////////////

void AnalysisManager::ProcessRun(const G4Run *run) {
	CloseROOT();
	G4cout << "Events entered detectors:" << numEventIn << G4endl;
	G4cout << "Events escaped from detectors:" << numEventOut << G4endl;
	G4cout << "Events captured:" << numEventIn - numEventOut << G4endl;
}

/// EventAction
void AnalysisManager::InitEvent(const G4Event *event) {

	if (DEBUG) {
		NEAR_SURFACE_R0 = 0.1 + 0.8 * G4UniformRand();
		NEAR_SURFACE_L = (5 + 3.5 * G4UniformRand()) * 1e-3;
	}
	for (G4int i = 0; i < NUM_CHANNELS; i++) {
		edepSum[i] = 0.0;
		collectedEnergySum[i] = 0;
		edepWithoutNoise[i] = 0;
		collectedEdepSumRealistic[i] = 0;
		sci[i] = -1;
	}
	for (int i = 0; i < 32; i++)
		nHits[i] = 0;

	for (int i = 0; i < MAX_TRACKS; i++) {
		pdg[i] = 0;
		energy[i] = 0;
		hitPixelID[i] = 0;
		parent[i] = 0;
		time[i] = 0;
		hitx[i] = 0;
		hity[i] = 0;
		hitz[i] = 0;
	}
	itrack = 0;
	totalNumSteps = 0;
	isNewEvent = true;
}
void AnalysisManager::ProcessEvent(const G4Event *event) {
	eventID = event->GetEventID();
	G4bool effectiveEvent= false;
	for (int i = 0; i < NUM_CHANNELS; i++) {
		if (edepSum[i] > 0) {
			hEdepSum->Fill(edepSum[i]);
			// hd[i]->Fill(edepSum[i]);
			effectiveEvent= true;
			G4double sigma = GetEnergyResolution(collectedEnergySum[i]);
			// std. Deviation of charge, in units of keV
			// randomized the energy

			edepWithoutNoise[i] = gRandom->Gaus(collectedEnergySum[i], sigma);
			// charge
			//	PAIR_CREATION_ENERGY /1000;
			//  convert back to keV, we randomize it to
			//  smear the energy resolution
			//
			collectedEdepSumRealistic[i] = gRandom->Gaus(edepWithoutNoise[i], ENOISE);
			// we asssue the electronics noise

			detectorID = i / 12;
			pixelID = i % 12;
			hpc->Fill(i);
			hdc->Fill(detectorID);

			if (collectedEdepSumRealistic[i] > threshold) {
				nHits[detectorID]++;
			}
			sci[i] = getScienceBin(collectedEdepSumRealistic[i]);

			hEdep[detectorID]->Fill(edepSum[i]);
			hReal[detectorID]->Fill(collectedEdepSumRealistic[i]);
			// not binned to stix
			hEdep[32]->Fill(edepSum[i]);
			hReal[32]->Fill(collectedEdepSumRealistic[i]);
			// histograms to store spectra of events of all detectors

			hEdepSci[detectorID]->Fill(edepSum[i]);
			hRealSci[detectorID]->Fill(collectedEdepSumRealistic[i]);
			//
			hEdepSci[32]->Fill(edepSum[i]);
			hRealSci[32]->Fill(collectedEdepSumRealistic[i]);

			if (detectorID != 8 && detectorID != 9 && pixelID < 8) {
				// big pixel except CFL and BKG
				hEdepSci[33]->Fill(edepSum[i]);
				hRealSci[33]->Fill(collectedEdepSumRealistic[i]);
			}
			// summed spectrum
		}
	}
	// iterate over detectors
	for (int i = 0; i < 32; i++) {
		if (nHits[i] != 1)
			continue;
		// only accepts single hit events
		for (int j = 0; j < 12; j++) {
			int ch = i * 12 + j;
			if (edepSum[ch] < 0)
				continue;

			// hEdep[detectorID]->Fill(edepSum[i]);
			// hReal[detectorID]->Fill(collectedEdepSumRealistic[i]);

			hEdepSingleHit[i]->Fill(edepSum[ch]);
			hRealSingleHit[i]->Fill(collectedEdepSumRealistic[ch]);
			hEdepSciSingleHit[i]->Fill(edepSum[ch]);
			hRealSciSingleHit[i]->Fill(collectedEdepSumRealistic[ch]);
			// stix energy bins

			hEdepSingleHit[32]->Fill(edepSum[ch]);
			hRealSingleHit[32]->Fill(collectedEdepSumRealistic[ch]);
			hEdepSciSingleHit[32]->Fill(edepSum[ch]);
			hRealSciSingleHit[32]->Fill(collectedEdepSumRealistic[ch]);
			// sum spectrum

			if (j < 8 && i != 8 && i != 9) {

				hEdepSingleHit[33]->Fill(edepSum[ch]);
				hRealSingleHit[33]->Fill(collectedEdepSumRealistic[ch]);
				hEdepSciSingleHit[33]->Fill(edepSum[ch]);
				hRealSciSingleHit[33]->Fill(collectedEdepSumRealistic[ch]);
			}
		}
	}

	///	hdc->Fill(detectorID);
	// toFill=true;

	if (numSourceTreeFilled < 100000)  {
		primTree->Fill();
		numSourceTreeFilled++;
	}
	//if (effectiveEvent)  evtTree->Fill();
	isNewEvent = true;
}
// Stepping Action
void AnalysisManager::UpdateParticleGunInfo() {
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

G4double AnalysisManager::GetNearSurfaceFactor(G4ThreeVector &pos) {
	G4double xx = pos.x() / mm;
	G4double z = abs(xx - CdTe_SURFACE_X); // depth from cathode
	G4double factor = 1 - NEAR_SURFACE_R0 * exp(-z / NEAR_SURFACE_L);
	hNS->Fill(factor);
	return factor;
}
G4double AnalysisManager::ComputeCollectionEfficiency(G4ThreeVector &pos) {
	// Hecht equation, see "Recent Progress in CdTe and CdZnTe Detectors"
	// Tadayuki Takahashi and Shin Watanabe
	// Oliver's paper
	// Energy must be in units of eV
	// Spectral signature of near-surface damage in CdTe X-ray detectors
	//
	G4double xx = pos.x() / mm;
	// G4cout<<xx<<G4endl;
	G4double z = abs(xx - CdTe_SURFACE_X); // depth from cathode
										   //
	hz->Fill(z);

	G4double freePathElectron = 1100 * 100 * 3e-6 * highVoltage;
	G4double freePathHoles = 100 * 100 * 2e-6 * highVoltage;

	G4double d = 1;

	G4double eff =
		(1 - exp((z - d) / freePathElectron)) * (freePathElectron / d) +
		(freePathHoles / d) * (1.0 - exp(-z / freePathHoles));

	hcol->Fill(eff);
	return eff;
}

G4double AnalysisManager::GetEnergyResolution(G4double edep) {
	//	Recent Progress in CdTe and CdZnTe Detectors
	// Tadayuki Takahashi and Shin Watanabe
	// energy in units of keV
	//
	/*	G4double rho = 0;
		if (edep > 0) {
		G4double numPairs = edep / PAIR_CREATION_ENERGY * 1000;
	// keV to ev
	//  see oliver's paper
	G4double numSigma2 = FANO_FACTOR * numPairs;
	// sigma^2 = fano * n  according to the text-book
	rho = sqrt(numSigma2) / numPairs;
	}
	return rho;
	The above is equv. to  delta_E=2.35 sqrt(FWE), see
https://www.sciencedirect.com/topics/neuroscience/fano-factor
*/
	//return 2.35 * sqrt(FANO_FACTOR * PAIR_CREATION_ENERGY * edep);
	return sqrt(FANO_FACTOR * PAIR_CREATION_ENERGY * edep);
}

////////////////////////////////////////////////////////////////////

inline void AnalysisManager::AddEnergy(G4int detId, G4double edep) {
	if (detId > NUM_CHANNELS || detId < 0) {
		G4cout << "invalid index" << G4endl;
		return;
	}
	edepSum[detId] += edep;
}
inline void AnalysisManager::AddCollectedEnergy(G4int detId, G4double dep) {
	if (detId > NUM_CHANNELS || detId < 0) {
		G4cout << "invalid index" << G4endl;
		return;
	}
	collectedEnergySum[detId] += dep;
}

AnalysisManager::~AnalysisManager() {
	if (fManager)
		delete fManager;
	// if (rootFile) delete rootFile;
	fManager = 0;
}

void AnalysisManager::CloseROOT() {
	rootFile->cd();
	evtTree->Write();
	inpTree->Write();
	primTree->Write();
	physTree->Write();
	TDirectory *cdhist = rootFile->mkdir("hist");
	cdhist->cd();
	c1->cd();
	c1->Divide(2, 5);

	for (int i = 0; i < 34; i++) {

		hReal[i]->Write();
		hEdep[i]->Write();

		normalizedEnergySpectrum(hRealSci[i]);
		normalizedEnergySpectrum(hEdepSci[i]);
		normalizedEnergySpectrum(hEdepSciSingleHit[i]);
		normalizedEnergySpectrum(hRealSciSingleHit[i]);
		// bin

		hRealSci[i]->Write();
		hEdepSci[i]->Write();
		hRealSciSingleHit[i]->Write();
		hEdepSciSingleHit[i]->Write();

		hRealSingleHit[i]->Write();
		hEdepSingleHit[i]->Write();
		hpat[i]->Write();

	}
	hEdepSum->Write();
	h2xy->Write();
	hdc->Write();
	hpc->Write();
	c1->Write();
	hz->Write();
	hcol->Write();
	hNS->Write();
	G4cout << ">> Number of event recorded:" << evtTree->GetEntries() << G4endl;
	G4cout << ">> Number of incident particles :" << inpTree->GetEntries()
		<< G4endl;
	G4cout << ">> Number of track killed:" << numKilled << G4endl;
	CopyMacrosToROOT(rootFile, macroFilename);

	rootFile->Close();
}

void AnalysisManager::ProcessStep(const G4Step *aStep) {

	UpdateParticleGunInfo();

	G4double px, py, pz;
	G4double edep;
	G4String volName;
	const G4Track *track = aStep->GetTrack();

	if (track->GetVolume())
		volName = track->GetVolume()->GetName();


	if (volName=="detector") { // don't store too many tracks
		FillDetectorIncidentParticle(aStep);
		aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		numKilled++;
	}
		//dummy detector, material as black hole
		//used to study grid effect

	






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
	h2xy->Fill(py - PY_ORIGIN, pz - PZ_ORIGIN);

	//if (numInpTreeFilled < MAX_NUM_TREE_TO_FILL) {
		G4ThreeVector inpV = track->GetMomentumDirection();
		inpVec[0] = inpV.x();
		inpVec[1] = inpV.y();
		inpVec[2] = inpV.z();
		inpTheta = asin(sqrt(inpVec[1] * inpVec[1] + inpVec[2] * inpVec[2])) * 180 /
			3.1415926;

		inpPos[0] = px;
		inpPos[1] = py - PY_ORIGIN;
		inpPos[2] = pz - PZ_ORIGIN;

		inpTree->Fill(); numInpTreeFilled++;
	//}
}

	// physics processes
	//
	//
	/*

ProcessType:

0 fNotDefined
1 fTransportation
2 fElectromagnetic
3 fOptical
4 fHadronic
5 fPhotolepton_hadron
6 fDecay
7 fGeneral
8 fParameterisation
9 fUserDefined
10 fParallel

#fComptonScattering
fCoulombScattering = 1,
fIonisation = 2,
fBremsstrahlung = 3,
fPairProdByCharged = 4,
fAnnihilation = 5,
fAnnihilationToMuMu = 6,
fAnnihilationToHadrons = 7,
fNuclearStopping = 8,
fElectronGeneralProcess = 9,

fMultipleScattering = 10,

fRayleigh = 11,
fPhotoElectricEffect = 12,
fComptonScattering = 13,
fGammaConversion = 14,
fGammaConversionToMuMu = 15,
fGammaGeneralProcess = 16,
fPositronGeneralProcess = 17,
fCerenkov = 21,
fScintillation = 22,
fSynchrotronRadiation = 23,
fTransitionRadiation = 24,
fSurfaceReflection = 25

*/
