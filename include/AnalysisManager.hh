//
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include <vector>

#include "globals.hh"

//#include "G4Step.hh"
//#include "G4Event.hh"
//#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "TString.h"

#define NUM_CHANNELS 12 
class G4Run;
class G4Event;
class G4Step;

class TCanvas;
class TH1F;
class TH2F;
class TFile;
class TTree;
const int MAX_TRACKS = 30;

class AnalysisManager {
public:
  static AnalysisManager *GetInstance();
  static void Dispose();

  void CloseROOT();
  void SetEventID(G4int eventid) { eventID = eventid; }
  void SetOutputFileName(TString filen) { outputFilename = filen; }

  void SetCommandLine(G4String s) { commandLine = s; }
  void InitEvent(const G4Event *event);
  void UpdateParticleGunInfo();
  void ProcessEvent(const G4Event *event);
  void ProcessStep(const G4Step *aStep);
  void InitRun(const G4Run *);
  void ProcessRun(const G4Run *);

  //void RegisterDummyDetector(){isDummyDetector=true;}
  void AddEnergy(G4int detId, G4double edep);
  void CopyMacrosToROOT(TFile *f, TString &);
  void SetMacroFileName(G4String &name) { macroFilename = name; }

void FillDetectorIncidentParticle(const G4Step *aStep);


  AnalysisManager();
  ~AnalysisManager();

private:
  static AnalysisManager *fManager;

  TString outputFilename;
  TString macroFilename;
  G4int numInpTreeFilled;
  G4int numSourceTreeFilled;
  G4int numPhysTreeFilled;
  G4int boundary;
  //G4bool isDummyDetector;

  G4bool isNewEvent;

  TFile *rootFile;
  TTree *evtTree;
  TTree *inpTree;
  TTree *physTree;
  long long totalHits;
  G4String commandLine;

  G4int itrack;
  G4int totalNumSteps;
  G4int pixelID;
  G4int detectorID;
  G4int pixelGlobalID;
  G4int inpPDG;
  G4int parentID;

  G4double inpPos[3];
  G4double inpVec[3];
  G4double inpEnergy, inpTheta;

  G4double gunPosition[3];
  G4double gunDirection[3];
  G4double gunEnergy;
  TCanvas *c1;
  TH1F *hEdepSum;
  TH1F *hz;
  TH1F *hx;
  TH1F *hNS;
  TH1F *hy;
  TH1F *hcol;
  TH1F *hpc;
  TH1F *hdc;

  TH2F *h2xy;
  TTree *primTree;
  G4int numKilled;

  G4int nHits[32];
  G4int eventID;
  G4bool killTracksEnteringGrids, killTracksEnteringDetectors;
  G4int processType, processSubtype;


  G4int numEventIn;
  G4int numEventOut;
};

#endif
