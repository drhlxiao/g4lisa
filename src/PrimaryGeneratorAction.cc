/***************************************************************
 * descriptions of partcile gun
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/

#include "PrimaryGeneratorAction.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "t2sim.h"
// For Random Generator
#include <CLHEP/Random/RandFlat.h>

#include <fstream>

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : sourceType(0), fTree(NULL), fFile(NULL), nEntries(0), iEntry(0)
{
  // fParticleGunMessenger = new PrimaryGeneratorMessenger(this);
  fHistoEnergy = NULL;

  // particle source
  particleSource = "";
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  // default particle kinematic
  particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition *particle =
      particleTable->FindParticle(particleName = "gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
  fParticleGun->SetParticleEnergy(30. * keV);
  // also defined in the DetectorConstruction
  fParticleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, 15 * cm));
  // particle source
  fParticleSource = new G4GeneralParticleSource();
  // fParticleSource->SetParticleEnergy(50*keV);
  histoEnergy = NULL;
}
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleSource;
  // delete fParticleGunMessenger;
  delete fParticleGun;
  if (fHistoEnergy)
    fHistoEnergy->Close();
}

void PrimaryGeneratorAction::InitParticleSpectrumFromROOT(G4String val)
{
}
void PrimaryGeneratorAction::GetGPS(G4ThreeVector &position,
                                    G4ThreeVector &direction,
                                    G4double &energy)
{
  // noted that this dosen't work with older versions of geant4
  if (particleSource == "")
  {
    position = fParticleSource->GetParticlePosition();
    direction = fParticleSource->GetParticleMomentumDirection();
    energy = fParticleSource->GetParticleEnergy() / keV;
  }
  else
  {
    position = fParticleGun->GetParticlePosition();
    direction = fParticleGun->GetParticleMomentumDirection();
    energy = fParticleGun->GetParticleEnergy() / keV;
  }
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  G4double energy;
  G4double rnd;
  G4double phi, theta, px, py, pz, radius;
  G4double sourcePlaneRadius = 157 / 2;
  //

  // particleTable->FindParticle("gamma");
  // fParticleGun->SetParticleDefinition(particle);

  if (particleSource == "fromROOT")
  {
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    fParticleSource->GeneratePrimaryVertex(anEvent);
  }
}

G4bool PrimaryGeneratorAction::InitFile() {}
