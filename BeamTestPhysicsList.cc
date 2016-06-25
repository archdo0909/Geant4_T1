// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// $Id:$
// GEANT4 tag $Name:$
//
// T. Aso        Original author
// Jane Tinslay  Minor modifications for Jefferson lab tutorial   
//
#include "BeamTestPhysicsList.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ComptonScattering.hh"
#include "G4Decay.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"
#include "G4GammaConversion.hh"
#include "G4hIonisation.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuPairProduction.hh"
#include "G4ParticleTypes.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ProcessManager.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"

BeamTestPhysicsList::BeamTestPhysicsList()
  :G4VUserPhysicsList()
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
}

BeamTestPhysicsList::~BeamTestPhysicsList() {}

void BeamTestPhysicsList::ConstructParticle()
{
  // Construct all the particles we may possibly need

  //  Baryons
  G4BaryonConstructor baryons;
  baryons.ConstructParticle();

  // Bosons (including geantinos)
  G4BosonConstructor bosons;
  bosons.ConstructParticle();

  // Ions
  G4IonConstructor ions;
  ions.ConstructParticle();

  // Leptons
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();

  // Construct all mesons
  G4MesonConstructor mesons;
  mesons.ConstructParticle();

  // Resonances and quarks
  G4ShortLivedConstructor shortLiveds;
  shortLiveds.ConstructParticle();
}

void BeamTestPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

void BeamTestPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // Tip:
    // G4ProcessManager::AddProcess(G4VProcess*, 
    //                              G4int ordAtRestDoIt, 
    //                              G4int ordAlongStepDoIt, 
    //                              G4int ordPostStepDoIt);
    // Integers in AddProcess() indicate the order in which the process 
    // is applied. If value is negative, the process is not added to 
    // the corresponding process vector.

    if (particleName == "gamma") {
      // Gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
    } else if (particleName == "e-") {
      // Electron
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);

    } else if (particleName == "e+") {
      // Positron
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      // Muon
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4, 4);

    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) &&
               (particle->GetParticleName() != "chargedgeantino")) {
      // All others charged particles except geantino
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2, 2);
    }
  }
}

void BeamTestPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();

  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      
      // Set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }    
  }
}

void BeamTestPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "BeamTestPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  // Set cut values for gamma first since some e+/e- processes
  // need gamma cut values
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  
  if (verboseLevel>0) DumpCutValuesTable();
}
