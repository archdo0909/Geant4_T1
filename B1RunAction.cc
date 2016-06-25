//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "T1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "g4root.hh"

#include <iostream>
std::ofstream ofs;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction()
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);        

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //

  //// Creating histograms
  //analysisManager->CreateH1("1","Edep in absorber", 100, 0., 800*MeV);
  //analysisManager->CreateH1("2","Edep in gap", 100, 0., 100*MeV);
  //analysisManager->CreateH1("3","trackL in absorber", 100, 0., 1*m);
  //analysisManager->CreateH1("4","trackL in gap", 100, 0., 50*cm);

  //// Creating ntuple
  ////
  //analysisManager->CreateNtuple("B4", "Edep and TrackL");
  //analysisManager->CreateNtupleDColumn("Eabs");
  //analysisManager->CreateNtupleDColumn("Egap");
  //analysisManager->CreateNtupleDColumn("Labs");
  //analysisManager->CreateNtupleDColumn("Lgap");
  //analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* B1RunAction::GenerateRun()
{
  return new B1Run; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "T1";
  analysisManager->OpenFile(fileName);

  ofs.open("a.dat", std::ios::out);
  if(! ofs.good()) {
	  G4String errorMessage= "*** fail to open a file (a.out).";
	  //G4Exception(errorMessage);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  //G4int nofEvents = run->GetNumberOfEvent();
  //if (nofEvents == 0) return;
  //
  //const B1Run* b1Run = static_cast<const B1Run*>(run);

  //// Compute dose
  ////
  //G4double edep  = b1Run->GetEdep();
  //G4double edep2 = b1Run->GetEdep2();
  //G4double rms = edep2 - edep*edep/nofEvents;
  //if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  //const B1DetectorConstruction* detectorConstruction
  // = static_cast<const B1DetectorConstruction*>
  //   (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  ////G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  //G4double dose = edep/1;
  //G4double rmsDose = rms/1;

  //// Run conditions
  ////  note: There is no primary generator action object for "master"
  ////        run manager for multi-threaded mode.
  ////const B1PrimaryGeneratorAction* generatorAction
  // //= static_cast<const B1PrimaryGeneratorAction*>
  //const BeamTestPrimaryGeneratorAction* generatorAction
  //= static_cast<const BeamTestPrimaryGeneratorAction*>
  //   (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  //G4String runCondition;
  //if (generatorAction)
  //{
  //  const G4ParticleGun* particleGun = generatorAction->Gun();
  //  runCondition += particleGun->GetParticleDefinition()->GetParticleName();
  //  runCondition += " of ";
  //  G4double particleEnergy = particleGun->GetParticleEnergy();
  //  runCondition += G4BestUnit(particleEnergy,"Energy");
  //}
  //      
  //// Print
  ////  
  //if (IsMaster()) {
  //  G4cout
  //   << G4endl
  //   << "--------------------End of Global Run-----------------------";
  //}
  //else {
  //  G4cout
  //   << G4endl
  //   << "--------------------End of Local Run------------------------";
  //}
  //
  //G4cout
  //   << G4endl
  //   << " The run consists of " << nofEvents << " "<< runCondition
  //   << G4endl
  //   << " Dose in scoring volume : " 
  //   << G4BestUnit(dose,"Dose") << " +- " << G4BestUnit(rmsDose,"Dose")
  //   << G4endl
  //   << "------------------------------------------------------------"
  //   << G4endl
  //   << G4endl;

  //// print histogram statistics
  ////
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //if ( analysisManager->GetH1(1) ) {
	 // G4cout << G4endl << " ----> print histograms statistic ";
	 // if(isMaster) {
		//  G4cout << "for the entire run " << G4endl << G4endl; 
	 // }
	 // else {
		//  G4cout << "for the local thread " << G4endl << G4endl; 
	 // }

	 // G4cout << " EAbs : mean = " 
		//  << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
		//  << " rms = " 
		//  << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

	 // G4cout << " EGap : mean = " 
		//  << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
		//  << " rms = " 
		//  << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;

	 // G4cout << " LAbs : mean = " 
		//  << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
		//  << " rms = " 
		//  << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

	 // G4cout << " LGap : mean = " 
		//  << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
		//  << " rms = " 
		//  << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
  //}

  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  ofs.close();
  G4cout << ">>> #events generated= " 
	  << run-> GetNumberOfEvent() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
