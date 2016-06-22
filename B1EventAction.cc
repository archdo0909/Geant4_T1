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
// $Id: B1EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "T1BarrelCalSD.hh"
#include "G4SystemOfUnits.hh"

extern std::ofstream ofs;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction()
: G4UserEventAction(),
  fEdep(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* anEvent)
{   
	G4cout << ">>> Event ID: " << anEvent-> GetEventID() << G4endl;


	G4SDManager* SDManager= G4SDManager::GetSDMpointer();

	// ====================================================================  
	// print out hit information
	// ====================================================================  
	// for barrel calorimeter
	T1BarrelCalSD* barrelcalSD = (T1BarrelCalSD*)
		SDManager->FindSensitiveDetector("TrackerChamberSD");
	barrelcalSD-> PrintAll();


	// ====================================================================  
	// output "Hit Collection of This Event" to a data file
	// ====================================================================  
	

	// get "Hit Collection of This Event"
	G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
	if(! HCTE) return;  // no hits in this events. nothing to do!

	// for barrel calorimeter
	// [ E0 E1 E2 E3 E4 E5 E6 E7 ] (deposit energy in each module in MeV)
	static G4int idcal= -1;
	if(idcal<0)
		idcal = SDManager->GetCollectionID("TrackerChamberSD_HC"); //idcal= 0; //SDManager-> GetCollectionID("/barrelCal_HC");

	T1BarrelCalHitsCollection* hccal = (T1BarrelCalHitsCollection*)HCTE-> GetHC(idcal);
	G4double edep[NCHANNEL_BCAL];
	
	G4int idx;
	for (idx=0; idx< NCHANNEL_BCAL; idx++) {
		edep[idx]= 0.;
	}

	if(hccal) {
		G4int nhits= hccal-> entries();
		for(idx=0; idx< nhits; idx++) {
			G4int ich= (*hccal)[idx]-> GetID();
			edep[ich]= (*hccal)[idx]-> GetEdep();
			G4ThreeVector position = (*hccal)[idx]->GetPos();
			ofs << position <<edep[ich] <<ich<<" ";

		}
	}
	ofs << G4endl;

	// output to a file
	for (idx=0; idx< NCHANNEL_BCAL; idx++) {
		ofs << edep[idx]/MeV << " ";
	}
	ofs << G4endl;

  // accumulate statistics in B1Run
  B1Run* run 
    = static_cast<B1Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEdep(fEdep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
