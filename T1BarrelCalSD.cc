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
/// \file persistency/P01/src/ExP01TrackerSD.cc
/// \brief Implementation of the ExP01TrackerSD class
//
//
// $Id: ExP01TrackerSD.cc 71111 2013-06-11 10:51:02Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "T1BarrelCalSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4ios.hh"

//#include "RootIO.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

T1BarrelCalSD::T1BarrelCalSD(G4String name)
	:G4VSensitiveDetector(name), fTrackerCollection(0), fHCID(0)
{
	G4String HCname = name + "_HC";
	collectionName.insert(HCname);
	G4cout << collectionName.size() << "   CalorimeterSD name:  " << name << " collection Name: " 
		<< HCname << G4endl;
	fHCID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

T1BarrelCalSD::~T1BarrelCalSD()
{ 
	//RootIO::GetInstance()->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void T1BarrelCalSD::Initialize(G4HCofThisEvent* HCE)
{
	//create hit collections
	fTrackerCollection = new T1BarrelCalHitsCollection
		(SensitiveDetectorName,collectionName[0]); 

	if (fHCID < 0) {
		G4cout << "CalorimeterSD::Initialize:  " << SensitiveDetectorName << "   " 
			<< collectionName[0] << G4endl;
		fHCID = GetCollectionID(0);
	}

	HCE->AddHitsCollection(fHCID, fTrackerCollection);

	//clear energy deposit buffer
	for (G4int i=0; i< NCHANNEL_BCAL; i++) edepbuf[i]=0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool T1BarrelCalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
#if 0
	//////////////////////////////////////////////////////////////////////////
	//Get energy deposited in this step
	// energy deposit
	G4double edep = aStep->GetTotalEnergyDeposit();

	// step length
	G4double stepLength = 0.;
	if ( aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
		stepLength = aStep->GetStepLength();
	}

	if ( edep==0. && stepLength == 0. ) return false;      

	G4TouchableHistory* touchable
		= (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

	// Get calorimeter cell id 
	G4int layerNumber = touchable->GetReplicaNumber(1);

	// Get hit accounting data for this cell
	T1BarrelCalHit* hit = (*fTrackerCollection)[layerNumber];
	if ( ! hit ) {
		G4ExceptionDescription msg;
		msg << "Cannot access hit " << layerNumber; 
		G4Exception("B4cCalorimeterSD::ProcessHits()",
			"MyCode0004", FatalException, msg);
	}         

	// Get hit for total accounting
	T1BarrelCalHit* hitTotal 
		= (*fTrackerCollection)[fTrackerCollection->entries()-1];

	// Add values
	hit->Add(edep, stepLength);
	hitTotal->Add(edep, stepLength); 
#endif
	//////////////////////////////////////////////////////////////////////////

#if 1

	G4double edep = aStep->GetTotalEnergyDeposit();

	if(edep==0.) return false;

	T1BarrelCalHit* newHit = new T1BarrelCalHit();
	newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
	newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()
		->GetReplicaNumber());
	newHit->SetEdep     (edep);
	newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
	newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
	fTrackerCollection->insert( newHit );

	//newHit->Print();
	//newHit->Draw();
#endif

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void T1BarrelCalSD::EndOfEvent(G4HCofThisEvent* HCTE)
{
	// storing the hits in ROOT file
	G4int NbHits = fTrackerCollection->entries();
	std::vector<T1BarrelCalHit*> hitsVector;

	{ 
		G4cout << "\n-------->Storing hits in the ROOT file: in this event there are " << NbHits 
			<< " hits in the tracker chambers: " << G4endl;
		for (G4int i=0;i<NbHits;i++) (*fTrackerCollection)[i]->Print();
	} 


	for (G4int i=0;i<NbHits;i++) 
		hitsVector.push_back((*fTrackerCollection)[i]);

	//RootIO::GetInstance()->Write(&hitsVector);

	//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

