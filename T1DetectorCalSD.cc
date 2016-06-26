//
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
//
// $Id:$
//
// Tatsumi Koi for HandsOn4 of MaGill Univ. Tutorial
//
#include "T1DetectorCalSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"    
#include "G4TouchableHistory.hh"
#include "G4Track.hh"

BeamTestEmCalorimeter::BeamTestEmCalorimeter(const G4String& name)
	:G4VSensitiveDetector(name)
{
	collectionName.insert( "CalorimeterCollection" );
	fHitsCollectionID = -1;
}

BeamTestEmCalorimeter::~BeamTestEmCalorimeter() {}

void BeamTestEmCalorimeter::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
	// HandsOn4: Creating hit collection
	// Create a new collection
	fHitsCollection =
		new BeamTestEmCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);

	if(fHitsCollectionID < 0)
		fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

	// Add collection to the event
	hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);

	// Initialise hits
	G4int i(0);

	for (i=0; i<128; i++) {
		BeamTestEmCalorimeterHit* aHit = new BeamTestEmCalorimeterHit(i);
		fHitsCollection->insert(aHit);
	}
}

G4bool BeamTestEmCalorimeter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
#if 0
	//// Limited by Geom Boundary?
	//if ( aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary )
	//{
	//	// Create Hit 
	//	BeamTestSiliconMonitorHit* aHit = new BeamTestSiliconMonitorHit();
	//	fHitsCollection->insert( aHit );

	//	// Get Transportaion Matrix
	//	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	//	G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
	//	aTrans.Invert();

	//	G4ThreeVector position = aTrans.NetRotation() * ( aStep->GetPreStepPoint()->GetPosition() - aTrans.NetTranslation() ); 
	//	G4ThreeVector momentumD = aTrans.NetRotation() * aStep->GetPreStepPoint()->GetMomentumDirection();

	//	G4ParticleDefinition* pd = aStep->GetTrack()->GetDefinition(); 
	//	G4double ke = aStep->GetPreStepPoint()->GetKineticEnergy();

	//	aHit->SetIncidenceDefinition(pd);
	//	aHit->SetIncidenceKineticEnergy(ke);
	//	aHit->SetIncidencePosition(position);
	//	aHit->SetIncidenceMomentumDirection(momentumD);
	//	
	//}
	//else 
	//{
	//	return true;
	//}

	//return true;
#endif
	// HandsOn4: Accumulating hit data
	// Get energy deposited in this step
	G4double depositedEnergy = aStep->GetTotalEnergyDeposit();
	if (0 == depositedEnergy) return true;


	// Get volume and copy number
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());

	G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
	G4int copyNo = thePhysical->GetCopyNo();

	// Get corresponding hit
	BeamTestEmCalorimeterHit* aHit = (*fHitsCollection)[copyNo];

	// Check to see if this is the first time the hit has been updated
	if (!(aHit->GetLogicalVolume())) {

		// Set volume information
		aHit->SetLogicalVolume(thePhysical->GetLogicalVolume());

		G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
		aTrans.Invert();

		aHit->SetRotation(aTrans.NetRotation());
		aHit->SetPosition(aTrans.NetTranslation());
	}

	// Accumulate energy deposition  
	aHit->AddDepositedEnergy(depositedEnergy);
	return true;

}

void BeamTestEmCalorimeter::EndOfEvent(G4HCofThisEvent*) {}
