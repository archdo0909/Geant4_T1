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
// Jane Tinslay - adapted from A01 example
// Tatsumi Koi - add Silicon Monitor Output 
//
#include "T1EventAction.hh"
#include "T1DetectorCalSDHit.hh"
////////////////////////////////////////////////////////////////////////
// HandsOn4: define Hit of Silicon Monitor
//#include "BeamTestSiliconMonitorHit.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

extern std::ofstream ofs;

BeamTestEventAction::BeamTestEventAction() {}

BeamTestEventAction::~BeamTestEventAction() {}

void BeamTestEventAction::BeginOfEventAction(const G4Event*) 
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  fHitsCollectionID = SDman->GetCollectionID("CalorimeterCollection");
  ////////////////////////////////////////////////////////////////////////
  // HandsOn4: Getting code for HitsCollection of Silicon Monitor
  //fHitsCollectionID_monitor = SDman->GetCollectionID("MonitorCollection");
}


void BeamTestEventAction::EndOfEventAction(const G4Event* event)
{

  G4HCofThisEvent *hitsCollectionOfThisEvent = event->GetHCofThisEvent();

  if (fHitsCollectionID >= 0)
  {

     BeamTestEmCalorimeterHitsCollection* hitsCollection = 
       dynamic_cast<BeamTestEmCalorimeterHitsCollection*>(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID));
  
     G4double totalEnergy = 0.;
  
     if ( 0 != hitsCollection ) {
        G4int i(0);
    
        for ( i=0 ; i<1200; i++ ) {    ////the number of i is crystal nums
           BeamTestEmCalorimeterHit* aHit = (*hitsCollection)[i];   
           totalEnergy += aHit->GetDepositedEnergy();

		if(0 != aHit->GetDepositedEnergy()){
		//ofs<< i <<", "<< aHit->GetDepositedEnergy() <<", "<<aHit->GetPosition()<< G4endl;
	
		//UP
			ofs<< i << ", " <<aHit->GetDepositedEnergy() <<", "
				<<((i-(1200/2-1))*10.0 - 5.0) << "," <<6000<<","<< 0 <<G4endl;  
		//DOWN
		//	ofs<< i << ", " <<aHit->GetDepositedEnergy() <<", "
		//		<<-((i-(1200/2-1))*10.0 - 5.0) << "," <<-6000<<","<< 0 <<G4endl;  
		////Left
		//	ofs<< i << ", " <<aHit->GetDepositedEnergy() <<", "
		//		<<-6000<< "," <<-((i-(1200/2-1))*10.0 - 5.0)<<","<< 0 <<G4endl;  
		////Right
		//	ofs<< i << ", " <<aHit->GetDepositedEnergy() <<", "
		//		<<6000<< "," <<((i-(1200/2-1))*10.0 - 5.0)<<","<< 0 <<G4endl;
          // aHit->Print();
		

		G4cout<<"Hit "<<i<<G4endl;
		}

        }
     }
     //G4cout<<"Energy deposited in calorimeter "<<totalEnergy/MeV<<" MeV"<<G4endl;
	 //ofs <<"Energy deposited in calorimeter "<<totalEnergy/MeV<<" MeV"<<G4endl;
  }

  ////////////////////////////////////////////////////////////////////////
  // HandsOn4: Output code of Silicon Monitor Hits
//  if ( fHitsCollectionID_monitor >= 0 )
//  {
//     BeamTestEmCalorimeterHit* hitsCollection_monitor = 
//        dynamic_cast<BeamTestEmCalorimeterHit*>(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID_monitor));
//      G4int numberOfHits = hitsCollection_monitor->GetSize();
//      for ( int i = 0 ; i < numberOfHits ; i++ )
//      {  
//         BeamTestEmCalorimeterHit* aHit = (*hitsCollection_monitor)[i];   
//         G4cout << "Information of " << i+1 << " th Silicon Monitor Hit of this event." << G4endl;
//		 ofs <<"Information of " << i+1 << " th Silicon Monitor Hit of this event."  <<G4endl; 
///*
//         G4cout << "Incidence Particle Name and Kinetic Energy " << aHit->GetIncidenceDefinition()->GetParticleName() 
//                << " " << aHit->GetIncidenceKineticEnergy()/MeV << " MeV" << G4endl; 
//         G4cout << "Insidence position in silicon monitor "<<aHit->GetIncidencePosition()/mm << " in mm" << G4endl; 
//         G4cout << "Incidence Direction " << aHit->GetIncidenceMomentumDirection() << G4endl; 
//*/
//         aHit->Print();
//      }
//   }
}

