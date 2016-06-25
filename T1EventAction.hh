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
//
#ifndef T1EventAction_HH
#define T1EventAction_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

class BeamTestEventAction : public G4UserEventAction {

public:

	// Constructor
	BeamTestEventAction();

	// Destructor
	~BeamTestEventAction();

	// Metohds
	void BeginOfEventAction(const G4Event* anEvent);
	void EndOfEventAction(const G4Event* anEvent);

	void AddEdep(G4double edep) {fEdep += edep;}

private:

	// Data member
	G4int fHitsCollectionID;
	G4int fHitsCollectionID_monitor;
	G4double  fEdep;

};

#endif


