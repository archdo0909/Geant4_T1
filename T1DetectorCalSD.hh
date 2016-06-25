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
// Tatsumi Koi
//
#ifndef T1DETECTORCALSD_HH
#define T1DETECTORCALSD_HH

#include "G4VSensitiveDetector.hh"
#include "T1DetectorCalSDHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class BeamTestSiliconMonitor : public G4VSensitiveDetector {

public:

	// Constructor
	BeamTestSiliconMonitor(const G4String& name);

	// Destructor
	virtual ~BeamTestSiliconMonitor();

	// Methods
	virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);

	virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

	virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);

private:

	// Data members
	BeamTestSiliconMonitorHitsCollection* fHitsCollection;
	G4int fHitsCollectionID;

};



#endif

