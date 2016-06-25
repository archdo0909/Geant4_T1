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
// GEANT4 tag $Name:$
// 
// T. Aso        Original author
// 
#ifndef T1PRIMARYGENERATORACTION_HH
#define T1PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class BeamTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

	// Constructor
	BeamTestPrimaryGeneratorAction();    

	// Destructor
	virtual ~BeamTestPrimaryGeneratorAction();

	// Methods
	void GeneratePrimaries(G4Event*);

	static G4ParticleGun* Gun() {return particleGun;}

private:

	// Data member
	static G4ParticleGun* particleGun;

};

#endif


