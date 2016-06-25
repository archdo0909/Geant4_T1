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
#ifndef T1DETECTORCALSDHIT_HH
#define T1DETECTORCALSDHIT_HH

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4ParticleDefinition.hh"
#include "G4VHit.hh"

class BeamTestSiliconMonitorHit : public G4VHit {

public:

	// Constructors
	BeamTestSiliconMonitorHit();
	BeamTestSiliconMonitorHit(G4int id);

	// Destructor
	virtual ~BeamTestSiliconMonitorHit();

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	// Methods
	virtual void Draw();

	virtual void Print();

	// Incidence Information 
	inline void SetIncidenceDefinition(G4ParticleDefinition* pd) {fIPD = pd;}
	inline G4ParticleDefinition* GetIncidenceDefinition() const {return fIPD;}

	inline void SetIncidenceKineticEnergy(G4double e) {fIKEnergy = e;}
	inline G4double GetIncidenceKineticEnergy() const {return fIKEnergy;}

	inline void AddDepositedEnergy(G4double energy) {fDepositedEnergy += energy;}
	inline G4double GetDepositedEnergy() const {return fDepositedEnergy;}

	inline void SetIncidencePosition(G4ThreeVector position) {fIPosition = position;}
	inline G4ThreeVector GetIncidencePosition() const {return fIPosition;}

	inline void SetIncidenceMomentumDirection(G4ThreeVector momentum) {fIMomentumD = momentum;}
	inline G4ThreeVector GetIncidenceMomentumDirection() const {return fIMomentumD;}

private:

	// Data members
	G4int fCellID;
	G4ParticleDefinition* fIPD;
	G4double fDepositedEnergy;
	G4double fIKEnergy;
	G4ThreeVector fIPosition;
	G4ThreeVector fIMomentumD;

};

typedef G4THitsCollection<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitsCollection;

extern G4Allocator<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitAllocator;

inline void* BeamTestSiliconMonitorHit::operator new(size_t)
{
	void* aHit;
	aHit = (void*)BeamTestSiliconMonitorHitAllocator.MallocSingle();
	return aHit;
}

inline void BeamTestSiliconMonitorHit::operator delete(void* aHit)
{
	BeamTestSiliconMonitorHitAllocator.FreeSingle((BeamTestSiliconMonitorHit*) aHit);
}

#endif
