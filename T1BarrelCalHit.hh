// ====================================================================
//   T1BarrelCalHit.hh
//
// ====================================================================
#ifndef T1_BARREL_CAL_HIT_H
#define T1_BARREL_CAL_HIT_H
 
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

enum {NCHANNEL_BCAL=32};

class T1BarrelCalHit : public G4VHit {
public:

	T1BarrelCalHit();
	~T1BarrelCalHit();
	T1BarrelCalHit(const T1BarrelCalHit&);

	//operators
	const T1BarrelCalHit& operator=(const T1BarrelCalHit&);
	G4int operator==(const T1BarrelCalHit&) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	//methods from base class
	virtual void Draw();
	virtual void Print();

	//methods to handle data from B4c
	void Add(G4double de, G4double dl);

	//get methods 
	G4double GetEdep() const;
	G4double GetTrackLength() const;

public:

	void SetTrackID  (G4int track)      { fTrackID = track; };
	void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };  
	void SetEdep     (G4double de)      { fEdep = de; };
	void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

	G4int GetID()    { return fTrackID; };
	G4int GetChamberNb()  { return fChamberNb; };
	G4double GetEdep()    { return fEdep; };      
	G4ThreeVector GetPos(){ return fPos; };

private:

	G4int         fTrackID;
	G4int         fChamberNb;
	G4double      fEdep;
	G4double      fTrackLength;
	G4ThreeVector fPos;
};

typedef G4THitsCollection<T1BarrelCalHit> T1BarrelCalHitsCollection;

extern G4Allocator<T1BarrelCalHit> T1TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* T1BarrelCalHit::operator new(size_t)
{
	void *aHit;
	aHit = (void *) T1TrackerHitAllocator.MallocSingle();
	return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void T1BarrelCalHit::operator delete(void *aHit)
{
	T1TrackerHitAllocator.FreeSingle((T1BarrelCalHit*) aHit);
}

inline void T1BarrelCalHit::Add(G4double de, G4double dl) {
	fEdep += de; 
	fTrackLength += dl;
}

inline G4double T1BarrelCalHit::GetEdep() const { 
	return fEdep; 
}

inline G4double T1BarrelCalHit::GetTrackLength() const { 
	return fTrackLength; 
}

#endif