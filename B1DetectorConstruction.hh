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
// $Id: B1DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExP01DetectorMessenger;

/// Detector construction class to define materials and geometry.

//class B1DetectorConstruction : public G4VUserDetectorConstruction
//{
//  public:
//    B1DetectorConstruction();
//    virtual ~B1DetectorConstruction();
//
//    virtual G4VPhysicalVolume* Construct();
//	virtual void ConstructSDandField();
//    
//    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
//
//	void DefineMaterials();
//
//  protected:
//    G4LogicalVolume*  fScoringVolume;
//
//private:
//	G4Box*             fSolidWorld;    // pointer to the solid envelope 
//	G4LogicalVolume*   fLogicWorld;    // pointer to the logical envelope
//	G4VPhysicalVolume* fPhysiWorld;    // pointer to the physical envelope
//
//	G4Box*             fSolidTarget;   // pointer to the solid Target
//	G4LogicalVolume*   fLogicTarget;   // pointer to the logical Target
//	G4VPhysicalVolume* fPhysiTarget;   // pointer to the physical Target
//
//	G4Box*             fSolidTracker;  // pointer to the solid Tracker
//	G4LogicalVolume*   fLogicTracker;  // pointer to the logical Tracker
//	G4VPhysicalVolume* fPhysiTracker;  // pointer to the physical Tracker
//
//	G4Box*             fSolidChamber;  // pointer to the solid Chamber
//	G4LogicalVolume*   fLogicChamber;  // pointer to the logical Chamber
//	G4VPhysicalVolume* fPhysiChamber;  // pointer to the physical Chamber
//
//	G4Material*         fTargetMater;  // pointer to the target  material
//	G4Material*         fChamberMater; // pointer to the chamber material
//
//	G4double fWorldLength;            // Full length of the world volume
//	G4double fTargetLength;           // Full length of Target
//	G4double fTrackerLength;          // Full length of Tracker
//	G4int    fNbOfChambers;            // Nb of chambers in the tracker region
//	G4double fChamberWidth;            // width of the chambers
//	G4double fChamberSpacing;          // distance between chambers
//};
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class B1DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	B1DetectorConstruction();
	virtual ~B1DetectorConstruction();

	virtual G4VPhysicalVolume* Construct();

	G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

protected:
	G4LogicalVolume*  fScoringVolume;

private:
	void constructVessels(G4LogicalVolume* pMotherVolume, G4Material* pMaterial);
	void constructWalls(G4LogicalVolume* pMotherVolume);
	void constructCalorimeter(G4LogicalVolume* pMotherVolume, G4Material* pMaterial, int nNum);    

	G4Material* makeVesselMaterial(void);
	G4Material* makePassingMaterial(void);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#endif

