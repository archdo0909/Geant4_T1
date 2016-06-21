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
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

static const double     pi  = 3.14159265358979323846;
static const double  twopi  = 2*pi;

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4MultiFunctionalDetector.hh"
//#include "G4Cons.hh"
//#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"

#include "T1BarrelCalSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//T1original

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{
	DefineMaterials();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::DefineMaterials()
{
	G4NistManager* man = G4NistManager::Instance();

	G4bool isotopes = false;

	G4Element*  O = man->FindOrBuildElement("O" , isotopes); 
	G4Element* Si = man->FindOrBuildElement("Si", isotopes);
	G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);  

	G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
	LSO->AddElement(Lu, 2);
	LSO->AddElement(Si, 1);
	LSO->AddElement(O , 5);  
}

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  //Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");

  // Gamma detector Parameters
  //
  G4double cryst_dX = 6*cm, cryst_dY = 6*cm, cryst_dZ = 3*cm;
  G4int X_cryst = 8;
  G4int Y_cryst = 8;
  G4int nb_cryst = X_cryst * Y_cryst;
  G4int nb_rings = 9;
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  // 
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;

  G4SDManager* SDman= G4SDManager::GetSDMpointer();
  G4String SDname="/barrelCal";

  // for barrel calorimeter
  T1BarrelCalSD* barrelCalSD= new T1BarrelCalSD(SDname);
  SDman-> AddNewDetector(barrelCalSD);
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  //G4double world_sizeXY = 1.2*env_sizeXY;
  //G4double world_sizeZ  = 1.2*env_sizeZ;
  G4double world_sizeXY = env_sizeXY;
  G4double world_sizeZ = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //Mother Volume
  G4VSolid* calorimeterS
	  =new G4Box("Calorimeter",
				 10*cm, 10*cm, 5*cm);

  G4LogicalVolume* calorLV
	  = new G4LogicalVolume(calorimeterS,
							default_mat,
							"Calorimeter");

	  new G4PVPlacement(  0,
	  G4ThreeVector(),
	  calorLV,
	  "Calorimeter",
	  logicWorld,
	  false,
	  0,
	  checkOverlaps);
      

  G4Box* solidCryst = new G4Box("crystal", 10*mm, 10*mm, 5*mm);

  G4LogicalVolume* logicCryst =
	  new G4LogicalVolume(solidCryst,
						  cryst_mat,
						  "CrystalLV");

   for(G4int icrys = 0; icrys < X_cryst; icrys++)
   {
	   for(G4int jcrys = 0; jcrys < Y_cryst; jcrys++)
	   {
		   G4double X = icrys * 10.0;
		   G4double Y = jcrys * 10.0;

		  // G4ThreeVector position = (icrys * 10., jcrys * 10., 0.);
		   G4ThreeVector position = G4ThreeVector(icrys * 10.,jcrys * 10.,0.);
		   G4int c_num = (icrys + 1)*(jcrys + 1);
		   new G4PVPlacement(
			                 0,
							 position,
							 logicCryst,
							 "crystal",
							 calorLV,
							 false,
							 c_num,
							 checkOverlaps);
	   }
   }
 G4Tubs* soildDetector =
	 new G4Tubs("Detector", )
                     
  ////
  //// ring
  ////
  //G4double gap = 0.5*mm;
  //G4double dX = cryst_dX - gap, dY = cryst_dY - gap;

  //G4Box* solidRing =
	 // new G4Box("Ring", dX*8, dY/2, cryst_dZ/2);

  //G4LogicalVolume* logicRing =                         
	 // new G4LogicalVolume(solidRing,           //its solid
	 // default_mat,         //its material
	 // "Ring");             //its name

  ////     
  //// define crystal
  ////
  //
  //G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);

  //G4LogicalVolume* logicCryst = 
	 // new G4LogicalVolume(solidCryst,          //its solid
	 // cryst_mat,           //its material
	 // "CrystalLV");        //its name

  //logicCryst->SetSensitiveDetector(barrelCalSD);

  //// place crystals within a ring 
  ////

  //for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
	 // G4double phi = 0;
	 // G4RotationMatrix rotm  = G4RotationMatrix();
	 // rotm.rotateY(0); 
	 // rotm.rotateZ(phi);
	 // G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
	 // G4ThreeVector position = (cryst_dZ)*uz*icrys - 0.5*cryst_dZ*uz*nb_cryst;
	 // G4Transform3D transform = G4Transform3D(rotm,position);

	 // new G4PVPlacement(transform,             //rotation,position
		//  logicCryst,            //its logical volume
		//  "crystal",             //its name
		//  logicRing,             //its mother  volume
		//  false,                 //no boolean operation
		//  icrys,                 //copy number
		//  checkOverlaps);       // checking overlaps 
  //}
 
  ////
  //// place detector in world
  ////                    

  //new G4PVPlacement(0,                       //no rotation
	 // G4ThreeVector(0,0,5*cm),         //at (0,0,0)
	 // logicRing,               //its logical volume
	 // "Detector",              //its name
	 // logicWorld,              //its mother  volume
	 // false,                   //no boolean operation
	 // 0,                       //copy number
	 // checkOverlaps);         // checking overlaps 

  //fScoringVolume = logicCryst;

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::ConstructSDandField()
{

}