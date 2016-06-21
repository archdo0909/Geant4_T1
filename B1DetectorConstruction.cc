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
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"
#include "G4PVReplica.hh"

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
  G4double cryst_dX = 10*mm, cryst_dY = 10*mm, cryst_dZ = 5*mm;
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
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
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
                     
  //
  // ring
  //
 // G4double gap = 0.5*mm;        //a gap for wrapping
 // G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
 // G4double dX = 80*mm, dY = 80*mm, dZ = 94*mm;

  //G4Box* solidRing =
	 // new G4Box("Ring", dX*8, dY/2, cryst_dZ/2);

  //G4LogicalVolume* logicRing =                         
	 // new G4LogicalVolume(solidRing,           //its solid
	 // default_mat,         //its material
	 // "Ring");             //its name

 // G4Box* solidBox =
//	  new G4Box("Container", dX, dY, dZ);

 // G4LogicalVolume* logicCon =
//	  new G4LogicalVolume(solidBox,			//its solid
	//					  default_mat,		//its material
	//					  "Container");		//its name

  //     
  // define crystal
  //
  
  //G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);

  //G4LogicalVolume* logicCryst = 
	 // new G4LogicalVolume(solidCryst,          //its solid
	 // cryst_mat,           //its material
	 // "CrystalLV");        //its name

  //logicCryst->SetSensitiveDetector(barrelCalSD);

//  G4Box* solidCryst = new G4Box("crystal", cryst_dX, cryst_dY, cryst_dZ);

//  G4LogicalVolume* logicCryst = 
//    new G4LogicalVolume(solidCryst,          //its solid
//                        cryst_mat,           //its material
//                       "CrystalLV");        //its name

//     logicCryst->SetSensitiveDetector(barrelCalSD);


  // place crystals within a ring 
  //

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

  //for (G4int icrys = 0; icrys < X_cryst ; icrys++) {
	 // for(G4int jcrys = 0; jcrys < Y_cryst; jcrys++){
	 // G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
	 // G4ThreeVector position = (cryst_dZ)*uz*icrys - 0.5*cryst_dZ*uz*nb_cryst;
	 // G4Transform3D transform = G4Transform3D(rotm,position);

	 // new G4PVPlacement(transform,             //rotation,position
		//  logicCryst,            //its logical volume
		//  "crystal",             //its name
		//  logicCon,             //its mother  volume
		//  false,                 //no boolean operation
		//  (icrys + jcrys),                 //copy number
		//  checkOverlaps);       // checking overlaps 
	 // }

  G4VSolid* hadCalorimeterSolid = new G4Box("HadCalorimeterBox",100*mm,80.*mm,80.*mm); 
  G4LogicalVolume* hadCalorimeterLogical = 
	  new G4LogicalVolume(hadCalorimeterSolid,
	                      default_mat,
						  "HadCalorimeterLogical");

  new G4PVPlacement(0,                                    //no rotation
	                G4ThreeVector(0.,0.,3.*cm),            //at (0,0,0)
	                hadCalorimeterLogical,                //its logical volume
					"HadCalorimeterPhysical",             //its name
					0,                                    //its mother  volume
					false,                                //no boolean operation
					0,                                    //copy number
					checkOverlaps);                       // checking overlaps 

  // hadron calorimeter column 
  G4VSolid* HadCalColumnSolid = new G4Box("HadCalColumnBox",5.*mm,80.*mm,80.*mm); 
  G4LogicalVolume* HadCalColumnLogical = 
	  new G4LogicalVolume(HadCalColumnSolid,
	                      default_mat,
	                      "HadCalColumnLogical"); 
  new G4PVReplica("HadCalColumnPhysical",
				  HadCalColumnLogical,
			      hadCalorimeterLogical,
		     	  kXAxis,
				  20,
				  5.*mm); 

  // hadron calorimeter cell 
  G4VSolid* HadCalCellSolid 
	  = new G4Box("HadCalCellBox",5.*mm,10.*mm,80.*cm); 
  G4LogicalVolume* HadCalCellLogical
	  = new G4LogicalVolume(HadCalCellSolid,
	                      default_mat,
						  "HadCalCellLogical"); 

  new G4PVReplica("HadCalCellPhysical",
					HadCalCellLogical,
					HadCalColumnLogical,
					kYAxis,
					8,
		     		10*mm); 

  // hadron calorimeter layers 
  G4VSolid* HadCalLayerSolid = new G4Box("HadCalLayerBox",5.*mm,10.*mm, 10.*mm); 
  G4LogicalVolume* HadCalLayerLogical
	  = new G4LogicalVolume(HadCalLayerSolid,
	                        cryst_mat,
	                        "HadCalLayerLogical"); 

  HadCalLayerLogical->SetSensitiveDetector(barrelCalSD);

  new G4PVPlacement(0,
	                G4ThreeVector(0., 0., 3.*cm),
					HadCalLayerLogical,
					"Calorimeter",
					HadCalCellLogical,
					false,
					8,
					checkOverlaps);

//  fScoringVolume = HadCalLayerLogical;
 
  //G4VSolid* HadCalLayerSolid = new G4Box("HadCalLayerBox",5.*mm,10.*mm, 10.*mm); 
  //G4LogicalVolume* HadCalLayerLogical = 
	 // new G4LogicalVolume(HadCalLayerSolid,
	 // cryst_mat,
	 // "HadCalLayerLogical"); 
  //new G4PVReplica("HadCalLayerPhysical",
	 // HadCalLayerLogical,
	 // HadCalCellLogical,
	 // kZAxis,
	 // 8,
	 // 10.*mm); 

  //// scintillator plates 
  //G4VSolid* HadCalScintiSolid = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm); 
  //fHadCalScintiLogical = 
  //    new G4LogicalVolume(HadCalScintiSolid,
		//				  scintillator,
		//				  "HadCalScintiLogical"); 
  //new G4PVPlacement(0,
	 //               G4ThreeVector(0.,0.,2.*cm),
	 //                fHadCalScintiLogical,
	 //               "HadCalScintiPhysical",
	 //                HadCalLayerLogical,
	 //                false,
	 //                0,
	 //                checkOverlaps);

  //
  // place detector in world
  //                    



  //new G4PVPlacement(0,                       //no rotation
	 // G4ThreeVector(0,0,15*cm),         //at (0,0,0)
	 // logicRing,               //its logical volume
	 // "Detector",              //its name
	 // logicWorld,              //its mother  volume
	 // false,                   //no boolean operation
	 // 0,                       //copy number
	 // checkOverlaps);         // checking overlaps 

  //new G4PVPlacement(0,                       //no rotation
  //                  pos2,                    //at position
  //                  logicShape2,             //its logical volume
  //                  "Shape2",                //its name
  //                  logicEnv,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking


  //    
  // Shape 1
  //  
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  //G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
  //      
  //// Conical section shape       
  //G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  //G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  //G4double shape1_hz = 3.*cm;
  //G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  //G4Cons* solidShape1 =    
  //  new G4Cons("Shape1", 
  //  shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //  shape1_phimin, shape1_phimax);
  //                    
  //G4LogicalVolume* logicShape1 =                         
  //  new G4LogicalVolume(solidShape1,         //its solid
  //                      shape1_mat,          //its material
  //                      "Shape1");           //its name
  //             
  //new G4PVPlacement(0,                       //no rotation
  //                  pos1,                    //at position
  //                  logicShape1,             //its logical volume
  //                  "Shape1",                //its name
  //                  logicEnv,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking


  ////     
  //// Shape 2
  ////
  //G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  //G4ThreeVector pos2 = G4ThreeVector(0, -100*cm, 7*cm);

  //// Trapezoid shape       
  //G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  //G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  //G4double shape2_dz  = 6*cm;      
  //G4Trd* solidShape2 =    
  //  new G4Trd("Shape2",                      //its name
  //            0.5*shape2_dxa, 0.5*shape2_dxb, 
  //            0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
  //              
  //G4LogicalVolume* logicShape2 =                         
  //  new G4LogicalVolume(solidShape2,         //its solid
  //                      shape2_mat,          //its material
  //                      "Shape2");           //its name
  //             
  //new G4PVPlacement(0,                       //no rotation
  //                  pos2,                    //at position
  //                  logicShape2,             //its logical volume
  //                  "Shape2",                //its name
  //                  logicEnv,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking
                
  // Set Shape2 as scoring volume
  //
  //fScoringVolume = logicCryst;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

