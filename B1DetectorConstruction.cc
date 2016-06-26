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
#include "T1DetectorCalSD.hh"
#include "T1CellParameterisation.hh"

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
#include "G4PVParameterised.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"

#include <iostream>
#include <string>
#define NUM_CRYSTAL 64


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//T1original

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fSolidWorld(0),
  fLogicWorld(0)
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
	G4double a, z;
	G4double density;
	G4int nel;

	//Air
	G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);

	G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
	Air->AddElement(N, 70*perCent);
	Air->AddElement(O, 30*perCent);

	//Lead
	G4Material* Pb = 
		new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

	// Print all the materials defined.
	//
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//--------- Sizes of the principal geometrical components (solids)  ---------

	fNbOfChambers = 5;
	fChamberWidth = 20*cm;
	fChamberSpacing = 80*cm;

	fTrackerLength = 100.0 * mm;						  //FULL Length of detector Length
	fTargetLength  = 1. * m;			                  // Full length of Target

	fTargetMater  = Air;

	fWorldLength= (fTargetLength+fTrackerLength);

	G4double targetSize  = 0.5*fTargetLength;    // Half length of the Target  
	G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	//------------------------------ 
	// World
	//------------------------------ 

	//fTrackerLength = 100mm + fTargetLength = 1.0m == 110cm
	G4double HalfWorldLength = 0.5*fWorldLength;

	//fSolidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	fSolidWorld= new G4Box("world",7*cm,7*cm,7*cm);
	fLogicWorld= new G4LogicalVolume( fSolidWorld,			//solid
	                                  Air,				    //Material
									  "World",				//name
									  0,					//FieldManager
									  0,					//Sensitive Detector
									  0);					//optimize

	//  Must place the World Physical volume unrotated at (0,0,0).
	// 
	fPhysiWorld = new G4PVPlacement(0,               // no rotation
		G4ThreeVector(), // at (0,0,0)
		fLogicWorld,      // its logical volume
		"World",         // its name
		0,               // its mother  volume
		false,           // no boolean operations
		0);              // copy number

	////------------------------------ 
	//// Target
	////------------------------------

	//G4ThreeVector positionTarget = G4ThreeVector(0,0,0);

	//fSolidTarget = new G4Box("target",targetSize,targetSize,targetSize);
	//fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target",0,0,0);
	//fPhysiTarget = new G4PVPlacement(0,               // no rotation
	//	positionTarget,  // at (x,y,z)
	//	fLogicTarget,     // its logical volume
	//	"Target",        // its name
	//	fLogicWorld,      // its mother  volume
	//	false,           // no boolean operations
	//	0);              // copy number 

	//G4cout << "Target is " << fTargetLength/cm << " cm of " 
	//	<< fTargetMater->GetName() << G4endl;
#if 0
	////------------------------------ 
	//// Tracker
	////------------------------------

	//G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

	//fSolidTracker = new G4Box("tracker",trackerSize,trackerSize,trackerSize);
	//fLogicTracker = new G4LogicalVolume(fSolidTracker , Air, "Tracker",0,0,0);  
	//fPhysiTracker = new G4PVPlacement(0,              // no rotation
	//	positionTracker, // at (x,y,z)
	//	fLogicTracker,    // its logical volume
	//	"Tracker",       // its name
	//	fLogicWorld,      // its mother  volume
	//	false,           // no boolean operations
	//	0);              // copy number 

	////------------------------------ 
	//// Tracker segments
	////------------------------------
	//// 
	//// An example of Parameterised volumes
	//// dummy values for G4Box -- modified by parameterised volume

	//fSolidChamber = new G4Box("chamber", 100*cm, 100*cm, 10*cm); 
	//fLogicChamber = new G4LogicalVolume(fSolidChamber, fChamberMater,"Chamber",0,0,0);

	//G4double firstPosition = -trackerSize + 0.5*fChamberWidth;
	////G4double firstPosition = 0.;
	//G4double firstLength = fTrackerLength/10;
	//G4double lastLength  = fTrackerLength;

	//G4VPVParameterisation* chamberParam = new ExP01ChamberParameterisation(  
	//	fNbOfChambers,          // NoChambers 
	//	firstPosition,         // Z of center of first 
	//	fChamberSpacing,        // Z spacing of centers
	//	fChamberWidth,          // Width Chamber 
	//	firstLength,           // lengthInitial 
	//	lastLength);           // lengthFinal

	//// dummy value : kZAxis -- modified by parameterised volume
	////
	//fPhysiChamber = new G4PVParameterised(
	//	"Chamber",       // their name
	//	fLogicChamber,    // their logical volume
	//	fLogicTracker,    // Mother logical volume
	//	kZAxis,          // Are placed along this axis 
	//	fNbOfChambers,    // Number of chambers
	//	chamberParam);   // The parametrisation

	//G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
	//	<< "The chambers are " << fChamberWidth/mm << " mm of " 
	//	<< fChamberMater->GetName() << "\n The distance between chamber is "
	//	<< fChamberSpacing/cm << " cm" << G4endl;

	//------------------------------------------------ 
	// Sensitive detectors
	//------------------------------------------------ 

//	G4SDManager* SDman = G4SDManager::GetSDMpointer();

//	G4String trackerChamberSDname = "TrackerChamberSD";
//	T1BarrelCalSD* aTrackerSD = new T1BarrelCalSD( trackerChamberSDname );
//	SDman->AddNewDetector( aTrackerSD );
//	fLogicChamber->SetSensitiveDetector( aTrackerSD );

	//--------- Visualization attributes -------------------------------

//	G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
//	fLogicWorld  ->SetVisAttributes(BoxVisAtt);  
//	fLogicTarget ->SetVisAttributes(BoxVisAtt);
//	fLogicTracker->SetVisAttributes(BoxVisAtt);

//	G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
//	fLogicChamber->SetVisAttributes(ChamberVisAtt);

	//--------- example of User Limits -------------------------------

	// below is an example of how to set tracking constraints in a given
	// logical volume(see also in N02PhysicsList how to setup the processes
	// G4StepLimiter or G4UserSpecialCuts).

	// Sets a max Step length in the tracker region, with G4StepLimiter
	//
//	G4double maxStep = 0.5*fChamberWidth; 
//	fLogicTracker->SetUserLimits(new G4UserLimits(maxStep));

	// Set additional contraints on the track, with G4UserSpecialCuts
	//
	// G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
	// logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
	//                                               minEkin));
#endif

	/////////////////////////////////////////////////////////////////
	//Mother Volume
	G4VSolid* calorimeterSolid = new G4Box("Calorimeter_Solid", // Name
											4.*cm,                // x half length
											4.*cm,                // y half length
											50.*mm) ;             // z half length
	////////////if there is any change, Go to parameterisation also /////////

	G4LogicalVolume* calorimeterLogical =
		new G4LogicalVolume(calorimeterSolid,       // Solid
							Air,                    // Material
							"Calorimeter_Logical"); // Name

	new G4PVPlacement(0,                          // Rotation matrix pointer
		              G4ThreeVector(0.,0.,10.*cm), // Translation vector
					  calorimeterLogical,         // Logical volume
					  "Calorimeter_Physical",     // Name
					  fLogicWorld,             // Mother volume
		              false,                      // Unused boolean
		              0);                         // Copy number     
   


	G4NistManager* nist = G4NistManager::Instance();
	G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");

	G4VSolid* cellSolid = new G4Box("Cell_Solid", // Name
									5.*mm,         // x half length
									5.*mm,         // y half length
									5.*mm);      // z half length

	G4LogicalVolume* cellLogical
		= new G4LogicalVolume(cellSolid,       // Solid
							  cryst_mat,       // Material
							  "Cell_Logical"); // Name

	G4VPVParameterisation* cellParam = new T1CellParameterisation(); 

	new G4PVParameterised("Cell_Physical",    // Name
		                  cellLogical,        // Logical volume
						  calorimeterLogical, // Mother volume
						  kXAxis,             // Axis    
					      128,                // Number of replicas
					      cellParam);         // Parameterisation

	//G4VSolid* cellSolid2 = new G4Box("Cell_Solid", // Name
	//	                             5.*mm,         // x half length
	//	                             5.*mm,         // y half length
	//	                             2.5*mm);      // z half length
	//G4LogicalVolume* celllogical2
	//	=new G4LogicalVolume(cellSolid2,
	//	                     cryst_mat,
	//						 "Cell_Logical2");

	//G4VPVParameterisation* cellParam2 = new T1CellParameterisation(100.0);

	//new G4PVParameterised("Cell_Physical",    // Name
	//					  celllogical2,        // Logical volume
	//					  calorimeterLogical, // Mother volume
	//					  kXAxis,             // Axis    
	//				      64,                // Number of replicas
	//					  cellParam2);         // Parameterisation


#if 0
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
  G4int X_cryst = 2;
  G4int Y_cryst = 2;
  G4int nb_cryst = X_cryst * Y_cryst;
  G4int nb_rings = 9;
  fNbOfChambers = X_cryst * Y_cryst;
  //
  //G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  //G4double cosdPhi = std::cos(half_dPhi);
  //G4double tandPhi = std::tan(half_dPhi);
  // 
  //G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  //G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;

  G4SDManager* SDman= G4SDManager::GetSDMpointer();
  G4String SDname="barrelCal";
  
  T1BarrelCalSD* barrelCalSD = new T1BarrelCalSD(SDname);
  SDman-> AddNewDetector(barrelCalSD);

  // for barrel calorimeter
  //T1BarrelCalSD* barrelCalSD[NUM_CRYSTAL] = { NULL };
   
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
	  G4ThreeVector(0., 0., 15*cm),
	  calorLV,
	  "Calorimeter",
	  logicWorld,
	  false,
	  0,
	  checkOverlaps);
      
  //! Crystalの形状とサイズを設定する
  G4Box* solidCryst = new G4Box("crystal", 10*mm, 10*mm, 5*mm);

  //! Crystalのボリュームを定義する
  //G4LogicalVolume* logicCryst[NUM_CRYSTAL] = { NULL };
  G4LogicalVolume* logicCryst = new G4LogicalVolume(solidCryst, cryst_mat, "CrystalLV");

  //! Crystalに検出器を定義する
  logicCryst->SetSensitiveDetector(barrelCalSD);

 // while (fNbOfChambers != 0){
   for(int icrys = 0; icrys < X_cryst; icrys++)
   {
	   for(int jcrys = 0; jcrys < Y_cryst; jcrys++)
	   {
		   G4double X = icrys * 10.0;
		   G4double Y = jcrys * 10.0;

		   //G4String SDNoName="BarrelCal" + std::to_string((long double)icrys+jcrys);

		   //barrelCalSD[icrys+jcrys] = new T1BarrelCalSD(SDNoName);
		   //SDman-> AddNewDetector(barrelCalSD[icrys+jcrys]);

		   //logicCryst[icrys+jcrys] = new G4LogicalVolume(solidCryst, cryst_mat, "CrystalLV");
		   //logicCryst[icrys+jcrys]->SetSensitiveDetector(barrelCalSD[icrys+jcrys]);

		  // G4ThreeVector position = (icrys * 10., jcrys * 10., 0.);
		   G4ThreeVector position = G4ThreeVector(icrys * 10.,jcrys * 10.,0.);
		   G4VPhysicalVolume* detector =
		        new G4PVPlacement(
			                       0,
					       	   	   position,
							       logicCryst,
							       "crystal",
							       calorLV,
							       false,
							       0,
							       checkOverlaps);
		  // fNbOfChambers--;
	   }
	 }
  //}
#endif
   ////////////////////////////////////////////////////////////////////////
   // HandsOn4: Defining sensitive detector
   // Create a new T1CellParameterisation sensitive detector
   G4VSensitiveDetector* detector = new BeamTestEmCalorimeter("Calorimeter"); 

   // Get pointer to detector manager
   G4SDManager* SDman = G4SDManager::GetSDMpointer();

   // Register detector with manager
   SDman->AddNewDetector(detector);

   // Attach detector to volume defining calorimeter cells
   cellLogical->SetSensitiveDetector(detector);

   ////////////////////////////////////////////////////////////////////////
   // Visualisation attributes

   // Invisible world volume.
   //fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);

   // HandsOn4: Calorimeter attributes 
   // Invisible calorimeter mother volume
   //calorimeterLogical->SetVisAttributes(G4VisAttributes::Invisible);

   /* G4VisAttributes* calorimeterMotherAttributes =
   new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 0.1));
   calorimeterMotherAttributes->SetForceSolid(true);
   calorimeterLogical->SetVisAttributes(calorimeterMotherAttributes);*/

   // Calorimeter cells - green with transparency
   G4VisAttributes* calorimeterAttributes =
	   new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.1));
   calorimeterAttributes->SetForceSolid(true);
   cellLogical->SetVisAttributes(calorimeterAttributes);

   return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::ConstructSDandField()
{

}