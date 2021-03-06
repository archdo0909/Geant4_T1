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
// Tatsumi Koi for HandsOn4 of McGill Univ. Tutorial
//
#include "T1DetectorCalSDHit.hh"
#include "G4SystemOfUnits.hh"
#include "T1PrimaryGeneratorAction.hh"
#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
// HandsOn5: Draw box
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4ParticleGun.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4Allocator<BeamTestEmCalorimeterHit> BeamTestEmCalorimeterHitAllocator;

BeamTestEmCalorimeterHit::BeamTestEmCalorimeterHit()
	:fCellID(-1)
	,fDepositedEnergy(0)
	,fPosition()
	,fRotation()
	,pLogicalVolume(0)
{}

BeamTestEmCalorimeterHit::BeamTestEmCalorimeterHit(G4int id)
	:fCellID(id)
	,fDepositedEnergy(0)
	,fPosition()
	,fRotation()
	,pLogicalVolume(0)
{}

BeamTestEmCalorimeterHit::~BeamTestEmCalorimeterHit() {}

void BeamTestEmCalorimeterHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

	if(pVVisManager && (fDepositedEnergy>0.)) {

		// HandsOn5: Draw a box with depth propotional to the energy deposition
		G4double scale = BeamTestPrimaryGeneratorAction::Gun()->GetParticleEnergy();
		G4double depth = (50.*cm)*(fDepositedEnergy*MeV)/(scale*MeV);

		// Back face of box is flat against front face of calorimeter cell
		double z = fPosition.z()  + 25.*cm;
		G4ThreeVector myPos(fPosition.x(), fPosition.y(), z+depth);

		G4Transform3D trans(fRotation.inverse(), myPos);
		G4VisAttributes attribs;

		// Magenta with transparency
		G4Colour colour(1., 0., 1., 0.6);
		attribs.SetColour(colour);
		attribs.SetForceSolid(true);

		G4Box box("MyBox", 5.*cm, 5.*cm, depth);

		pVVisManager->Draw(box, attribs, trans);
	}
}

const std::map<G4String,G4AttDef>* BeamTestEmCalorimeterHit::GetAttDefs() const
{
	G4bool isNew;
	std::map<G4String,G4AttDef>* store 
		= G4AttDefStore::GetInstance("BeamTestEmCalorimeterHit",isNew);

	if (isNew) {
		G4String HitType("HitType");
		(*store)[HitType] = G4AttDef(HitType,"Hit Type", "Bookkeeping", "", "G4String");

		G4String ID("ID");
		(*store)[ID] = G4AttDef(ID, "ID", "Bookkeeping", "", "G4int");

		G4String Column("Column");
		(*store)[Column] = G4AttDef(Column, "Column ID", "Bookkeeping", "", "G4int");

		G4String Row("Row");
		(*store)[Row] = G4AttDef(Row, "Row ID", "Bookkeeping", "", "G4int");

		G4String Energy("Energy");
		(*store)[Energy] = G4AttDef(Energy, "Energy Deposited", "Physics", 
			"G4BestUnit", "G4double");

		G4String Pos("Pos");
		(*store)[Pos] = G4AttDef(Pos, "Position", "Physics", "G4BestUnit", "G4ThreeVector");

		G4String LVol("LVol");
		(*store)[LVol] = G4AttDef(LVol, "Logical Volume", "Bookkeeping", "", "G4String");
	}

	return store;
}

std::vector<G4AttValue>* BeamTestEmCalorimeterHit::CreateAttValues() const
{
	std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

	values->push_back(G4AttValue("HitType", "EmCalorimeterHit", ""));

	values->push_back
		(G4AttValue("ID", G4UIcommand::ConvertToString(fCellID), ""));

	values->push_back(G4AttValue("Column", " ", ""));

	values->push_back(G4AttValue("Row", " ", ""));

	values->push_back(G4AttValue("Energy", G4BestUnit(fDepositedEnergy, "Energy"), ""));

	values->push_back(G4AttValue("Pos", G4BestUnit(fPosition,"Length"), ""));

	if (pLogicalVolume) values->push_back(G4AttValue("LVol", pLogicalVolume->GetName(), ""));
	else values->push_back(G4AttValue("LVol", " ", ""));

	return values;
}

void BeamTestEmCalorimeterHit::Print()
{
	G4cout << "  Cell[" << fCellID << "] " << fDepositedEnergy/MeV << " (MeV)" << G4endl;


}
