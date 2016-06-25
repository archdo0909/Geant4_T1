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
// Jane Tinslay - adapted from A01 example
//
#include "T1CellParameterisation.hh"

#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

T1CellParameterisation::T1CellParameterisation()
{
	// Initialise
	G4int i(0);

	for(i=0; i<64; i++) {
		G4int column = i / 8;
		G4int row = i % 8;
		xCell.push_back((column-3)*10.*mm - 5*mm);
		yCell.push_back((row-3)*10*mm - 5*mm);
	}
}

T1CellParameterisation::~T1CellParameterisation() {}

void T1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
	physVol->SetTranslation(G4ThreeVector(xCell[copyNo],yCell[copyNo],0.));
}

