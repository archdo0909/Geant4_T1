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

//T1CellParameterisation::T1CellParameterisation()
//{
//	// Initialise
//	G4int i(0);
//	G4int j(0);
//
//	for(i=0; i<128; i++) {                                                                                                                                              
//		G4int j = i / 64;
//		if(j != 1){
//		G4int column = i / 8;
//		G4int row = i % 8;
//		xCell.push_back((column-3)*10.*mm - 5*mm);
//		yCell.push_back(-(row-3)*10*mm + 5*mm);
//		zCell.push_back(-50.*mm);
//		}
//		else{
//			G4int colum_2 = (i - 64) / 8;
//			G4int row_2= (i - 64) % 8;
//			xCell.push_back((colum_2-3)*10.*mm - 5*mm);
//			yCell.push_back(-(row_2-3)*10*mm + 5*mm);
//			zCell.push_back(50.*mm);
//		}
//	}
//}

T1CellParameterisation::T1CellParameterisation(int nRow, int nColumn)///////go to Eventaction to check Nums 
{
	// Initialise
	int nNumCell = nRow * nColumn;

	for(int i=0; i<nNumCell; i++)
	{
		int column = i / nRow;
		int row = i % nRow;

		yCell.push_back(0.*mm);        ///5mm
		xCell.push_back((column-(nColumn/2-1))*10.*mm - 5.*mm);
		zCell.push_back(0.*mm);
	}


	/*for(int i=0; i<nNumCell; i++)
	{
	int column = i % nColumn;
	int row = i % nRow;

	xCell.push_back((row-(nRow/2-1))*10.*mm - 5.*mm);
	yCell.push_back((column-(nColumn/2-1))*10.*mm - 5.*mm);
	zCell.push_back(10.*mm);
	}*/

}

T1CellParameterisation::~T1CellParameterisation() {}

void T1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
	physVol->SetTranslation(G4ThreeVector(xCell[copyNo],yCell[copyNo],zCell[copyNo]));
}

