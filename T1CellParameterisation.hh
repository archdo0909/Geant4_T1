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
#ifndef T1CELLPARAMETERISATION_HH
#define T1CELLPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;

class T1CellParameterisation : public G4VPVParameterisation {

public:

	// Constructor
	//T1CellParameterisation();
	T1CellParameterisation(int nRow, int nColumn);


	// Destructor
	virtual ~T1CellParameterisation();

	// Method
	virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;

private:

	// Data members
	std::vector<G4double> xCell;
	std::vector<G4double> yCell;
	std::vector<G4double> zCell;

};

#endif


