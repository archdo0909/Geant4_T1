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
// GEANT4 tag $Name:$
//
// T. Aso        Original author
// Jane Tinslay  Minor modifications for Jefferson lab tutorial   
// 
#ifndef BEAMTESTPHYSICSLIST_HH
#define BEAMTESTPHYSICSLIST_HH

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class BeamTestPhysicsList: public G4VUserPhysicsList {

public:

  // Constructor
  BeamTestPhysicsList();

  // Destructor
  virtual ~BeamTestPhysicsList();
  
protected:

  // Construct particles and physics processes
  void ConstructParticle();
  void ConstructProcess();
  
  void SetCuts();
  
private:
  
  // Helper methods
  void ConstructGeneral();
  void ConstructEM();

};

#endif



