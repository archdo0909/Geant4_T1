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

G4Allocator<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitAllocator;

BeamTestSiliconMonitorHit::BeamTestSiliconMonitorHit()
	:fDepositedEnergy(0)
	,fCellID(-1)
{
   fIPD = 0;
   fIKEnergy = 0.0;
   fIPosition(0.0);
   fIMomentumD(0.0);
}
BeamTestSiliconMonitorHit::BeamTestSiliconMonitorHit(G4int id)
	:fCellID(id)
	,fDepositedEnergy(0)
{}

BeamTestSiliconMonitorHit::~BeamTestSiliconMonitorHit() {}



void BeamTestSiliconMonitorHit::Draw() {}



void BeamTestSiliconMonitorHit::Print()
{
   G4cout << "Incidence Particle Name and Kinetic Energy " << fIPD->GetParticleName() 
          << " " << fIKEnergy/MeV << " MeV" << G4endl; 
   G4cout << "Insidence position in silicon monitor " << fIPosition/mm << " in mm" << G4endl; 
   G4cout << "Incidence Direction " << fIMomentumD << G4endl; 
   G4cout << "  Cell[" << fCellID << "] " << fDepositedEnergy/MeV << " (MeV)" << G4endl;
}
