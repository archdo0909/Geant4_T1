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
// GEANT4 tag $Name:$
//
// T. Aso        Original author
// 
#include "T1PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

G4ParticleGun* BeamTestPrimaryGeneratorAction::particleGun(0);

BeamTestPrimaryGeneratorAction::BeamTestPrimaryGeneratorAction()
{
	G4int n_particle = 1;
	particleGun  = new G4ParticleGun(n_particle);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	//G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	particleGun->SetParticleEnergy(15.*MeV);
	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.*cm));
}

BeamTestPrimaryGeneratorAction::~BeamTestPrimaryGeneratorAction()
{
	delete particleGun;
}

void BeamTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//particleGun->GeneratePrimaryVertex(anEvent);

	G4ThreeVector g1direction =  G4RandomDirection();
	G4ThreeVector g2direction = G4RandomDirection();

	particleGun->SetParticleEnergy(1.33*MeV);
	particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
	particleGun->SetParticleMomentumDirection(g1direction);
	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleEnergy(1.17*MeV);
	particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
	particleGun->SetParticleMomentumDirection(g2direction);
	particleGun->GeneratePrimaryVertex(anEvent);
}
