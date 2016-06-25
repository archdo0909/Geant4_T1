#include "T1ActionInitialization.hh"
#include "T1PrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "T1EventAction.hh"
#include "T1SteppingAction.hh"

T1ActionInitialization::T1ActionInitialization()
	:G4VUserActionInitialization()
{}

T1ActionInitialization::~T1ActionInitialization()
{}

void T1ActionInitialization::BuildForMaster()const
{
	SetUserAction(new B1RunAction);
}

void T1ActionInitialization::Build()const
{
	SetUserAction(new BeamTestPrimaryGeneratorAction);
	SetUserAction(new B1RunAction);

	BeamTestEventAction* eventAction = new BeamTestEventAction;
	SetUserAction(eventAction);

}