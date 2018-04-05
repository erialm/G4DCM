#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
StackingAction::StackingAction()
{
}
StackingAction::~StackingAction()
{
}
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
	//Kill all electrons outside the CT phantom
	if (aTrack->GetParticleDefinition()==G4Electron::ElectronDefinition() && aTrack->GetVolume()->GetName()=="World") return fKill;
	else return fWaiting;
	
}
void StackingAction::NewStage()
{
}
void StackingAction::PrepareNewEvent()
{
}
