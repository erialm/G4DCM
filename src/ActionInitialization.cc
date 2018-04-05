
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(G4String FilePath)
 : G4VUserActionInitialization(), PlanPath{FilePath}
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void ActionInitialization::BuildForMaster() const
{
  RunAction* runAction = new RunAction;
  SetUserAction(runAction);
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  RunAction* TheRun = new RunAction(PlanPath);
  SetUserAction(TheRun); 
  SetUserAction(new StackingAction);
  SetUserAction(new SteppingAction(TheRun));
  SetUserAction(new PrimaryGeneratorAction(TheRun));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
