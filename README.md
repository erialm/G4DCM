# G4DCM
Will add more here at some point. 

The beam model is currently hard coded for Skandionkliniken. If I find it useful at some point I might add some messenger classes
so that some fiddling can more easily be done with text input. Currently, the necessary input files are all hardcoded as stored in ../../INPUTDATA/ relative the executable.

Multithreading is enabling if Geant4 is compiled with mulithreading support. The RunAction class exists as a single instance and 
stores the dose distribution. The RunAction::AddDose member function is protected by a mutex lock. The SteppingAction class is thread local,
and calls RunAction::AddDose for every step.

The StackingAction class ensures that no electrons are tracked outside the CT phantom volume.


