#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>
#include "ProtonPlan.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
class RunAction : public G4UserRunAction
{
  public:
    RunAction(G4String);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    struct OpticalParameters {
	G4double A0X;
	G4double A1X;
	G4double A2X;

	G4double A0Y;
	G4double A1Y;
	G4double A2Y;
    };
    struct BaseSpotParameters {
	G4ThreeVector Position;
	G4ThreeVector Direction;

    };
    
    void ComputeOpticalParameters();	
    void ComputeStartingEnergies();	
    void ComputeProtonsPerSpot(G4int);
    void ComputeBaseSpotParameters();
  //void GetStoppingPowers();
    void ComputeEnergySpread();

    G4double GetBaseStartingEnergy(size_t LayerNum) const {return SimStartEnergies[LayerNum];}
    G4double GetA0X(size_t LayerNum) const {return Optical[LayerNum].A0X;}
    G4double GetA1X(size_t LayerNum) const {return Optical[LayerNum].A1X;}
    G4double GetA2X(size_t LayerNum) const {return Optical[LayerNum].A2X;}

    G4double GetA0Y(size_t LayerNum) const {return Optical[LayerNum].A0Y;}
    G4double GetA1Y(size_t LayerNum) const {return Optical[LayerNum].A1Y;}
    G4double GetA2Y(size_t LayerNum) const {return Optical[LayerNum].A2Y;}

    G4double GetBaseX(size_t LayerNum, size_t SpotNum) const {return BaseSpot[LayerNum][SpotNum].Position.getX();}
    G4double GetBaseTheta(size_t LayerNum, size_t SpotNum) const {return BaseSpot[LayerNum][SpotNum].Direction.getX();}
    G4double GetBaseY(size_t LayerNum, size_t SpotNum) const {return BaseSpot[LayerNum][SpotNum].Position.getY();}
    G4double GetBasePhi(size_t LayerNum, size_t SpotNum) const {return BaseSpot[LayerNum][SpotNum].Direction.getY();}
    G4double GetBaseZ(G4int LayerNum, G4int SpotNum) const {return BaseSpot[LayerNum][SpotNum].Position.getZ();}
    G4double GetBasePsi(G4int LayerNum, G4int SpotNum) const {return BaseSpot[LayerNum][SpotNum].Direction.getZ();}
    G4double GetNozzleExit() const {return NozzleExit;}
    G4double GetEnergySpread(size_t LayerNum) const {return EnergySpreads[LayerNum];}

    size_t GetNoLayers() const {return NoProtons.size();}
    size_t GetNoSpots(size_t LayerNum) const {return NoProtons[LayerNum].size();}
    size_t GetNoProtons(size_t LayerNum, size_t SpotNum) const {return NoProtons[LayerNum][SpotNum];}
    void DecrementNoProtons(size_t LayerNum, size_t SpotNum) {--NoProtons[LayerNum][SpotNum];}
    ProtonPlan* GetThePlan() {return ThePlan;}

    void AddDose(G4int, G4int, G4int, G4double);
 private:
    static constexpr G4double NozzleExit=538.91;
    static constexpr G4double SADX=2203;
    static constexpr G4double SADY=1830;
    static constexpr G4double MeV2J=1.60217662e-13;
    G4double TotalNoProtons;
    G4double SimFraction;
    ProtonPlan* ThePlan;
    G4EmCalculator EmCalc;
	
    std::vector<G4double> SimStartEnergies;
    std::vector<G4double> StoppingPowers;
    std::vector<OpticalParameters> Optical;
    std::vector<std::vector<BaseSpotParameters>> BaseSpot;
    std::vector<std::vector<G4int>> NoProtons;
    std::vector<G4double> EnergySpreads;

    G4double NormFactor;
    G4double**** DoseSpectrum; 	//4D array; (x,y,z,D)
};

#endif
