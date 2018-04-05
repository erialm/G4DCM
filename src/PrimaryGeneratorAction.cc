#include "PrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* ARun)
: G4VUserPrimaryGeneratorAction(),
  TheParticleGun{nullptr}, TheRun{ARun}
{
	ThePlan=TheRun->GetThePlan();
	TheParticleGun = new G4ParticleGun;
  	TheParticleGun->SetParticleDefinition(G4Proton::ProtonDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  	delete TheParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::SampleSpotParameters(G4int LayerNumber, G4int SpotNumber)
{
	using std::exp;
	using std::sqrt;
	using std::pow;
	using std::tan;
	using std::cos;
	using std::tan;
	
	G4double GaussVal,Rand,B,MaxGauss;
	G4double A0X,A1X,A2X,X,Theta,NTheta,BaseTheta, BasePsi, NPsi;
	G4double GantryAngle=ThePlan->GetGantryAngle(LayerNumber);
	A0X=TheRun->GetA0X(LayerNumber);
	A1X=TheRun->GetA1X(LayerNumber);
	A2X=TheRun->GetA2X(LayerNumber);
	B=A0X*A2X-pow(A1X,2);
	MaxGauss=1/(2*Pi*sqrt(B));
	while (true)
	{
		X=CLHEP::RandFlat::shoot(-2*A2X,2*A2X);
		Theta=CLHEP::RandFlat::shoot(-2*A0X,2*A0X);
		GaussVal=1/(2*Pi*sqrt(B))*exp(-0.5*(A0X*pow(X,2)-2*A1X*X*Theta+A2X*pow(Theta,2))/B);
		Rand=CLHEP::RandFlat::shoot(0.,MaxGauss);
		if (Rand<GaussVal) 
		{
			SampledParameters.X=TheRun->GetBaseX(LayerNumber,SpotNumber)*mm+X*cos(GantryAngle);
			SampledParameters.Z=TheRun->GetBaseZ(LayerNumber,SpotNumber)*mm+X*sin(GantryAngle);
			BaseTheta=TheRun->GetBaseTheta(LayerNumber,SpotNumber);
			BasePsi=TheRun->GetBasePsi(LayerNumber,SpotNumber);	
			NTheta=(BaseTheta-tan(Theta*cos(GantryAngle)*1e-3))/(1+BaseTheta*tan(Theta*cos(GantryAngle)*1e-3));
			NPsi=(BasePsi-tan(Theta*sin(GantryAngle)*1e-3))/(1+BasePsi*tan(Theta*sin(GantryAngle)*1e-3));
			SampledParameters.Theta=NTheta;
			SampledParameters.Psi=NPsi;
			break;
		}
	}	


	G4double A0Y,A1Y,A2Y,Y,Phi,NPhi,BasePhi;
	A0Y=TheRun->GetA0Y(LayerNumber);
	A1Y=TheRun->GetA1Y(LayerNumber);
	A2Y=TheRun->GetA2Y(LayerNumber);
	B=A0Y*A2Y-pow(A1Y,2);
	MaxGauss=1/(2*Pi*sqrt(B));
	while (true) 
	{
		Y=CLHEP::RandFlat::shoot(-2*A2Y,2*A2Y);
		Phi=CLHEP::RandFlat::shoot(-2*A0Y,2*A0Y);
		GaussVal=1/(2*Pi*sqrt(B))*exp(-0.5*(A0Y*pow(Y,2)-2*A1Y*Y*Phi+A2Y*pow(Phi,2))/B);
		Rand=CLHEP::RandFlat::shoot(0.,MaxGauss);
		if (Rand<GaussVal) 
		{	
			SampledParameters.Y=TheRun->GetBaseY(LayerNumber,SpotNumber)*mm+Y*mm;
			BasePhi=TheRun->GetBasePhi(LayerNumber,SpotNumber);
			NPhi=(BasePhi-tan(Phi*1e-3))/(1+BasePhi*tan(Phi*1e-3));
			SampledParameters.Phi=NPhi;
			break;
		}
	}	
	
	G4double EnergySpread=TheRun->GetEnergySpread(LayerNumber);
	G4double NozzleEnergy=TheRun->GetBaseStartingEnergy(LayerNumber);
	G4double E=CLHEP::RandGauss::shoot(NozzleEnergy,EnergySpread);
	SampledParameters.Energy=E*MeV;
	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::array<G4int,2> PrimaryGeneratorAction::SampleSpot()
{
	G4double P;
	G4double R=CLHEP::RandFlat::shoot(0.,1.);
	std::array<G4int,2> LayerSpot;
	for (size_t i=0;i<ThePlan->GetCDFSize();++i)
	{
		P=ThePlan->GetCDFProb(i);
		if (R<=P)
		{
			LayerSpot[0]=ThePlan->GetCDFLayer(i);
			LayerSpot[1]=ThePlan->GetCDFSpot(i);
			break;
		}
	}
	return LayerSpot;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	std::array<G4int,2> LayerSpot=SampleSpot();
	SampleSpotParameters(LayerSpot[0], LayerSpot[1]);
	TheParticleGun->SetParticlePosition(G4ThreeVector(SampledParameters.X,SampledParameters.Y,SampledParameters.Z));
	TheParticleGun->SetParticleMomentumDirection((G4ThreeVector(SampledParameters.Theta,SampledParameters.Phi,SampledParameters.Psi)));
	TheParticleGun->SetParticleEnergy(SampledParameters.Energy);
	TheParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

