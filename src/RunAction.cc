#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DicomDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <fstream>
#include "c2_factory.hh"
#include <cmath>
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4ParticleDefinition.hh"
#include <stdexcept>
static c2_factory<G4double> c2; 
typedef c2_ptr<G4double> c2p;
#define BINS 3
#define PI 3.14159265359
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction(G4String PlanPath)
: G4UserRunAction(), ThePlan{nullptr}, TotalNoProtons{0.}, SimFraction{0.}
{
	ThePlan=new ProtonPlan(PlanPath);
	
	size_t NoLayers=ThePlan->GetNoLayers();	
	BaseSpot.resize(NoLayers);
	NoProtons.resize(NoLayers);
	for (size_t i=0;i<NoLayers;++i)
	{
		BaseSpot[i].resize(ThePlan->GetNoSpots(i));	
		NoProtons[i].resize(ThePlan->GetNoSpots(i));
	}
	ComputeStartingEnergies();	
	ComputeEnergySpread();
	ComputeOpticalParameters();	
	ComputeBaseSpotParameters();
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::ComputeOpticalParameters()
{

	//THIS USES INTERPOLATION OF THE PARAMETERS, RATHER THAN FITTING; CHANGE AT SOME POINT OR KEEP?
	using std::vector;
	std::ifstream ParametersFile{"../../INPUTDATA/MergedPara.txt"};
	if (!ParametersFile) throw std::runtime_error("Could not open optical parameters data file!");
	G4double EEnt,A0XEnt,A1XEnt,A2XEnt,A0YEnt,A1YEnt,A2YEnt;
	vector<G4double> Energy;
	vector<G4double> A0X;	//A0 = variance of angular distribution
	vector<G4double> A1X;	//A1 = correlation of angular and spatial distributions
	vector<G4double> A2X;	//A2 = variance of spatial distribution
	vector<G4double> A0Y;
	vector<G4double> A1Y;
	vector<G4double> A2Y;

	while(ParametersFile>>EEnt>>A0XEnt>>A1XEnt>>A2XEnt>>A0YEnt>>A1YEnt>>A2YEnt)
	{
		Energy.push_back(EEnt);		//MeV at nozzle exit
		A0X.push_back(A0XEnt*1e6);	//mrad^2
		A1X.push_back(A1XEnt*1e3);	//mm*mrad
		A2X.push_back(A2XEnt);		//mm	
		A0Y.push_back(A0YEnt*1e6);
		A1Y.push_back(A1YEnt*1e3);
		A2Y.push_back(A2YEnt);
		
	}

	c2p A0XInt=c2.interpolating_function().load(Energy, A0X,true,0,true,0,false);//false=linear interpolation
	c2p A1XInt=c2.interpolating_function().load(Energy, A1X,true,0,true,0,false);
	c2p A2XInt=c2.interpolating_function().load(Energy, A2X,true,0,true,0,false);
	c2p A0YInt=c2.interpolating_function().load(Energy, A0Y,true,0,true,0,false);
	c2p A1YInt=c2.interpolating_function().load(Energy, A1Y,true,0,true,0,false);
	c2p A2YInt=c2.interpolating_function().load(Energy, A2Y,true,0,true,0,false);


	OpticalParameters Entry;
	G4double LayerEnergy;
	for (size_t i=0;i<SimStartEnergies.size();++i)
	{
		LayerEnergy=SimStartEnergies[i];	
		Entry.A0X=A0XInt(LayerEnergy);
		Entry.A1X=A1XInt(LayerEnergy);
		Entry.A2X=A2XInt(LayerEnergy);
		Entry.A0Y=A0YInt(LayerEnergy);
		Entry.A1Y=A1YInt(LayerEnergy);
		Entry.A2Y=A2YInt(LayerEnergy);
		Optical.push_back(Entry);	
	}
	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::ComputeStartingEnergies()
{
	std::ifstream EnergyFile{"../../INPUTDATA/RoundedBestEnergies.txt"};
	if (!EnergyFile) throw std::runtime_error("Could not open energies data file!");
	
	std::vector<G4double> NozzleEnergy;
	std::vector<G4double> IsoEnergy;
	G4double EnergyEntry;
	while (EnergyFile>>EnergyEntry)
	{
		NozzleEnergy.push_back(EnergyEntry);
		IsoEnergy.push_back(std::floor(EnergyEntry));
	}

	c2p Iso2NozzleExit=c2.interpolating_function().load(IsoEnergy,NozzleEnergy,true,0,true,0,true);
	for (size_t i=0;i<ThePlan->GetNoLayers();++i) SimStartEnergies.push_back(Iso2NozzleExit(ThePlan->GetLayerEnergy(i)));

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::ComputeEnergySpread()
{
	using std::pow;
	G4double E,EnergySpread;
	for (size_t i=0;i<SimStartEnergies.size();++i)
	{
		E=SimStartEnergies[i];
		EnergySpread=-3.0147e-07*pow(E,3)+7.5117e-05*pow(E,2)+1.2467e-03*E+2.0982e-01;
		EnergySpreads.push_back(EnergySpread);
	}
	
		
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::ComputeBaseSpotParameters()
{
	//FIXED COORDINATE SYSTEM: LOOKING TOWARD THE GANTRY, AT GANTRY ANGLE 0, X AXIS POINTS TO THE RIGHT, Y TOWARD THE GANTRY AND Z UPWARDS IN THE OPPOSITE DIRECTION OF THE BEAM	
	//GANTRY ROTATES CLOCKWISE
	using std::cos;
	using std::sin;
	using std::abs;
	G4double GantryAngle, KY, Y;
	G4ThreeVector SourceX, IsoX, DiffX;
	for (size_t i=0;i<ThePlan->GetNoLayers();++i)
	{
		GantryAngle=ThePlan->GetGantryAngle(i)*PI/180; //convert to radians
		SourceX=G4ThreeVector{SADX*sin(GantryAngle),0,SADX*cos(GantryAngle)}; //Position at virtual source position X in fixed coordinates
		for (size_t j=0;j<ThePlan->GetNoSpots(i);++j)
		{					
			IsoX=G4ThreeVector{ThePlan->GetSpotX(i,j)*cos(GantryAngle),0,ThePlan->GetSpotX(i,j)*sin(-GantryAngle)}; //Fixed X coordinate at isoplane
			DiffX=G4ThreeVector{(IsoX.x()-SourceX.x())/SADX,(IsoX.y()-SourceX.y())/SADX,(IsoX.z()-SourceX.z())/SADX};
			BaseSpot[i][j].Position.setX(SourceX.x()-DiffX.x()*(NozzleExit-SADX)); //X-position at nozzle exit
			BaseSpot[i][j].Position.setZ(SourceX.z()-DiffX.z()*(NozzleExit-SADX)); //Z-position at nozzle exit
			Y=ThePlan->GetSpotY(i,j); //Fixed Y coordinate at isoplane
			KY=Y/SADY;
			BaseSpot[i][j].Position.setY(KY*(SADY-NozzleExit)); //Y-position at nozzle exit

			BaseSpot[i][j].Direction.setX(DiffX.x()); //X-direction
			BaseSpot[i][j].Direction.setY(KY); //Y-direction
			BaseSpot[i][j].Direction.setZ(DiffX.z()); //Z-direction		
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::AddDose(G4int x, G4int y, G4int z,G4double D)
{
	using std::pow;
	DoseSpectrum[x][y][z][0]+=D;	//Dose
	DoseSpectrum[x][y][z][1]+=pow(D,2); //Dose squared
	++DoseSpectrum[x][y][z][2];	//Number of depositions
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
	delete ThePlan;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
	G4int NoSimPro=run->GetNumberOfEventToBeProcessed(); //need access to current G4Run to get NoEvent
	NormFactor=ThePlan->GetNoProtons()/NoSimPro;
	G4RunManager* RunManager=G4RunManager::GetRunManager();
	const DicomDetectorConstruction* Detector = static_cast<const DicomDetectorConstruction*> (RunManager->GetUserDetectorConstruction()); //get DetectorConstruction; need to recast since RunManager will return G4VUserDetectorConstruction pointer only		
	size_t XVoxNum=Detector->GetNoVoxelsX();
	size_t YVoxNum=Detector->GetNoVoxelsY();
	size_t ZVoxNum=Detector->GetNoVoxelsZ();
	DoseSpectrum = new G4double***[XVoxNum];
	for (size_t i=0;i<XVoxNum;++i)
	{
		DoseSpectrum[i]= new G4double**[YVoxNum];
		for (size_t j=0;j<YVoxNum;++j)
		{
			DoseSpectrum[i][j]=new G4double*[ZVoxNum];
			for (size_t k=0;k<ZVoxNum;++k) 
			{
				DoseSpectrum[i][j][k]=new G4double[BINS]; //TotalEnergyDep,TotalEnergyDep^2,NoEnergyDeps	
				for (size_t l=0;l<BINS;++l) DoseSpectrum[i][j][k][l]=0.;
			}
		}
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run*)
{
	G4cout << "Theta: " << GetBaseTheta(0,0) << " Phi: " << GetBasePhi(0,0) << " Psi: " << GetBasePsi(0,0) << G4endl;
	G4cout << "X: " << GetBaseX(0,0) << " Y: " << GetBaseY(0,0) << " Z: " << GetBaseZ(0,0) << G4endl;
	using std::ofstream;
	using std::pow;
	using std::sqrt;
	G4RunManager* RunManager=G4RunManager::GetRunManager();
	const DicomDetectorConstruction* Detector = static_cast<const DicomDetectorConstruction*> (RunManager->GetUserDetectorConstruction());
	size_t XVoxNum=Detector->GetNoVoxelsX();
	size_t YVoxNum=Detector->GetNoVoxelsY();
	size_t ZVoxNum=Detector->GetNoVoxelsZ();

	ofstream DoseMatrix{"../../OUTPUTDATA/DoseMatrix.dat",ofstream::binary};
	//ofstream DoseUncertainty{"/home/erik/PhDMain/Geant4/PROJECTS/DICOMSKANDION/OUTPUTDATA/DoseUncertainty.dat",ofstream::binary};
	
	DoseMatrix.write(reinterpret_cast<char*>(&XVoxNum),sizeof(size_t));
	DoseMatrix.write(reinterpret_cast<char*>(&YVoxNum),sizeof(size_t));
	DoseMatrix.write(reinterpret_cast<char*>(&ZVoxNum),sizeof(size_t));
	
	//DoseUncertainty.write(reinterpret_cast<char*>(&XVoxNum),sizeof(size_t));
	//DoseUncertainty.write(reinterpret_cast<char*>(&YVoxNum),sizeof(size_t));
	//DoseUncertainty.write(reinterpret_cast<char*>(&ZVoxNum),sizeof(size_t));
	G4double VoxelDose=-1, VoxelMass=0;
	for (size_t i=0;i<ZVoxNum;++i)
	{
		for (size_t j=0;j<YVoxNum;++j)
		{
			for (size_t k=0;k<XVoxNum;++k)
			{	
				VoxelDose=DoseSpectrum[k][j][i][0];
				//if (VoxelDose>1) G4cout << VoxelDose << G4endl;				
				//DoseSquare=DoseSpectrum[k][j][i][1];
				//N=DoseSpectrum[k][j][i][2];
				//S=(DoseSquare/N-pow(VoxelDose/N,2))*N/(N-1);	//sample population variance
				//S=sqrt(S)/sqrt(N);	//standard deviation
				//S=N*S*MeV2J/VoxelMass;
				//DoseUncertainty.write(reinterpret_cast<char*>(&S),sizeof(G4double));
				
	
				VoxelMass=Detector->GetVoxelMass(k,j,i);
				VoxelDose*=MeV2J/VoxelMass*NormFactor;
				DoseMatrix.write(reinterpret_cast<char*>(&VoxelDose),sizeof(G4double));
			}
		}
	}
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
