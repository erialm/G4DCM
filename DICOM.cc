//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DICOM.cc 101109 2016-11-07 08:14:53Z gcosmo $
//
/// \file medical/DICOM/DICOM.cc
/// \brief Main program of the medical/DICOM example
//
// The code was written by :
//        *Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
// *******************************************************


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "G4Threading.hh"
#else
#include "G4RunManager.hh"
#endif

#include "globals.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4GenericPhysicsList.hh"

#include "DicomRegularDetectorConstruction.hh"
#include "DicomNestedParamDetectorConstruction.hh"
#include "DicomPartialDetectorConstruction.hh"

#include "ActionInitialization.hh"

#ifdef G4_DCMTK
#include "DicomFileMgr.hh"
#else
#include "DicomHandler.hh"
#endif
   
#include "DicomIntersectVolume.hh"
#include "PhysicsList.hh"
#include "G4tgrMessenger.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "Shielding.hh"

//=================================================================================

int main(int argc,char** argv)
{
  
  new G4tgrMessenger;
  char* part = getenv( "DICOM_PARTIAL_PARAM" );
  G4bool bPartial = FALSE;
  if( part && G4String(part) == "1" ) {
    bPartial = TRUE;
  }
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
#ifdef G4MULTITHREADED
  
  G4MTRunManager* runManager = new G4MTRunManager;
  G4int nthreads=G4Threading::G4GetNumberOfCores();
  runManager->SetNumberOfThreads(nthreads);
  
  G4cout << "\n\n\tDICOM running in multithreaded mode with " << nthreads 
         << " threads\n\n" << G4endl;
  
  
#else
  G4RunManager* runManager = new G4RunManager;
  G4cout << "\n\n\tDICOM running in serial mode\n\n" << G4endl;
  
#endif
  
  DicomDetectorConstruction* theGeometry = nullptr;
  
#ifdef G4_DCMTK
  DicomFileMgr* theFileMgr = nullptr;
#else
  DicomHandler* dcmHandler = nullptr;
#endif
  
  if( !bPartial ){
#ifdef G4_DCMTK
    
    theFileMgr = DicomFileMgr::GetInstance();
    theFileMgr->Convert("Data.dat");
    
#else
    // Treatment of DICOM images before creating the G4runManager
    dcmHandler = new DicomHandler;
    dcmHandler->CheckFileFormat();
#endif
    
    // Initialisation of physics, geometry, primary particles ...
     char* nest = getenv( "DICOM_NESTED_PARAM" );
     if( nest && G4String(nest) == "1" ) {
      theGeometry = new DicomNestedParamDetectorConstruction();
	G4cout << "Doing Nested!" << G4endl;
    } else {
      theGeometry = new DicomRegularDetectorConstruction();
	G4cout << "Doing Regular!" << G4endl;
    }
  } else {
    theGeometry = new DicomPartialDetectorConstruction();
	G4cout << "Doing Partial!" << G4endl;
  }    
  runManager->SetUserInitialization(theGeometry);
  
 PhysicsList* phys = new PhysicsList(); 
 runManager->SetUserInitialization(phys);
  
  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization("../../INPUTDATA/Plan.txt"));
   
  new DicomIntersectVolume();
  
#ifdef G4VIS_USE
  // visualisation manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  
  if (argc==1)
    {
      runManager->Initialize();
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  else
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      G4int NoProtons=runManager->GetNumberOfEventsToBeProcessed();
      G4cout << "Running a total number of : " << NoProtons << " protons" << G4endl;
      runManager->Initialize();
      UImanager->ApplyCommand(command+fileName);
      runManager->BeamOn(NoProtons);
    }
  
  delete runManager;
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  
  if( !bPartial ) {
#ifdef G4_DCMTK
    delete theFileMgr;
#else
    delete dcmHandler;
#endif
  }
  
  return 0;
}

