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
// $Id: GarfieldSteppingAction.cc 999990 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldSteppingAction.cc
/// \brief Implementation of the GarfieldSteppingAction class

#include "GarfieldSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "GarfieldDetectorConstruction.hh"
#include "GarfieldEventAction.hh"

#include "GarfieldPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldSteppingAction::GarfieldSteppingAction(
    const GarfieldDetectorConstruction* detectorConstruction,
    GarfieldEventAction* eventAction)
    : G4UserSteppingAction(),
      fDetConstruction(detectorConstruction),
      fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldSteppingAction::~GarfieldSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldSteppingAction::UserSteppingAction(const G4Step* step) {
  // Collect energy and track length step by step

    G4StepPoint* preStepPoint = step->GetPreStepPoint();//获取数据接口
    G4String VolNamePre = preStepPoint->GetPhysicalVolume()->GetName();
    G4ThreeVector PosPre = preStepPoint->GetPosition();
    G4Track* aTrack = step->GetTrack();//获取数据接口
    G4ParticleDefinition* theparticle = aTrack->GetDefinition();
    G4String PName = theparticle->GetParticleName();
    G4double EkPre = preStepPoint->GetKineticEnergy();//

  // get volume of the current step
  G4VPhysicalVolume* volume =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();

    GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();
  //  if(garfieldPhysics->CoutMessageToDebug()){
  //       G4cout<<"VolNamePre: "<<VolNamePre<<" PNameG4: "<<PName<<
  //           " EkPre:"<<EkPre*1e6<<" "<<" energy deposit:"<<edep*1e6<<" xyz:"
  //         <<PosPre.x()<<" "<<PosPre.y()<<" "<<PosPre.z()<<G4endl;
  //       static G4double ill=0;
  //       if(edep!=garfieldPhysics->GetEnergyDeposit_MeV()){
  //           ill+=edep*1e6;
  //           G4cout<<"edep!="<<edep*1e6<<" dif:"<<(garfieldPhysics->GetEnergyDeposit_MeV())*1e6<<" "<<ill<<G4endl;
  //           G4cout<<"VolNamePre: "<<VolNamePre<<" PNameG4: "<<PName<<G4endl;
  //           G4cout<<" xyz:"<<PosPre.x()<<" "<<PosPre.y()<<" "<<PosPre.z()<<G4endl;
  //       }
  //  }

  //  static G4double ill=0;
  //       if(edep!=garfieldPhysics->GetEnergyDeposit_MeV()){
  //           ill+=edep*1e6;
  //           G4cout<<"edep!="<<edep*1e6<<" dif:"<<(garfieldPhysics->GetEnergyDeposit_MeV())*1e6<<" "<<ill<<G4endl;
  //           G4cout<<"VolNamePre: "<<VolNamePre<<" PNameG4: "<<PName<<G4endl;
  //           G4cout<<" xyz:"<<PosPre.x()<<" "<<PosPre.y()<<" "<<PosPre.z()<<G4endl;
  //         /// 获取相互作用过程信息
  //         const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  //         G4String interactionType = (process) ? process->GetProcessName() : "No Process";
  //         // 输出相互作用过程信息
  //         G4cout << "Interaction Process: " << interactionType << G4endl;
  //  }



  // step length
  G4double stepLength = 0.;
  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
    stepLength = step->GetStepLength();
  }

  if (volume == fDetConstruction->GetAbsorberPV()) {
    fEventAction->AddAbs(edep, stepLength);
  }

  if (volume == fDetConstruction->GetChamberPV()) {
    fEventAction->AddGas(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
