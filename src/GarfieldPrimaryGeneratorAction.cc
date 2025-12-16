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
// $Id: GarfieldPrimaryGeneratorAction.cc 999998 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldPrimaryGeneratorAction.cc
/// \brief Implementation of the GarfieldPrimaryGeneratorAction class

#include "GarfieldPrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleMomentum.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ChargedGeantino.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldPrimaryGeneratorAction::GarfieldPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

//  G4ParticleDefinition* particleDefinition
//    = G4ParticleTable::GetParticleTable()->FindParticle("proton");
//  fParticleGun->SetParticleDefinition(particleDefinition);
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
//  fParticleGun->SetParticleEnergy(300*56.*MeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldPrimaryGeneratorAction::~GarfieldPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // This function is called at the begining of event


        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle = particleTable->FindParticle("alpha");
        if(particle)
            fParticleGun->SetParticleDefinition(particle);
        else
            G4cout<<"##Null pp in PrimaryGeneratorAction::SetParticleGun()"<<G4endl;

        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
        fParticleGun->SetParticlePosition(G4ThreeVector(87./2.,147/2.,69.9)); //初始位置  注意还有阴极板厚度1mm
        fParticleGun->SetParticleEnergy(5.153*MeV);// Pu-239 5.153    PU-238 5.486MeV


//    //声明离子
//    G4ParticleDefinition *particle=G4IonTable::GetIonTable()->GetIon(26,56,0.0); //Fe, z,a,ex
//    fParticleGun->SetParticleDefinition(particle);
//    fParticleGun->SetParticleCharge(0.);                          //电荷态，EM 物理过程自动修正
//    fParticleGun->SetParticleEnergy(3000.*MeV);                 //动能                            //时间
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));//动量方向


//    G4double worldZHalfLength = 0;
//    G4LogicalVolume* worlLV
//      = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
//    G4Box* worldBox = 0;
//    if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
//    if ( worldBox ) {
//      worldZHalfLength = worldBox->GetZHalfLength();
//    }
//    else  {
//      G4ExceptionDescription msg;
//      msg << "World volume of box not found." << G4endl;
//      msg << "Perhaps you have changed geometry." << G4endl;
//      msg << "The gun will be place in the center.";
//      G4Exception("GarfieldPrimaryGeneratorAction::GeneratePrimaries()",
//        "MyCode0002", JustWarning, msg);
//    }
//    // Set gun position
//    fParticleGun
//      ->SetParticlePosition(G4ThreeVector(50., 50., 0));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
