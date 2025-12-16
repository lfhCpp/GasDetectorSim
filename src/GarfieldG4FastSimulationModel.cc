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
// $Id: GarfieldG4FastSimulationModel.cc 999994 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldG4FastSimulationModel.cc
/// \brief Implementation of the GarfieldG4FastSimulationModel class

#include "GarfieldG4FastSimulationModel.hh"

#include <iostream>

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GDMLParser.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

#include "GarfieldPhysics.hh"

GarfieldG4FastSimulationModel::GarfieldG4FastSimulationModel(G4String modelName,
                                                             G4Region* envelope)
    : G4VFastSimulationModel(modelName, envelope) { //初始化列表：基类构造函数
  fGarfieldPhysics = GarfieldPhysics::GetInstance();//实例化：new GarfieldPhysics();
  fGarfieldPhysics->InitializePhysics();            //初始化：Garfield：ediumMagboltz、ComponentAnalyticField、Sensor
}

GarfieldG4FastSimulationModel::GarfieldG4FastSimulationModel(G4String modelName)
    : G4VFastSimulationModel(modelName) {
  fGarfieldPhysics = GarfieldPhysics::GetInstance();
  fGarfieldPhysics->InitializePhysics();
}

GarfieldG4FastSimulationModel::~GarfieldG4FastSimulationModel() {}

void GarfieldG4FastSimulationModel::WriteGeometryToGDML(
    G4VPhysicalVolume* physicalVolume) {
  G4GDMLParser* parser = new G4GDMLParser();
  remove("garfieldGeometry.gdml");
  parser->Write("garfieldGeometry.gdml", physicalVolume, false);
  delete parser;
}

G4bool GarfieldG4FastSimulationModel::IsApplicable(//判断粒子在Grafield中是否可用->系统调用
    const G4ParticleDefinition& particleType) {
  G4String particleName = particleType.GetParticleName();
  if (fGarfieldPhysics->FindParticleName(particleName, "garfield")) {
    return true;
  }
  return false;
}

G4bool GarfieldG4FastSimulationModel::ModelTrigger(//模式触发，确定模拟程序->系统调用
    const G4FastTrack& fastTrack) {
  double ekin_MeV = fastTrack.GetPrimaryTrack()->GetKineticEnergy() / MeV;
  G4String particleName =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
  if (fGarfieldPhysics->FindParticleNameEnergy(particleName, ekin_MeV,
                                               "garfield")) {
    return true;
  }
  return false;
}

void GarfieldG4FastSimulationModel::DoIt(const G4FastTrack& fastTrack,//Garfield模拟
                                         G4FastStep& fastStep) {

  G4ThreeVector localdir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localpos = fastTrack.GetPrimaryTrackLocalPosition();

  double ekin_MeV = fastTrack.GetPrimaryTrack()->GetKineticEnergy() / MeV;
  double globalTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();

  G4String particleName =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);

  if (particleName == "kaon+") {
    particleName = "K+";
  } else if (particleName == "kaon-") {
    particleName = "K-";
  } else if (particleName == "anti_proton") {
    particleName = "anti-proton";
  }

  static int ProgressBar=0;  //进度条
  if(fGarfieldPhysics->CoutMessageToDebug()){
     G4cout << " ===================================================== "<<G4endl;
     G4cout<<"BeforeInDoIt:Ek x y z "<<ekin_MeV<<" "<<localpos.x()+87./2
          <<" "<< localpos.y()+147/2.<<" "<< localpos.z()+69.9/2<<G4endl;
  }
  else{
  if(ProgressBar%2000==0){
     G4cout<<ProgressBar<<" DeltaElectrons has been simulated."<<G4endl;
    }ProgressBar++;
  }

  fGarfieldPhysics->DoIt(
      particleName, ekin_MeV, globalTime, (localpos.x()+87./2) / CLHEP::cm,
      (localpos.y()+147/2.) / CLHEP::cm, (localpos.z()+69.9/2) / CLHEP::cm,
      localdir.x(), localdir.y(), localdir.z());

  fastStep.SetTotalEnergyDeposited(fGarfieldPhysics->GetEnergyDeposit_MeV());

  //是否在G4中创建次级粒子
  fGarfieldPhysics->EnableCreateSecondariesInGeant4(false);

  if (!fGarfieldPhysics->GetCreateSecondariesInGeant4()) return;
  std::vector<GarfieldParticle*>* secondaryParticles =
      fGarfieldPhysics->GetSecondaryParticles();
     G4cout<<"------------->fGarfieldPhysics->EnableCreateSecondariesInGeant4(true);"<< G4endl;

  if (secondaryParticles->empty()) return;
  fastStep.SetNumberOfSecondaryTracks(secondaryParticles->size());
  G4cout<<"------------->secondaryParticles->Not empty"<< G4endl;

  G4double totalEnergySecondaries_MeV = 0;

  for (auto it = secondaryParticles->begin(); it != secondaryParticles->end(); ++it) {
    G4double eKin_MeV = (*it)->getEkin_MeV();
    G4double time = (*it)->getTime();
    G4ThreeVector momentumDirection((*it)->getDX(), (*it)->getDY(), 
                                    (*it)->getDZ());
    G4ThreeVector position((*it)->getX_mm(), (*it)->getY_mm(),
                           (*it)->getZ_mm());

    if ((*it)->getParticleName() == "e-") {

        G4cout<<"------------->e-  X:"<<(*it)->getX_mm()
                                <<"Y:"<<(*it)->getY_mm()
                                <<"Z:"<<(*it)->getZ_mm()
                                <<G4endl;

      G4DynamicParticle particle(G4Electron::ElectronDefinition(),
                                 momentumDirection, eKin_MeV);
     // G4DynamicParticle particle(G4Proton::ProtonDefinition(),
     //                            momentumDirection,eKin_MeV*10000);
      fastStep.CreateSecondaryTrack(particle, position, time, true);
      totalEnergySecondaries_MeV += eKin_MeV;
    } else if ((*it)->getParticleName() == "gamma") {
       G4DynamicParticle particle(G4Gamma::GammaDefinition(),
                                 momentumDirection, eKin_MeV);
      fastStep.CreateSecondaryTrack(particle, position, time, true);
      totalEnergySecondaries_MeV += eKin_MeV;
    }
  }
  
}
