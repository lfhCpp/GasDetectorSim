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
// $Id: GarfieldPhysicsList.cc 999997 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldPhysicsList.cc
/// \brief Implementation of the GarfieldPhysicsList class

#include "GarfieldPhysicsList.hh"

#include "G4EmConfigurator.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4LossTableManager.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "GarfieldPhysics.hh"
#include "QGSP_BERT_HP.hh"
#include "QBBC.hh"


#include "G4EmParameters.hh"//添加的
#include "G4EmProcessOptions.hh"
#include "G4IonFluctuations.hh"
#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"


#include "G4PhysListFactory.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldPhysicsList::GarfieldPhysicsList() : G4VModularPhysicsList() {
  G4int verb = 0;
  SetVerboseLevel(verb);
  defaultCutValue = 1 * CLHEP::mm;

    // 设置电磁过程选项
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(10 * eV);         // 设置最小能量     default 100 eV
  emOptions.SetMaxEnergy(10 * TeV);        // 设置最大能量     default 100 TeV
  emOptions.SetDEDXBinning(12 * 10);       // 设置能量损失bin  default=12*7
  emOptions.SetLambdaBinning(12 * 10);     // 设置能量损失bin  default=12*7


  //QGSP_BIC_EMY：这是一个经过验证的物理列表，适用于低能重离子的模拟。它结合了 QGSP 和 BIC（Binary Cascade）模型，对电磁过程采用了标准的 Geant4 模型。
  //FTFP_BERT_EMV：这个物理列表结合了 Fritiof（FTF）模型和 Bertini（BERT）模型，并使用 Geant4 标准的电磁物理模型。它适用于高能重离子的模拟。
  //QGSP_BERT_HP默认列表 QGSP_BERT_HP* physicsList = new QGSP_BERT_HP;
  //QGSP_BERT_HP* physicsList = new QGSP_BERT_HP;


  //G4VModularPhysicsList* physicsList = new QBBC;

  // 创建物理列表工厂
   G4PhysListFactory factory;
  // 指定使用FTFP_BERT_EMV物理列表
   G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QGSP_BIC_EMY");

  for (G4int i = 0;; ++i) {
    G4VPhysicsConstructor* elem =
        const_cast<G4VPhysicsConstructor*>(physicsList->GetPhysics(i));
    if (elem == NULL) break;
    G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
    RegisterPhysics(elem);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldPhysicsList::~GarfieldPhysicsList() {}

void GarfieldPhysicsList::AddParameterisation() {
  GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();

  std::string ionizationModel = garfieldPhysics->GetIonizationModel();


  auto fastSimProcess_garfield = new G4FastSimulationManagerProcess("G4FSMP_garfield");

  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4EmConfigurator* config = G4LossTableManager::Instance()->EmConfigurator();
    G4LossTableManager::Instance()->SetVerbose(1);

    auto particleName = particle->GetParticleName();
    if (garfieldPhysics->FindParticleName(particleName, "garfield")) {
      pmanager->AddDiscreteProcess(fastSimProcess_garfield);
    }

    if (garfieldPhysics->FindParticleName(particleName, "geant4")) {
      double eMin = MeV * garfieldPhysics->GetMinEnergyMeVParticle(
          particleName, "geant4");
      double eMax = MeV * garfieldPhysics->GetMaxEnergyMeVParticle(
          particleName, "geant4");
      if (ionizationModel == "PAI") {
        G4PAIModel* pai = new G4PAIModel(particle, "G4PAIModel");
        if (particleName == "e-" || particleName == "e+") {
          config->SetExtraEmModel(particleName, "eIoni", pai,
                                  "RegionGarfield", eMin, eMax, pai);
        } else if (particleName == "mu-" || particleName == "mu+") {
          config->SetExtraEmModel(particleName, "muIoni", pai,
                                  "RegionGarfield", eMin, eMax, pai);
        } else if (particleName == "proton" ||
                   particleName == "pi+" || particleName == "pi-") {
          config->SetExtraEmModel(particleName, "hIoni", pai,
                                  "RegionGarfield", eMin, eMax, pai);
        } 
        else if (particleName == "alpha" || particleName == "He3" ||
                   particleName == "GenericIon") {
          config->SetExtraEmModel(particleName, "ionIoni", pai,
                                  "RegionGarfield", eMin, eMax, pai);
        }
        // else if (particleName == "alpha" || particleName == "He3" || particleName == "GenericIon") {
        //       G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
        //       G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
        //       G4double eth = 2. * MeV * particle->GetPDGMass() / CLHEP::proton_mass_c2;
        //       config->SetExtraEmModel(particleName, "ionIoni", mod1, "", 0.0, eth,
        //                                 new G4IonFluctuations());

        //       config->SetExtraEmModel(particleName, "ionIoni", mod2, "", eth, 100 * TeV,
        //                                 new G4UniversalFluctuation());
        //G4cout<<"--------------------------------------------------------eMin:"<<eMin<<" eth:"<<eth<<G4endl;
        // }
      } else if (ionizationModel == "PAIPhot") {
        G4PAIPhotModel* paiPhot = new G4PAIPhotModel(particle, "G4PAIModel");
        if (particleName == "e-" || particleName == "e+") {
          config->SetExtraEmModel(particleName, "eIoni", paiPhot,
                                  "RegionGarfield", eMin, eMax, paiPhot);
        } else if (particleName == "mu-" || particleName == "mu+") {
          config->SetExtraEmModel(particleName, "muIoni", paiPhot, 
                                  "RegionGarfield", eMin, eMax, paiPhot);
        } else if (particleName == "proton" ||
                   particleName == "pi+" || particleName == "pi-") {
          config->SetExtraEmModel(particleName, "hIoni", paiPhot,
                                  "RegionGarfield", eMin, eMax, paiPhot);
        } 
        else if (particleName == "alpha" || particleName == "He3" ||
                   particleName == "GenericIon") {
          config->SetExtraEmModel(particleName, "ionIoni", paiPhot, 
                                  "RegionGarfield", eMin, eMax, paiPhot);
        }
        // else if (particleName == "alpha" || particleName == "He3" || particleName == "GenericIon") {
        //       G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
        //       G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
        //       G4double eth = 2. * MeV * particle->GetPDGMass() / CLHEP::proton_mass_c2;
        //       config->SetExtraEmModel(particleName, "ionIoni", mod1, "", 0.0, eth,
        //                                 new G4IonFluctuations());
        //       config->SetExtraEmModel(particleName, "ionIoni", mod2, "", eth, 100 * TeV,
        //                                 new G4UniversalFluctuation());
        //G4cout<<"--------------------------------------------------------eMin:"<<eMin<<" eth:"<<eth<<G4endl;
        // }
      }
    }
  }
}

void GarfieldPhysicsList::SetCuts() {
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(21 * eV,//21  def：100eV
                                                                  100. * TeV);

  G4cout << "PhysicsList::SetCuts:";
  G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;

  SetCutsWithDefault();

  G4Region* region = G4RegionStore::GetInstance()->GetRegion("RegionGarfield");
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e-"));//1
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e+"));
  if (region) {
    region->SetProductionCuts(cuts);
  }

  //LEE-limit 默认100 eV
  G4EmParameters * emParams = G4EmParameters ::Instance();
  emParams->SetLowestElectronEnergy(30 * eV);//30

  DumpCutValuesTable();
}

void GarfieldPhysicsList::ConstructParticle() {
  G4VModularPhysicsList::ConstructParticle();
}

void GarfieldPhysicsList::ConstructProcess() {
  G4VModularPhysicsList::ConstructProcess();
  AddParameterisation();
}


