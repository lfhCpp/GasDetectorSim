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
/// \file GarfieldPhysics.cc
/// \brief Implementation of the GarfieldPhysics class
#include "GarfieldPhysics.hh"
#include "GarfieldAnalysis.hh"

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/DriftLineRKF.hh"

/////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"

#include <TApplication.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/MediumConductor.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"
#include "Garfield/TrackSrim.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/Plotting.hh"
#include <TSystem.h>
#include <cstdlib>
#include <iostream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLatex.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"


using namespace Garfield;
using namespace std;

/////////////////////////////////////////


void GarfieldPhysics::ModifyLineData(std::string fileName, int lineNum, char* lineData){//修正读入的COMSOL文件

        cout<<"ModifyLineData"<<endl;
        ifstream in;
        in.open(fileName.c_str());
        std::string strFileData = "";
        int line = 1;
        char tmpLineData[1024] = {0};
        while(in.getline(tmpLineData, sizeof(tmpLineData)))
        {
                if (line == lineNum)
                {
                        strFileData += string(lineData);
                        strFileData += "\n";
                }
                else
                {
                        strFileData += string(tmpLineData);
                        strFileData += "\n";
                }
                line++;
        }
        in.close();
        //写入文件
        ofstream out;
        out.open(fileName.c_str());
        out.flush();
        out<<strFileData;
        out.close();
}


GarfieldPhysics* GarfieldPhysics::fGarfieldPhysics = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GarfieldPhysics* GarfieldPhysics::GetInstance() {
  if (!fGarfieldPhysics) {
    fGarfieldPhysics = new GarfieldPhysics();
  }
  return fGarfieldPhysics;
}

void GarfieldPhysics::Dispose() {//程序结束执行
  delete fGarfieldPhysics;
  fGarfieldPhysics = 0;
}

GarfieldPhysics::GarfieldPhysics() {
  fSecondaryParticles = new std::vector<GarfieldParticle*>();
  fMediumMagboltz = 0;
  fSensor = 0;
  fComponentComsol = 0;
  fTrackHeed = 0;
  fGarfieldEDep = 0;//统计沉积能量
  NumberOfElectron=0;//统计电子数
  createSecondariesInGeant4 = false;
  CoutMessage=false;
  Draw=false;
  fIonizationModel = "PAIPhot";

  GarfieldPhysics::SetIonizationModel("PAI",true);
  cDrift = new TCanvas("cDrift", "", 600, 600);
  Signal_pCathode  = new TCanvas("Cathode Signal", "", 600, 600);
  Signal_pGrid = new TCanvas("Grid Signal", "", 600, 600);
  Signal_pAnode = new TCanvas("Anode Signal", "", 600, 600);
  cCell=new TCanvas("cCell", "", 600, 600);
  CellView=new Garfield::ViewCell();
  fieldViewXY = new  Garfield::ViewField();
  fieldViewXZ = new  Garfield::ViewField();

}

GarfieldPhysics::~GarfieldPhysics() {
  DeleteSecondaryParticles();
  delete fSecondaryParticles;
  delete fMediumMagboltz;
  delete fSensor;
  delete fComponentComsol;
  delete fTrackHeed;
  delete MCdrift;
  delete RKFdrift;
  delete avalanche;

  delete cDrift;
  delete Signal_pCathode;
  delete Signal_pGrid;
  delete Signal_pAnode;
  delete CellView;
  delete fieldViewXY;
  delete fieldViewXZ;

  std::cout << "Deconstructor GarfieldPhysics" << std::endl;
  G4cout<<"More information to debug -> #define CoutMessageToDebug in GarfieldPhysics.cc"<<G4endl;
}

std::string GarfieldPhysics::GetIonizationModel() { return fIonizationModel; }

void GarfieldPhysics::SetIonizationModel(std::string model, bool useDefaults) {//信使
  if (model != "PAIPhot" && model != "PAI" && model != "Heed") {
    std::cout << "Unknown ionization model " << model << std::endl;
    std::cout << "Using PAIPhot as default model!" << std::endl;
    model = "PAIPhot";
  }
  fIonizationModel = model;

  fIonizationModel="PAIPhot";//lfh

  if (fIonizationModel == "PAIPhot" || fIonizationModel == "PAI") {
    if (useDefaults == true) {
      // Particle types and energies for which the G4FastSimulationModel with
      // Garfield++ is valid
      this->AddParticleName("e-", 1e-6, 1., "garfield");//默认：this->AddParticleName("e-", 1e-6, 1e-3, "garfield");
      this->AddParticleName("gamma", 1e-6, 1e+8, "garfield");

      // Particle types and energies for which the PAI or PAIPhot model is valid
      this->AddParticleName("e-", 0, 1e+8, "geant4");
      this->AddParticleName("e+", 0, 1e+8, "geant4");
      this->AddParticleName("mu-", 0, 1e+8, "geant4");
      this->AddParticleName("mu+", 0, 1e+8, "geant4");
      this->AddParticleName("proton", 0, 1e+8, "geant4");
      this->AddParticleName("pi+", 0, 1e+8, "geant4");
      this->AddParticleName("pi-", 0, 1e+8, "geant4");
      this->AddParticleName("alpha", 0, 1e+8, "geant4");
      this->AddParticleName("He3", 0, 1e+8, "geant4");
      this->AddParticleName("GenericIon", 0, 1e+8, "geant4");
    }

  } else if (fIonizationModel == "Heed") {
    if (useDefaults == true) {
      // Particle types and energies for which the G4FastSimulationModel with
      // Garfield++ is valid
      this->AddParticleName("gamma", 1e-6, 1e+8, "garfield");
      this->AddParticleName("e-", 6e-2, 1e+7, "garfield");
      this->AddParticleName("e+", 6e-2, 1e+7, "garfield");
      this->AddParticleName("mu-", 1e+1, 1e+8, "garfield");
      this->AddParticleName("mu+", 1e+1, 1e+8, "garfield");
      this->AddParticleName("pi-", 2e+1, 1e+8, "garfield");
      this->AddParticleName("pi+", 2e+1, 1e+8, "garfield");
      this->AddParticleName("kaon-", 1e+1, 1e+8, "garfield");
      this->AddParticleName("kaon+", 1e+1, 1e+8, "garfield");
      this->AddParticleName("proton", 9.e+1, 1e+8, "garfield");
      this->AddParticleName("anti_proton", 9.e+1, 1e+8, "garfield");
      this->AddParticleName("deuteron", 2.e+2, 1e+8, "garfield");
      this->AddParticleName("alpha", 4.e+2, 1e+8, "garfield");
    }
  }
}

void GarfieldPhysics::AddParticleName(const std::string particleName,
                                      double ekin_min_MeV, double ekin_max_MeV,
                                      std::string program) {
  if (ekin_min_MeV >= ekin_max_MeV) {
    std::cout << "Ekin_min=" << ekin_min_MeV
              << " keV is larger than Ekin_max=" << ekin_max_MeV << " keV"
              << std::endl;
    return;
  }

  if (program == "garfield") {
    std::cout << "Garfield model (Heed) is applicable for G4Particle "
              << particleName << " between " << ekin_min_MeV << " MeV and "
              << ekin_max_MeV << " MeV" << std::endl;

    fMapParticlesEnergyGarfield.insert(std::make_pair(
        particleName, std::make_pair(ekin_min_MeV, ekin_max_MeV)));
  } else {
    std::cout << fIonizationModel << " is applicable for G4Particle "
              << particleName << " between " << ekin_min_MeV << " MeV and "
              << ekin_max_MeV << " MeV" << std::endl;
    fMapParticlesEnergyGeant4.insert(std::make_pair(
        particleName, std::make_pair(ekin_min_MeV, ekin_max_MeV)));
  }
}

bool GarfieldPhysics::FindParticleName(std::string name, std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) return true;
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) return true;
  }
  return false;
}

bool GarfieldPhysics::FindParticleNameEnergy(std::string name, double ekin_MeV,
                                             std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      if (range.first <= ekin_MeV && range.second >= ekin_MeV) {
        return true;
      }
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      if (range.first <= ekin_MeV && range.second >= ekin_MeV) {
        return true;
      }
    }
  }
  return false;
}

double GarfieldPhysics::GetMinEnergyMeVParticle(std::string name,
                                                std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      return range.first;
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      return range.first;
    }
  }
  return -1;
}

double GarfieldPhysics::GetMaxEnergyMeVParticle(std::string name,
                                                std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      return range.second;
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      return range.second;
    }
  }
  return -1;
}

void GarfieldPhysics::InitializePhysics() {
#define ar_ch4               //使用ar_ch4气体，注意应同步更改GarfieldDetectorConstruction.cc的材料
#define MC                 //使用Garfield::AvalancheMC模拟
//#define plot               //用于调试，画出电场、电势、读出信号图
//#define CoutMessageToDebug //输出更多的模拟过程信息用于调试
#ifdef  ar_ch4
    //Define the gas mixture.
    fMediumMagboltz = new Garfield::MediumMagboltz();
    fMediumMagboltz->SetComposition("ar", 90.,"CH4", 10.);
    fMediumMagboltz->SetTemperature(293.15);
    fMediumMagboltz->SetPressure(760.);
    fMediumMagboltz->Initialise(true);
    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    fMediumMagboltz->LoadGasFile("./GasFile/ar_90_ch4_10_760Torr.gas");         //加载气体表
#else
    // Define the gas mixture.
    fMediumMagboltz = new Garfield::MediumMagboltz();
    fMediumMagboltz->SetComposition("ar", 70., "co2", 30.);
    fMediumMagboltz->SetTemperature(293.15);
    fMediumMagboltz->SetPressure(760.);
    fMediumMagboltz->Initialise(true);
    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    fMediumMagboltz->LoadGasFile("ar_70_co2_30_1000mbar.gas");
#endif

    const std::string path = std::getenv("GARFIELD_INSTALL");
    fMediumMagboltz->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    G4cout<<"\tW------->Medium:"<<fMediumMagboltz->GetW()<<G4endl;


    //
    // ---- Load the field map ----
    //

    std::string file_pAnode1, file_pCathode,file_mesh, file_field;



    file_pAnode1 =   "./ComsolFile/pAnode1.txt";
    file_pCathode =   "./ComsolFile/pCathode.txt";
    file_field=   "./ComsolFile/wfield.txt";
    file_mesh=   "./ComsolFile/mesh.mphtxt";

    this->GarfieldPhysics::ModifyLineData(file_pAnode1.c_str(),
    9,
    "% x                       y                        z                        V (V)");



    this->GarfieldPhysics::ModifyLineData(file_pCathode.c_str(),
    9,
    "% x                       y                        z                        V (V)");

    this->GarfieldPhysics::ModifyLineData(file_field.c_str(),
    9,
    "% x                       y                        z                        V (V)");

    //
    // ---- Load the field map ----

    fComponentComsol = new Garfield::ComponentComsol();
    fComponentComsol->Initialise(file_mesh.c_str(), "./GasFile/dielectrics.dat", file_field.c_str(), "mm");
    fComponentComsol->SetWeightingPotential(file_pAnode1.c_str(),label_pAnode1);
    fComponentComsol->SetWeightingPotential(file_pCathode.c_str(),label_pCathode);
    fComponentComsol->EnableMirrorPeriodicityX();
    fComponentComsol->EnableMirrorPeriodicityY();
    fComponentComsol->PrintRange();

    // Associate the gas with the corresponding field map material.
    const unsigned int nMaterials = fComponentComsol->GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i) {
            const double eps = fComponentComsol->GetPermittivity(i);
            if (eps == 1.) fComponentComsol->SetMedium(i, fMediumMagboltz);
    }
    fComponentComsol->PrintMaterials();



    //
    // ---- Signal时间窗 ----
    //
    constexpr double tMin = 0.;
    constexpr double tMax = 3000.;   // ns
    constexpr double tStep = 1;
    constexpr int nTimeBins = int((tMax-tMin)/tStep);


  // Create the sensor.
  fSensor = new Garfield::Sensor();
  fSensor->Clear();
  fSensor->AddComponent(fComponentComsol);
  fSensor->SetArea(0., 0., 0.,20., 20.,73);
  fSensor->AddElectrode(fComponentComsol, label_pCathode);
  fSensor->AddElectrode(fComponentComsol,label_pAnode1);
  fSensor->SetTimeWindow(0., tStep, nTimeBins);

  fTrackHeed = new Garfield::TrackHeed();
  fTrackHeed->SetSensor(fSensor);
  fTrackHeed->EnableDeltaElectronTransport();
  //fTrackHeed->EnableDebugging();

  //Garfield::AvalancheMC MCdrift; in Garfield.hh
  MCdrift= new AvalancheMC();
  MCdrift->SetSensor(fSensor);
  //MCdrift->SetDistanceSteps(1.e-4);//过小会限制模拟速度，视模拟精度，建议使用默认
  MCdrift->EnableAttachmentMap();
  MCdrift->EnableTownsendMap();
  MCdrift->EnableSignalCalculation(true);
  MCdrift->UseWeightingPotential(true);
  MCdrift->SetTimeWindow(tMin, tMax);

  //Garfield::DriftLineRKF RKFdrift; in Garfield.hh
  RKFdrift=new DriftLineRKF();
  RKFdrift->SetSensor(fSensor);
  RKFdrift->EnableSignalCalculation(true);

  //Garfield::AvalancheMicroscopic avalanche; in Garfield.hh
  avalanche=new AvalancheMicroscopic();
  avalanche->SetSensor(fSensor);
  avalanche->EnableSignalCalculation(true);



  //    // ---- 几何、track、电子轨迹 可视化 ----


  #ifdef  plot
  Draw=true;
            driftView.SetClusterMarkerSize(0.1);
            driftView.SetCollisionMarkerSize(0.1);
            driftView.SetColourElectrons(kRed);
            driftView.SetColourTracks(kGreen + 3);
            driftView.SetPlaneXZ();
            MCdrift->EnablePlotting(&driftView);
            fTrackHeed->EnablePlotting(&driftView);
            //driftView.SetPlane(0, 0, 0, 0, 0, 0);
            //driftView.SetArea(-Dim_x,-Dim_y,Dim_z_down,
            //                   Dim_x, Dim_y,Dim_z_up);
            //driftView.SetArea();
            driftView.SetCanvas(cDrift);  

            Signalview_pCathode.SetSensor(fSensor);
            // fSensor->IntegrateSignal(label_pAnode1);
            Signalview_pCathode.SetCanvas(Signal_pCathode);
            //Signalview_pCathode.SetLabelY("Induced Charge [fC]");

            Signalview_pGrid.SetSensor(fSensor);
            // fSensor->IntegrateSignal(label_pGrid);
            Signalview_pGrid.SetCanvas(Signal_pGrid);
            //Signalview_pGrid.SetLabelY("Induced Charge [fC]");

            Signalview_pAnode.SetSensor(fSensor);
            //fSensor->IntegrateSignal(label_pCathode);
            Signalview_pAnode.SetCanvas(Signal_pAnode);
            //Signalview_pAnode.SetLabelY("Induced Charge [fC]");


            cFieldXZ = new TCanvas("cFieldXZ", "", 800, 800);
            fieldViewXZ->SetComponent(fComponentComsol);
            //fieldViewXZ->SetPlane(0, -0.1, 0, 0, 0.05, 0);
            //fieldViewXZ->SetArea(-0.1, -1, 0.3,0.1, 1,0.7);
            fieldViewXZ->SetPlaneXZ();
            fieldViewXZ->SetSensor(fSensor);
            fieldViewXZ->SetNumberOfSamples2d(500,500);
            //fieldViewXZ->SetElectricFieldRange(0,2000);
            fieldViewXZ->SetCanvas(cFieldXZ);
            fieldViewXZ->Plot("e","colz");

            cFieldXY = new TCanvas("cFieldXY", "", 900, 800);
            fieldViewXY->SetComponent(fComponentComsol);
            //fieldViewXY->SetPlane(0, 0, -0.1,  0, 0, 0.05);
            //fieldViewXY->SetArea(-0.1, -10, 0,0.1, 10,7.3);
            fieldViewXY->SetPlaneXZ();
            fieldViewXY->SetSensor(fSensor);
            fieldViewXY->SetNumberOfSamples2d(500,500);
            //fieldViewXY->SetElectricFieldRange(0,2000);
            fieldViewXY->SetCanvas(cFieldXY);
            fieldViewXY->PlotWeightingField(label_pAnode1,"e","colz");
  #endif

}//InitializePhysics() end here!


void GarfieldPhysics::DoIt(std::string particleName, double ekin_MeV,
                           double time, double x_cm, double y_cm, double z_cm,
                           double dx, double dy, double dz) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  fEnergyDeposit = 0;
  DeleteSecondaryParticles();




  // 簇信息
  double eKin_eV = ekin_MeV * 1e+6;
  double xc = 0., yc = 0., zc = 0., tc = 0., ee1 = 0.;
  // Number of electrons produced in a collision
  int nc = 0;
  // Energy loss in a collision
  double ec = 0.;
  // Dummy variable (not used at present)
  double extra = 0.;
  if (fIonizationModel != "Heed" || particleName == "gamma") {
    if (particleName == "gamma") {//输运光子与Delta电子
      fTrackHeed->TransportPhoton(x_cm, y_cm, z_cm, time, eKin_eV, dx, dy, dz,
                                  nc);
    } else {
      fTrackHeed->TransportDeltaElectron(x_cm, y_cm, z_cm, time, eKin_eV, dx,
                                         dy, dz, nc);
      fEnergyDeposit=eKin_eV;
      fGarfieldEDep =fGarfieldEDep + eKin_eV;
    }
    

   //统计电子数
   NumberOfElectron+=nc;

#ifdef CoutMessageToDebug
    CoutMessage=true;
    G4cout<<"InDoIt:ek,x_cm,y_cm,z_cm "<<ekin_MeV * 1e+6
         <<" "<<x_cm<<" "<<" "<< y_cm<<" "<< z_cm<<G4endl;
    G4cout<<"InDoIt:dx,dy,dz "
         <<" "<<dx<<" "<<" "<< dy<<" "<< dz<<G4endl;
    G4cout << " Number of Electron: " << nc <<G4endl; //当前cluster的电子个数
#endif

    for (int cl = 0; cl < nc; cl++) {//遍历团簇中的电子

      double xe, ye, ze, te;
      double ee, dxe, dye, dze;
      fTrackHeed->GetElectron(cl, xe, ye, ze, te, ee, dxe, dye, dze);

#ifdef CoutMessageToDebug
      G4cout << "Electron Start Position of Clusters & te " << xe <<","<< ye <<","<< ze <<","<<te<<G4endl;
      G4cout << "beginning to calculate the electron: "<< G4endl;
#endif
        nsum++;
        if (particleName == "gamma") {
          fEnergyDeposit += fTrackHeed->GetW();
          fGarfieldEDep+= fTrackHeed->GetW();
          //G4cout<<"particleName ==gamma"<<G4endl;
        }
       analysisManager->FillH3(1, ze * 10, xe * 10, ye * 10);
        if (createSecondariesInGeant4) {
          double newTime = te;
          if (newTime < time) {
            newTime += time;
          }
          fSecondaryParticles->push_back(new GarfieldParticle(
              "e-", ee, newTime, xe, ye, ze, dxe, dye, dze));
        }


#ifdef MC

        MCdrift->DriftElectron(xe, ye, ze, te);//开始计算电子
        // electrons start points and stop point
        double xe1, ye1, ze1, te1;
        double xe2, ye2, ze2, te2;
        int status;
        MCdrift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);

     #ifdef CoutMessageToDebug
        G4cout << " After drift Position of electron: " << xe1 <<","<< ye1 <<","<< ze1 <<G4endl;
        G4cout << " drift stop Position of electron: " << xe2 <<","<< ye2 <<","<< ze2 <<G4endl;
       // G4cout << " fEnergyDeposit: " << fEnergyDeposit <<G4endl;
     #endif

#else
        RKFdrift->DriftElectron(xe, ye, ze, te);//开始计算电子
        double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0.;
        int status = 0;
        RKFdrift->GetEndPoint(xe1, ye1, ze1, te1, status);
        ye1-=1.e-10;
        RKFdrift->DriftElectron(xe1, ye1, ze1, te1);
     #ifdef CoutMessageToDebug
        G4cout << "RKFdrift stop Position of electron:      " << xe1 <<","<< ye1 <<","<< ze1 <<","<< te1 <<G4endl;
     #endif
#endif

      }//遍历团簇中的电子*/
  } else {
    fTrackHeed->SetParticle(particleName);
    fTrackHeed->SetKineticEnergy(eKin_eV);
    fTrackHeed->NewTrack(x_cm, y_cm, z_cm, time, dx, dy, dz);
     G4cout << "SetKineticEnergy(eKin_eV) " <<eKin_eV <<G4endl;
    // loop clusters
    while (fTrackHeed->GetCluster(xc, yc, zc, tc, nc, ec, extra)) { //获取当前团簇信息
        nsum += nc;
        fEnergyDeposit += ec;
        fGarfieldEDep  += ec;
        G4cout << "e_ Number of Electron: " << nc <<G4endl; //当前cluster的电子个数
        for (int cl = 0; cl < nc; cl++) {//遍历团簇中的电子
          double xe, ye, ze, te;
          double ee, dxe, dye, dze;
          fTrackHeed->GetElectron(cl, xe, ye, ze, te, ee, dxe, dye, dze);
            analysisManager->FillH3(1, ze * 10, xe * 10, ye * 10);
            if (createSecondariesInGeant4) {
              double newTime = te;
              if (newTime < time) {
                newTime += time;
              }
              fSecondaryParticles->push_back(new GarfieldParticle(
                  "e-", ee, newTime, xe, ye, ze, dxe, dye, dze));
            }
            G4cout << "beginning to calculate the electron: e_"<< G4endl;
            MCdrift->DriftElectron(xe, ye, ze, te);

            // electrons start points and stop point
            double xe1, ye1, ze1, te1;
            double xe2, ye2, ze2, te2;

            int status;
            MCdrift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2,
                                      te2, status);

          }
        }
  }
  fGain = fAvalancheSize / nsum;
  constexpr bool plotDrift = false;// ---- 绘图控制 ----
  if (plotDrift){
    driftView.Plot2d(true);        //画漂移径迹，耗时长！
  }

    #ifdef  plot
         // Plot the induced current.
         Signalview_pCathode.PlotSignal(label_pCathode);
         Signalview_pGrid.PlotSignal(label_pGrid);
         Signalview_pAnode.PlotSignal(label_pAnode1);
    #endif

}//DoIt() end here!


std::vector<GarfieldParticle*>* GarfieldPhysics::GetSecondaryParticles() {
  return fSecondaryParticles;
}

void GarfieldPhysics::DeleteSecondaryParticles() {
  if (!fSecondaryParticles->empty()) {
    fSecondaryParticles->erase(fSecondaryParticles->begin(),
                               fSecondaryParticles->end());
  }
}
