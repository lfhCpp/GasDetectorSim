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
/// \file GarfieldPhysics.hh
/// \brief Definition of the GarfieldPhysics class
/// \author D. Pfeiffer
//
#ifndef GarfieldPhysics_h
#define GarfieldPhysics_h 1

#include <iostream>
#include <map>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

#include <TCanvas.h>
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"



typedef std::pair<double, double> EnergyRange_MeV;
typedef std::map<const std::string, EnergyRange_MeV> MapParticlesEnergy;


class GarfieldParticle {
 public:
  GarfieldParticle(std::string particleName, double ekin_eV, double time,
                   double x_cm, double y_cm, double z_cm, double dx, double dy,
                   double dz)
      : fParticleName(particleName),
        fEkin_MeV(ekin_eV / 1000000),
        fTime(time),
        fx_mm(10 * x_cm),
        fy_mm(10 * y_cm),
        fz_mm(10 * z_cm),
        fdx(dx),
        fdy(dy),
        fdz(dz) {}
  ~GarfieldParticle() {}

  std::string getParticleName() { return fParticleName; }
  double getX_mm() { return fx_mm; }
  double getY_mm() { return fy_mm; }
  double getZ_mm() { return fz_mm; }
  double getEkin_MeV() { return fEkin_MeV; }
  double getTime() { return fTime; }
  double getDX() { return fdx; }
  double getDY() { return fdy; }
  double getDZ() { return fdz; }

 private:
  std::string fParticleName;
  double fEkin_MeV, fTime, fx_mm, fy_mm, fz_mm, fdx, fdy, fdz;
};

class GarfieldPhysics {

 public:
    const std::string label_pAnode1 =   "readout_pAnode1";         //读出极标签
    const std::string label_pGrid =    "readout_pGrid";
    const std::string label_pCathode = "readout_pCathode";

  TCanvas* cDrift= nullptr;                                    //画板
  TCanvas* Signal_pAnode =   nullptr;//电流信号
  TCanvas* Signal_pGrid =    nullptr;
  TCanvas* Signal_pCathode = nullptr;
  TCanvas* cCell = nullptr;
  TCanvas* cFieldXY = nullptr;
  TCanvas* cFieldXZ = nullptr;

  Garfield::Sensor* fSensor;                                   //recommond private:
  Garfield::AvalancheMC*  MCdrift;
  Garfield::DriftLineRKF* RKFdrift;
  Garfield::AvalancheMicroscopic* avalanche;
  Garfield::ViewDrift  driftView;
  Garfield::ViewCell* CellView;
  Garfield::ViewField* fieldViewXY;
  Garfield::ViewField* fieldViewXZ;
  Garfield::ViewSignal Signalview_pCathode;                    //  pCathode 信号可视化
  Garfield::ViewSignal Signalview_pGrid;
  Garfield::ViewSignal Signalview_pAnode;                      //  pAnode 信号可视化

public:

  static GarfieldPhysics* GetInstance();

  static void Dispose();

  void InitializePhysics();

  void ModifyLineData(std::string fileName, int lineNum, char* lineData);//修正读入的COMSOL文件


  void DoIt(std::string particleName, double ekin_MeV, double time, double x_cm,
            double y_cm, double z_cm, double dx, double dy, double dz);

  bool CoutMessageToDebug(){return CoutMessage; }

  bool WhetherToDraw(){return Draw; }

  void AddParticleName(const std::string particleName, double ekin_min_MeV,
                       double ekin_max_MeV, std::string program);
  bool FindParticleName(const std::string name,
                        std::string program = "garfield");
  bool FindParticleNameEnergy(std::string name, double ekin_MeV,
                              std::string program = "garfield");
  double GetMinEnergyMeVParticle(std::string name,
                                 std::string program = "garfield");
  double GetMaxEnergyMeVParticle(std::string name,
                                 std::string program = "garfield");
  void SetIonizationModel(std::string model, bool useDefaults = true);

  std::string GetIonizationModel();

  std::vector<GarfieldParticle*>* GetSecondaryParticles();

  void DeleteSecondaryParticles();

  inline void EnableCreateSecondariesInGeant4(bool flag) {
    createSecondariesInGeant4 = flag;
  }

  inline bool GetCreateSecondariesInGeant4() {
    return createSecondariesInGeant4;
  }

  double GetGarfieldEDep() { return fGarfieldEDep / 1000000; }

  int GetNumberOfElectron() { return NumberOfElectron; }

   //Geant4统计能量沉积的，勿改动
  inline double GetEnergyDeposit_MeV() { return fEnergyDeposit / 1000000; }

  inline double GetAvalancheSize() { return fAvalancheSize; }

  inline double GetGain() { return fGain; }

  inline void Clear() {
    fEnergyDeposit = 0;
    fGarfieldEDep = 0;
    fAvalancheSize = 0;
    fGain = 0;
    nsum = 0;
    NumberOfElectron=0;//统计电子数
  }

 private:
  GarfieldPhysics();
  ~GarfieldPhysics();

  std::string fIonizationModel;

  static GarfieldPhysics* fGarfieldPhysics;
  MapParticlesEnergy fMapParticlesEnergyGeant4;
  MapParticlesEnergy fMapParticlesEnergyGarfield;
  Garfield::MediumMagboltz* fMediumMagboltz;
  //Garfield::Sensor* fSensor;
  Garfield::TrackHeed* fTrackHeed;
  Garfield::ComponentComsol* fComponentComsol;

  std::vector<GarfieldParticle*>* fSecondaryParticles;

  bool CoutMessage;
  bool Draw;
  bool createSecondariesInGeant4;
  double fEnergyDeposit; //Geant4统计能量沉积的，勿改动
  double fGarfieldEDep;
  double fAvalancheSize;
  double fGain;
  int nsum;
  int NumberOfElectron;//统计电子数





};
#endif /* GARFIELDMODELCONFIG_HH_ */
