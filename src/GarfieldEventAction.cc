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
// $Id: GarfieldEventAction.cc 999993 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldEventAction.cc
/// \brief Implementation of the GarfieldEventAction class

#include "GarfieldEventAction.hh"

#include <iomanip>

#include <TFile.h>

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "GarfieldAnalysis.hh"
#include "GarfieldPhysics.hh"
#include "GarfieldRunAction.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldEventAction::GarfieldEventAction()
    : G4UserEventAction(), fEnergyAbs(0.), fEnergyGas(0.), fTrackLAbs(0.) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldEventAction::~GarfieldEventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double GarfieldEventAction::transfer(double t) {
  const G4double tau = 2000; // peaking time in ns
  const G4double fC_to_mV = 1/3.; //is equlal to 1/Cf  ---> mV per fC
  // Impulse response of the amplifier.
  return -fC_to_mV *exp(- t / tau);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldEventAction::BeginOfEventAction(const G4Event* /*event*/) {
  // initialisation per event
  fEnergyAbs = 0;
  fEnergyGas = 0;
  fTrackLAbs = 0;
  fAvalancheSize = 0;
  fGain = 0;

  GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();
  garfieldPhysics->fSensor->ClearSignal();
  garfieldPhysics->Clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldEventAction::EndOfEventAction(const G4Event* event) {
  // Accumulate statistics
  //
  GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  fAvalancheSize = garfieldPhysics->GetAvalancheSize();
  fGain = garfieldPhysics->GetGain();

  fThisEventID=event->GetEventID();//  Returns the event ID

  if(garfieldPhysics->WhetherToDraw()){
          TFile *hist = new TFile("Results.root", "RECREATE");
          garfieldPhysics->Signal_pAnode->Write("Anode1 Signal");
          garfieldPhysics->Signal_pGrid->Write("Grid Signal");
          garfieldPhysics->Signal_pCathode->Write("Cathode signal");
          garfieldPhysics->cCell->Write("cCell");
          garfieldPhysics->cDrift->Write("Drift");
          garfieldPhysics->cFieldXY->Write("cFieldXY");
          garfieldPhysics->cFieldXZ->Write("cFieldXZ");
          hist->Close();

          static std::ofstream out;//电流信号输出到文本检查
          out.open("CurrentSignalInEvent.txt", std::ios::out);
          out <<"Time(ns) AnodeSignal CathodeSignal GridSignal(fC)"<< "\n";
          for (unsigned int i = 0; i < 3000; ++i) {
            const double t = (i + 0.5) ;
            const double a = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pAnode1, i);
            const double c = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pCathode, i);
            const double g = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pGrid, i);
            out << t << "  " << a << "  " <<  c<< "  " <<g  << "\n";
          }
  }

  for(int i=0;i<3000;i++){//电流信号
      fAnodeSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pAnode1, i);
      fCathodeSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pCathode, i);
      fGridSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pGrid, i);
  }


  if(garfieldPhysics->WhetherToDraw()){//卷积电流信号并输出到文本
      garfieldPhysics->fSensor->SetTransferFunction(transfer);//卷积电流信号
      garfieldPhysics->fSensor->ConvoluteSignals(true);
      for(int i=0;i<3000;i++){//电压信号
          fAnodeVoltageSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pAnode1, i);
          fCathodeVoltageSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pCathode, i);
          fGridVoltageSignal[i]=garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pGrid, i);
      }
          static std::ofstream out2;//电压信号输出到文本检查
          out2.open("VoltageSignalInEvent.txt", std::ios::out);
          out2 <<"Time(ns) AnodeSignal CathodeSignal GridSignal"<< "\n";
          for (unsigned int i = 0; i < 3000; ++i) {
            const double t = (i + 0.5) ;
            const double a = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pAnode1, i);
            const double c = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pCathode, i);
            const double g = garfieldPhysics->fSensor->GetSignal(garfieldPhysics->label_pGrid, i);
            out2 << t << "  " << a << "  " <<  c<< "  " <<g  << "\n";
          }
  }

  // fill histograms
  analysisManager->FillH1(1, fEnergyAbs);
  analysisManager->FillH1(2, fTrackLAbs);
  analysisManager->FillH1(3, fEnergyGas);
  analysisManager->FillH1(4, fAvalancheSize);
  analysisManager->FillH1(5, fGain);

// fill ntuple
//  analysisManager->FillNtupleDColumn(1,0, fEnergyAbs);
//  analysisManager->FillNtupleDColumn(1,1, fTrackLAbs);
//  analysisManager->FillNtupleDColumn(1,2, fEnergyGas);
//  analysisManager->FillNtupleDColumn(1,3, fAvalancheSize);
//  analysisManager->FillNtupleDColumn(1,4, fGain);
//  analysisManager->AddNtupleRow(1);  //相当于 Fill


  for(int i=0;i<3000;i++){
    analysisManager->FillNtupleDColumn(0, fThisEventID);
    analysisManager->FillNtupleDColumn(1, fAnodeSignal[i]);
    analysisManager->FillNtupleDColumn(2, fCathodeSignal[i]);
    analysisManager->FillNtupleDColumn(3, fGridSignal[i]);
    analysisManager->AddNtupleRow();  //相当于 Fill
  }



  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ((printModulo > 0) && (eventID % printModulo == 0)) {
    G4cout << "---> End of event: " << eventID << G4endl;

    G4cout << "   Absorber: total energy: " << std::setw(7)
           << G4BestUnit(fEnergyAbs, "Energy")
           << "       total track length: " << std::setw(7)
           << G4BestUnit(fTrackLAbs, "Length") << G4endl;

    G4cout << "        Gas: total energy: " << std::setw(7)
           << G4BestUnit(fEnergyGas, "Energy")
           << "       avalanche size: " << fAvalancheSize
           << "       gain: " << fGain << G4endl;
    G4cout << "GarfieldEDe: "<<
          garfieldPhysics->GetGarfieldEDep()<<" MeV"<<G4endl;
    G4cout << "SimulationW: "<<
          3000*5.3160*1e6/garfieldPhysics->GetNumberOfElectron()/26.49<<G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
