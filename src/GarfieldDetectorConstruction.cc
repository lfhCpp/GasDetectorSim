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
// $Id: GarfieldDetectorConstruction.cc 999992 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldDetectorConstruction.cc
/// \brief Implementation of the GarfieldDetectorConstruction class

#include "GarfieldDetectorConstruction.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "GarfieldG4FastSimulationModel.hh"
#include "GarfieldMessenger.hh"

#include "CADMesh.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldDetectorConstruction::GarfieldDetectorConstruction()
    : G4VUserDetectorConstruction(),fCheckOverlaps(true) {
  fGarfieldMessenger = new GarfieldMessenger(this);//初始化列表：信使
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldDetectorConstruction::~GarfieldDetectorConstruction() {
  delete fGarfieldMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GarfieldDetectorConstruction::Construct() {
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldDetectorConstruction::DefineMaterials() {
  G4bool isotopes = false;
  G4String name, symbol;
  G4int ncomponents, natoms;
  G4double density, fractionmass;

  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();

  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Cu");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Au");
  nistManager->FindOrBuildMaterial("G4_W");

  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Galactic");

  G4Element* H = nistManager->FindOrBuildElement("H", isotopes);
  G4Element* N = nistManager->FindOrBuildElement("N", isotopes);
  G4Element* C = nistManager->FindOrBuildElement("C", isotopes);
  G4Element* O = nistManager->FindOrBuildElement("O", isotopes);
  G4Element* Ar = nistManager->FindOrBuildElement("Ar", isotopes);

  //Co2Density:1.997g/L Ch4Density: 0.717g/L ArDensity:1.7800g/L STP  p10 1.59mg/cm3

  G4Material* CO2 = new G4Material(
      "CO2", density = 1.977 * CLHEP::mg / CLHEP::cm3, ncomponents = 2);
  CO2->AddElement(C, natoms = 1);
  CO2->AddElement(O, natoms = 2);

  G4Material* ArCO2_70_30 =
      new G4Material("ArCO2_70_30", density = 1.8223 * CLHEP::mg / CLHEP::cm3,
                     ncomponents = 2, kStateGas);
  ArCO2_70_30->AddElement(Ar, fractionmass = 0.70);  //fractionmass warnring!!!!!!!!!!!!!!
  ArCO2_70_30->AddMaterial(CO2, fractionmass = 0.30);//fractionmass warnring!!!!!!!!!!!!!!

  G4Material* ArCh4_90_10 =
        new G4Material("ArCh4_90_10", density = 1.56 * CLHEP::mg / CLHEP::cm3,
                       ncomponents = 3);//, kStateGas);
  ArCh4_90_10->AddElement(Ar, 90);
  ArCh4_90_10->AddElement(C, 10);
  ArCh4_90_10->AddElement(H, 40);


  density = 1.413 * CLHEP::g / CLHEP::cm3;
  G4Material* Kapton =
      new G4Material(name = "Kapton", density, ncomponents = 4);
  Kapton->AddElement(O, 5);
  Kapton->AddElement(C, 22);
  Kapton->AddElement(N, 2);
  Kapton->AddElement(H, 10);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GarfieldDetectorConstruction::DefineVolumes() {
  // Geometry parameters
  G4double worldSizeXYZ = 2000 * mm;

  G4double GridWireRadius = 0.075/2 * mm;
  G4double DistanceBetweenGrid = 0.75 * mm;
  G4double ChamberHalfLengthZ = 69.9/2 * mm;
  G4double ChamberHalfLengthX = 87./2 * mm;
  G4double ChamberHalfLengthY = 147/2. * mm;
  G4double gap = 6.5 * mm;
  G4double AnodeThickness = 2 * mm;
  G4double AnodeSegmentZ = 45 * mm;
  G4double AnodeSegmentX = 100 * mm;

  // Get materials
  //G4Material* defaultMaterial = G4Material::GetMaterial("G4_Galactic");
  G4Material* defaultMaterial = G4Material::GetMaterial("G4_Pb");
  fAbsorberMaterial = G4Material::GetMaterial("G4_Cu");
  G4Material* gasMaterial_backup = G4Material::GetMaterial("ArCO2_70_30");
  G4Material* gasMaterial = G4Material::GetMaterial("ArCh4_90_10");
  G4Material* cathodeMaterial = G4Material::GetMaterial("G4_Al");
  G4Material* wireMaterial = G4Material::GetMaterial("G4_W");

  if (!defaultMaterial || !fAbsorberMaterial || !gasMaterial ||
      !cathodeMaterial || !wireMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("GarfieldDetectorConstruction::DefineVolumes()",
                "exampleGarfield", FatalException, msg);
  }



  //
  // World
  //
  G4VSolid* worldS = new G4Box("World",  // its name
                               0.5 * worldSizeXYZ, 0.5 * worldSizeXYZ,
                               0.5 * worldSizeXYZ);  // its size

  G4LogicalVolume* worldLV =
      new G4LogicalVolume(worldS,           // its solid
                          gasMaterial,  // its material
                          "World");         // its name

  G4VPhysicalVolume* worldPV =
      new G4PVPlacement(0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        worldLV,          // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps

  //
  //ChamberCAD
  //
  auto Cage = CADMesh::TessellatedMesh::FromSTL("./ComsolFile/F4MUSIC.STL");
  Cage->SetScale(1.0);// 放缩，默认单位是1mm
  G4VSolid* solid = Cage->GetSolid();
  
  G4LogicalVolume* ChamberCAD =
      new G4LogicalVolume(solid,       // its solid
                          fAbsorberMaterial,  // its material
                          "ChamberCAD");     // its name

  fChamberPV = new G4PVPlacement(
      0,
      G4ThreeVector(0*cm, 0*cm, 0*cm),                 // its position
      ChamberCAD,                                       // its logical volume
      "ChamberCAD",                                       // its name
      worldLV,                                         // its mother  volume
      false,                                           // no boolean operation
      0,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps

  //
  //ChamberSD
  //
  G4VSolid* ChamberBox = new G4Box("ChamberBox",         // its name
                                   ChamberHalfLengthX, ChamberHalfLengthX,
                                   ChamberHalfLengthZ);  // its size

  G4LogicalVolume* ChamberLV =
      new G4LogicalVolume(ChamberBox,       // its solid
                          gasMaterial,  // its material
                          "Chamber");     // its name


  fChamberPV = new G4PVPlacement(
      0,
      G4ThreeVector(ChamberHalfLengthX+gap, ChamberHalfLengthY+gap,
                    ChamberHalfLengthZ+0.05*mm),               // its position
      ChamberLV,                                       // its logical volume
      "Chamber",                                       // its name
      worldLV,                                         // its mother  volume
      false,                                           // no boolean operation
      0,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps



  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* VisAttBlue = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5));
  G4VisAttributes* VisAttGreen = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  G4VisAttributes* VisAttRed = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  //G4VisAttributes* VisAttWhite = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));

  G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.8, 0.5)); // RGBA（A=透明度，0=全透明）
  //visAtt->SetForceWireframe(true);  // 强制显示为线框
  visAtt->SetForceAuxEdgeVisible(true); // 显示辅助边
  visAtt->SetLineStyle(G4VisAttributes::dashed); // 设置线型（dashed/dotted）
  ChamberCAD->SetVisAttributes(visAtt);

  //VisAttRed->SetVisibility(false);//不显示边框
  //VisAttRed->SetForceWireframe(true);
  //VisAttRed->SetForceLineSegmentsPerCircle(360);
  VisAttRed->SetForceSolid(true);
  ChamberLV->SetVisAttributes(VisAttGreen);

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  G4Region* regionGarfield = new G4Region("RegionGarfield");               ///
  regionGarfield->AddRootLogicalVolume(ChamberLV);                         ///
                                                                           ///
  //G4Region* regionWire = new G4Region("RegionWire");                     ///
  //regionWire->AddRootLogicalVolume(wireLV);                              ///
                                                                           ///
  fGarfieldG4FastSimulationModel = new GarfieldG4FastSimulationModel(      ///
      "GarfieldG4FastSimulationModel", regionGarfield);                    ///
                                                                           ///
  fGarfieldG4FastSimulationModel->WriteGeometryToGDML(fChamberPV);         ///
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* GarfieldDetectorConstruction::AbsorberMaterialWithSingleIsotope(
    G4String name, G4String symbol, G4double density, G4int Z, G4int A) {
  // define a material from an isotope
  //
  G4int ncomponents;
  G4double abundance, massfraction;

  G4Isotope* isotope = new G4Isotope(symbol, Z, A);

  G4Element* element = new G4Element(name, symbol, ncomponents = 1);
  element->AddIsotope(isotope, abundance = 100. * perCent);

  G4Material* material = new G4Material(name, density, ncomponents = 1);
  material->AddElement(element, massfraction = 100. * perCent);

  return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldDetectorConstruction::SetAbsorberMaterial(
    G4String materialChoice) {
  // search the material by its name
  G4Material* newMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (newMaterial) {
    if (fAbsorberMaterial != newMaterial) {
      fAbsorberMaterial = newMaterial;
      if (fAbsorberLV) {
        fAbsorberLV->SetMaterial(fAbsorberMaterial);
      }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from GarfieldDetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
