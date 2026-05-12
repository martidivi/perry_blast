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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

using namespace B2;

namespace B2a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{
  fMessenger = new DetectorMessenger(this);
// GLAST
  fNbOfLayers = 10;
  fLogicLayer = new G4LogicalVolume*[fNbOfLayers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete[] fLogicLayer;
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Silicon defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_Si");

  // Lead defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4Material* air = G4Material::GetMaterial("G4_AIR");
  G4Material* silicon = G4Material::GetMaterial("G4_Si");
  G4Material* lead = G4Material::GetMaterial("G4_Pb");

  // Sizes of the principal geometrical components (solids)

  G4double zCoords[10] = { 0, 41, 41+24, 2*41+24, 2*41+2*24, 3*41+2*24, 3*41+3*24, 4*41+3*24, 4*41+4*24, 5*41+4*24 };

  G4double worldLength = 1.5 * m;

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm << " mm" << G4endl;

  auto worldS = new G4Box("world",  // its name
                          worldLength / 2, worldLength / 2, worldLength / 2);  // its size
  auto worldLV = new G4LogicalVolume(worldS,  // its solid
                                     air,  // its material
                                     "World");  // its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
                                   G4ThreeVector(),  // at (0,0,0)
                                   worldLV,  // its logical volume
                                   "World",  // its name
                                   nullptr,  // its mother  volume
                                   false,  // no boolean operations
                                   0,  // copy number
                                   fCheckOverlaps);  // checking overlaps

  // Visualization attributes

  G4VisAttributes boxVisAtt(G4Colour::White());
  G4VisAttributes chamberVisAtt(G4Colour::Yellow());
  G4VisAttributes targetVisAtt(G4Colour::Blue());

  worldLV->SetVisAttributes(boxVisAtt);

  // Tracker segments
    
// GLAST
  G4double EdgeWidth = 1 * mm;
  G4double StripPitch = .228 * mm;
  G4double LadderSeparation = .220 * mm;
  G4double StripLength = 89.5 * 4 * mm;

  G4cout << "There are " << fNbOfLayers << " layers in the tracker region. " << G4endl
         << "The layers are " << StripLength / cm << " cm of Si. The distance between layers is " << (zCoords[1]-zCoords[0]) / cm << " cm or " << (zCoords[2]-zCoords[1]) / cm << " cm"<< G4endl;
    
  for (G4int copyNo = 0; copyNo < fNbOfLayers; copyNo++) {
    G4double Zposition = (zCoords[copyNo] -150) * mm;

    auto layerS = new G4Box("Plane_solid", 179.0 * mm, 179.0 * mm, 0.2 * mm);

    fLogicLayer[copyNo] =
      new G4LogicalVolume(layerS, silicon, "Chamber_LV", nullptr, nullptr, nullptr);

    fLogicLayer[copyNo]->SetVisAttributes(chamberVisAtt);

    new G4PVPlacement(nullptr,  // no rotation
                      G4ThreeVector(0, 0, Zposition),  // at (x,y,z)
                      fLogicLayer[copyNo],  // its logical volume
                      "Chamber_PV",  // its name
                      worldLV,  // its mother  volume
                      false,  // no boolean operations
                      copyNo,  // copy number
                      fCheckOverlaps);  // checking overlaps
  }
    
    G4double Ztargetposition = (zCoords[9] - 150 + 350 ) * mm;
    auto targetS = new G4Box("Target_solid", 179.0 * mm, 179.0 * mm, 3 * mm);
    fLogicTarget =
      new G4LogicalVolume(targetS, lead, "Target_LV", nullptr, nullptr, nullptr);
    fLogicTarget->SetVisAttributes(targetVisAtt);
    new G4PVPlacement(nullptr,  // no rotation
                      G4ThreeVector(0, 0, Ztargetposition),  // at (x,y,z)
                      fLogicTarget,  // its logical volume
                      "Target_PV",  // its name
                      worldLV,  // its mother  volume
                      false,  // no boolean operations
                      1,  // copy number
                      fCheckOverlaps);  // checking overlaps

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "/TrackerChamberSD";
  auto trackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD);
  // Setting trackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", trackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit) && (maxStep > 0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B2a
