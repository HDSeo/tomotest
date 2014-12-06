#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  G4double world_size=2.*m;
  G4bool checkOverlaps = true;

//////////////////////////////////////////////////////////////

	G4Element* Ca = new G4Element("Calcium","Ca",20,40.078*g/mole);
	G4Element* C = new G4Element("Carbon","C",6,12.0107*g/mole);
	G4Element* O = new G4Element("Oxygen","O",8,15.9994*g/mole);

	G4Material* CaCO3 =
	new G4Material("CaCO3",2.83*g/cm3,3);

	CaCO3->AddElement(Ca,20.*perCent);
	CaCO3->AddElement(C,20.*perCent);
	CaCO3->AddElement(O,60.*perCent);

//////////////////////////////////////////////////////////////


  //     
  // World
  //
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_WATER");

  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_size,0.5*world_size, 0.5*world_size);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        env_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  

  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*world_size, 0.5*world_size, 0.5*world_size); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  //     
  // Detector
  //  
	
  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_Galactic");

	G4double detector_sizeY=30.*cm;
	G4double detector_sizeX=30.*cm;
	G4double detector_sizeZ=2.*cm;
	G4ThreeVector detector_pos= G4ThreeVector(0,0,3.1*cm);

  G4Box* detector =    
    new G4Box("Detector",                    //its name
        0.5*detector_sizeX, 0.5*detector_sizeY, 0.5*detector_sizeZ); //its size
      
  G4LogicalVolume* detector_logic =                         
    new G4LogicalVolume(detector,            //its solid
                        detector_mat,             //its material
                        "Detector");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    detector_pos,         //at (0,0,0)
                    detector_logic,                //its logical volume
                    "Detector",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
										checkOverlaps);

	// breast
	//
/*
	G4Material* breast_mat = nist->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");

	G4Material* breasto_mat = nist->FindOrBuildMaterial("G4_Pb");
	G4ThreeVector breast_pos = G4ThreeVector((0+10)*cm,0,0);

	G4RotationMatrix rotm = G4RotationMatrix();
	rotm.rotateY(90.*deg);
	
	G4Transform3D transform = G4Transform3D(rotm,breast_pos);




	G4double xsemi =2.2*cm;
	G4double ysemi = 6.*cm;
	G4double zsemi = 6.*cm;
	G4double zbottomcut = 0.*cm;
	G4double ztopcut = 6.*cm;
	

	G4Ellipsoid* breast_2 =
		new G4Ellipsoid("Breast2",
										xsemi,
										ysemi,
										zsemi,
										zbottomcut,
										ztopcut);

	G4LogicalVolume* breast2_logic =
		new G4LogicalVolume(breast_2,breast_mat,"Breast2");

	new G4PVPlacement(transform,
										breast2_logic,
										"Breast2",
										logicEnv,
										false,
										0,
										checkOverlaps);


	G4double inner =0*cm;
	G4double outer = 0.07978845608028654*cm;
	//G4double outer = 2.57978845608028654*cm;
	G4double SPhi = 0*deg;
	G4double DPhi = 360*deg;
	G4double STheta = 0*deg;
	G4double DTheta = 360*deg;


	G4Sphere* breast =
		new G4Sphere("Breast",
								inner,
								outer,
								SPhi,
								DPhi,
								STheta,
								DTheta);

	G4LogicalVolume* breast_logic =
		//new G4LogicalVolume(breast,CaCO3,"Breast");
		new G4LogicalVolume(breast,breasto_mat,"Breast");
	new G4PVPlacement(0,
										G4ThreeVector(0*cm,0*cm,(1.043+3)*cm),
										breast_logic,
										"Breast",
										breast2_logic,
										false,
										0,
										checkOverlaps);
*/

  fScoringVolume = logicEnv;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
