#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(30.*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
        G4double cone_beam_angle = 10.*M_PI/180.;
        G4double source_size = 1e-2*cm;

        //G4double source_size = 1*cm;

	G4double dz=61.*cm;
	G4double boundary = dz*tan(cone_beam_angle);

	G4double x0=0.;
	G4double y0=0.;
	G4double z0=-60.*cm;

	G4double px,py,pz;

	G4double theta;
	G4double pi;
	
	G4int flag =1;
	while(flag)
	{
		theta = M_PI/2.+2.*(cone_beam_angle)*(G4UniformRand()-0.5);
		pi = M_PI/2.+2.*(cone_beam_angle)*(G4UniformRand()-0.5);

		G4double x = dz*sin(theta)*cos(pi);
		G4double y = dz*cos(theta);
		G4double z = dz*sin(theta)*sin(pi);

		if(x*x+y*y<=((dz*sin(cone_beam_angle))*dz*sin(cone_beam_angle)))
		{
			px = x;
			py = y;
			pz = z;

			flag=0;
		}
	}



        G4ThreeVector source_pos = G4ThreeVector(x0,y0,z0);
        //G4ThreeVector source_pos = G4ThreeVector(xf,yf,zf);

  fParticleGun->SetParticlePosition(source_pos);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

