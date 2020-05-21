//******************************************************************************
// PrimaryGeneratorAction.cc
//
// 1.00 JMV, LLNL, Jan-2007:  First version.
//******************************************************************************
//

#include <stdlib.h>     /* getenv */
#include <iomanip>
#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"

using namespace std;


#include "G4Event.hh"

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction(const char *inputfile)
{
  //  sk debug
  std::cout<<"PrimaryGeneratorAction is called.... (SK)"<<std::endl;
  std::cout<<"  inputfile="<<inputfile<<std::endl;
  // define a particle gun
  particleGun = new G4ParticleGun();

  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(inputfile,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    std::cout<<"(SK)  input file fail... "<<std::endl;
    if( *inputfile !=0)  //....only complain if a filename was given
      G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    std::string  dataDir="../data";
    char* pPath;
    pPath = getenv ("CRYDATAPATH");
    if (pPath!=NULL) { string ss(pPath); dataDir=ss;}

    std::cout<<"CRY::PrimaryGenratorAction  dataDIR="<<dataDir<<std::endl;
    // CRYSetup *setup=new CRYSetup(setupString,"../data");
    CRYSetup *setup=new CRYSetup(setupString,dataDir);

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
  // create a vector to store the CRY particle properties
  vect=new std::vector<CRYParticle*>;

  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();

  // Create the messenger file
  gunMessenger = new PrimaryGeneratorMessenger(this);
}

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::InputCRY()
{
  InputState=1;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::UpdateCRY(std::string* MessInput)
{

    std::string  dataDir="../data";
    char* pPath;
    pPath = getenv ("CRYDATAPATH");
    if (pPath!=NULL) { string ss(pPath); dataDir=ss;}

    std::cout<<"CRY::PrimaryGenratorAction(Update)  dataDIR="<<dataDir<<std::endl;

    //   momentum cut...
    momentumCut=0.0;   // in MeV
    char* envPcut;
    envPcut=getenv("G4CRYPCUT");
    if(envPcut!=NULL) {int ppp=atoi(envPcut); momentumCut=ppp;}

    std::cout<<"G4CRY momentum cut:  "<<momentumCut<<" MeV."<<std::endl;

    //   X,Y,Z offset...
    XSHIFT=0;
    char* envXS=0;
    envXS=getenv("G4CRYXSHIFT");
    if(envXS!=NULL) {int XS=atoi(envXS); XSHIFT=XS;}
    std::cout<<"XSHIFT is :  "<<XSHIFT<<std::endl;

    YSHIFT=0;
    char* envYS=0;
    envYS=getenv("G4CRYYSHIFT");
    if(envYS!=NULL) {int YS=atoi(envYS); YSHIFT=YS;}
    std::cout<<"YSHIFT is :  "<<YSHIFT<<std::endl;

    ZSHIFT=0;
    char* envZS=0;
    envZS=getenv("G4CRYZSHIFT");
    if(envZS!=NULL) {int ZS=atoi(envZS); ZSHIFT=ZS;}
    std::cout<<"ZSHIFT is :  "<<ZSHIFT<<std::endl;

    std::cout<<" beam X, Y, Z shift: "<<XSHIFT<<"  "<<YSHIFT<<"  "<<ZSHIFT<<"  "<<std::endl;


  // CRYSetup *setup=new CRYSetup(*MessInput,"../data");
  CRYSetup *setup=new CRYSetup(*MessInput,dataDir);

  gen = new CRYGenerator(setup);

  // set random number generator
  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState=0;

}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(newValue,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    G4cout << "Failed to open input file " << newValue << G4endl;
    G4cout << "Make sure to define the cry library on the command line" << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup=new CRYSetup(setupString,"../data");

    gen = new CRYGenerator(setup);

  // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
}
//This part of the code is what needs change!! - sas Nov 25
//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  if (InputState != 0) {
    G4String* str = new G4String("CRY library was not successfully initialized");
    //G4Exception(*str);
    G4Exception("PrimaryGeneratorAction", "1",
                RunMustBeAborted, *str);
  }
  G4String particleName;

 //  loop until we find a particle with momentum above threshold...
 //  this is practically infinit loop...
 for(int ijk=0; ijk<100000000; ijk++) {

  vect->clear();
  gen->genEvent(vect);

  std::cout<<"CRY::Simulation time in seconds: "<<gen->timeSimulated()<<std::endl; //SAS 2O20

  // apply momentum cut...
  int acceptThis=0;
  for ( unsigned j=0; j<vect->size(); j++) {
      if( (*vect)[j]->ke() > momentumCut) { acceptThis=1;} 
  }
  if(acceptThis==0) continue;

  //  high momentum particle(s) are found
  // G4double xshift=XSHIFT;
  // G4double yshift=YSHIFT;
  // G4double zshift=ZSHIFT;    // in mm

  //....debug output
/*  G4cout << "\nEvent=" << anEvent->GetEventID() << " "
         << "CRY generated nparticles=" << vect->size()
         << G4endl;
*/

  for ( unsigned j=0; j<vect->size(); j++) {
    particleName=CRYUtils::partName((*vect)[j]->id());

    //....debug output  
/*    cout << "  "          << particleName << " "
         << "charge="      << (*vect)[j]->charge() << " "
         << setprecision(4)
         << "energy (MeV)=" << (*vect)[j]->ke()*MeV << " "
         << "pos (m)"
         << G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), (*vect)[j]->z())
         << " " << "direction cosines "
         << G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w())
         << " " << endl;
*/

    // G4double xvtx=(*vect)[j]->x()*m + xshift;  
    // G4double yvtx=(*vect)[j]->y()*m + yshift;
    // G4double zvtx=(*vect)[j]->z()*m + zshift;

    // calculate the starting point of particles which end up near the detector

    G4double xvtx=(*vect)[j]->x()*m;  
    G4double yvtx=(*vect)[j]->y()*m;
    G4double zvtx=(*vect)[j]->z()*m;

  // sizes are copyr from B4DetectorConstruction.cc...
  auto worldSizeX = 6000.0*cm ;
  auto worldSizeY = 3000.0*cm ;
  auto worldSizeZ = 6000.0*cm ;  
    
    // move particle back to the world boundary...
    double deltaL=worldSizeZ/1000.0;
    for (int i=0; i<100000; i++) {
      double x= xvtx-(*vect)[j]->u() * deltaL;
      double y= yvtx-(*vect)[j]->v() * deltaL;
      double z= zvtx-(*vect)[j]->w() * deltaL;
      if(abs(x) > worldSizeX ) break;
      if(abs(y) > worldSizeY ) break;
      if(abs(z) > worldSizeZ ) break;
      xvtx=x; yvtx=y; zvtx=z;
    }

    // std::cout<<"  xvtx "<<xvtx<<"  yvtx  "<<yvtx<<"  zvtx  "<<zvtx<<std::endl;

    particleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
    particleGun->SetParticleEnergy((*vect)[j]->ke()*MeV);
    // particleGun->SetParticlePosition(G4ThreeVector((*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m));
    particleGun->SetParticlePosition(G4ThreeVector(xvtx, yvtx, zvtx));
    particleGun->SetParticleMomentumDirection(G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w()));
    particleGun->SetParticleTime((*vect)[j]->t());
    particleGun->GeneratePrimaryVertex(anEvent);
    delete (*vect)[j];
  }

  //  exist from infinit loop...
  break;

 }   // end of ijk loop...

}
