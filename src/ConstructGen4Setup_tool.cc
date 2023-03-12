///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "art/Utilities/ToolMacros.h"

// Mu2e includes.

// #include "Mu2eG4/inc/ConstructMaterials.hh"
// #include "Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/InitEnvToolBase.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"
// #include "Mu2eG4/inc/findMaterialOrThrow.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4NistManager.hh"
#include "Geant4/G4UserLimits.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4RotationMatrix.hh"

#include "Geant4/G4ParticleTable.hh"
// #include "Geant4/G4ProcessManager.hh"
// #include "Geant4/G4StepLimiter.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class ConstructGen4Setup: public InitEnvToolBase {
  public:
    ConstructGen4Setup (const fhicl::ParameterSet& PSet);
    ~ConstructGen4Setup();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructGen4Setup::ConstructGen4Setup(const fhicl::ParameterSet& PSet) {
    _name = "Gen4Setup";
  }

//-----------------------------------------------------------------------------
  ConstructGen4Setup::~ConstructGen4Setup() {
  }

//-----------------------------------------------------------------------------
// geometry information passed through the config file
  int ConstructGen4Setup::construct(VolumeInfo const& parentVInfo, SimpleConfig const& Config) {
    // 2016-11-06 P.Murat: if needed, add GaAs as material

    G4NistManager* nist = G4NistManager::Instance();
    nist->FindOrBuildMaterial("G4_Al",true,true);
    nist->FindOrBuildMaterial("G4_Au",true,true);
    
    G4Material* aluminum = G4Material::GetMaterial("G4_Al",false);
    //    G4Material* gold     = G4Material::GetMaterial("G4_Au",false);

    G4Material* gaas = G4Material::GetMaterial("GaAs",false);

    if (gaas == NULL) {
      gaas = new G4Material("GaAs", 5.3*CLHEP::g/CLHEP::cm3, 2 );

      G4Element* ga = nist->FindOrBuildElement("Ga",true);
      gaas->AddElement( ga, 1);

      G4Element* as = nist->FindOrBuildElement("As",true);
      gaas->AddElement( as, 1);
    }

    // polyethylene

    G4Material* poly = G4Material::GetMaterial("Polyethylene",false);
    if (poly == NULL) {
      poly = new G4Material("Polyethylene", 0.956*CLHEP::g/CLHEP::cm3, 2);
      G4Element* carbon = nist->FindOrBuildElement("C",true);
      poly->AddElement(carbon, 1);
      G4Element* hydrogen = nist->FindOrBuildElement("H",true);
      poly->AddElement(hydrogen, 2);
    }

    G4Material* steel = G4Material::GetMaterial("StainlessSteel",false);
    if (steel == NULL) {
      steel = new G4Material("StainlessSteel", 8.02*CLHEP::g/CLHEP::cm3, 5);

      steel->AddMaterial(nist->FindOrBuildMaterial("G4_Mn",true,true), 0.02);
      steel->AddMaterial(nist->FindOrBuildMaterial("G4_Si",true,true), 0.01);
      steel->AddMaterial(nist->FindOrBuildMaterial("G4_Cr",true,true), 0.19);
      steel->AddMaterial(nist->FindOrBuildMaterial("G4_Ni",true,true), 0.10);
      steel->AddMaterial(nist->FindOrBuildMaterial("G4_Fe",true,true), 0.68);
    }

    const bool forceAuxEdgeVisible = Config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = Config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    G4Colour  orange  (.75, .55, .0);
//-----------------------------------------------------------------------------
// construct source holder and the source - several tubes
// source holder - an external polyethylene tube
//-----------------------------------------------------------------------------
    G4bool arodVisible    = Config.getBool("arod.visible",true);
    G4bool arodSolid      = Config.getBool("arod.solid",true);
//-----------------------------------------------------------------------------
// aluminum rod
//-----------------------------------------------------------------------------
    double arod_halfLength = Config.getDouble("arod.halfLength");

    TubsParams arodParams( Config.getDouble("arod.rIn"),
                           Config.getDouble("arod.rOut"),
                           arod_halfLength,
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );
//-----------------------------------------------------------------------------
// source position is defined in the .txt file, nominally - the particleGun
// input file, thus 'particleGun.point'
//-----------------------------------------------------------------------------
    const G4ThreeVector source_pos(Config.getHep3Vector("particleGun.point"));

    const G4ThreeVector arodCenterInWorld(source_pos.x(),source_pos.y(),source_pos.z()+arod_halfLength);
    
    VolumeInfo arodVInfo(nestTubs( "ARod",
                                   arodParams,
                                   aluminum,
                                   0,                                // no rotation,
                                   arodCenterInWorld,
                                   parentVInfo,
                                   Config.getInt("arod.copyNumber"), // non-0 for tracking purposes
                                   arodVisible,
                                   G4Color::Blue(),
                                   arodSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// the source itself is an infinitely thin, coating layer
//-----------------------------------------------------------------------------
    TubsParams coatParams( Config.getDouble("coat.rIn"),
                           Config.getDouble("coat.rOut"),
                           Config.getDouble("coat.halfThickness"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    const G4ThreeVector coatCenterInWorld(arodCenterInWorld[0],
					  arodCenterInWorld[1],
					  arodCenterInWorld[2]-arod_halfLength-coatParams.zHalfLength());
  
    G4Material* coatMaterial = G4Material::GetMaterial(Config.getString("coat.material"),false);
  
    VolumeInfo coatVInfo(nestTubs( "Coat",
                                   coatParams,
                                   coatMaterial,
                                   0,                                // rotation,
                                   coatCenterInWorld,
                                   parentVInfo,
                                   Config.getInt("coat.copyNumber"), // non-0 for tracking purposes
                                   arodVisible,
                                   G4Color::Blue(),
                                   arodSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// wall: 0.5 mm deep hole
//-----------------------------------------------------------------------------
    TubsParams wallParams( Config.getDouble("wall.rIn"),
                           Config.getDouble("wall.rOut"),
                           Config.getDouble("wall.halfLength"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    const G4ThreeVector wallCenterInWorld(arodCenterInWorld[0],
					  arodCenterInWorld[1],
					  arodCenterInWorld[2]-arod_halfLength-wallParams.zHalfLength());

    VolumeInfo wallVInfo(nestTubs( "Wall",
                                   wallParams,
                                   aluminum,
                                   0,                                // no rotation,
                                   wallCenterInWorld,
                                   parentVInfo,
                                   Config.getInt("wall.copyNumber"), // non-0 for tracking purposes
                                   arodVisible,
                                   G4Color::Blue(),
                                   arodSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// "collimator" 
//-----------------------------------------------------------------------------
    TubsParams collParams( Config.getDouble("coll.rIn"),
                           Config.getDouble("coll.rOut"),
                           Config.getDouble("coll.halfLength"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    const G4ThreeVector collCenterInWorld(arodCenterInWorld[0],
					  arodCenterInWorld[1],
					  wallCenterInWorld[2]-wallParams.zHalfLength()-collParams.zHalfLength());

    VolumeInfo collVInfo(nestTubs( "Coll",
                                   collParams,
                                   poly,
                                   0,                                // rotation,
                                   collCenterInWorld,
                                   parentVInfo,
                                   Config.getInt("coll.copyNumber"), // non-0 for tracking purposes
                                   arodVisible,
                                   G4Color::Blue(),
                                   arodSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// virtual detector just below zero
//-----------------------------------------------------------------------------
    TubsParams vd1Params( Config.getDouble("vd1.rIn"),
			  Config.getDouble("vd1.rOut"),
			  Config.getDouble("vd1.halfLength"),
			  0.*CLHEP::degree,
			  360.*CLHEP::degree );

    const G4ThreeVector vd1CenterInWorld(Config.getHep3Vector("vd1.centerInWorld"));

    G4Material* vd1Material = G4Material::GetMaterial(Config.getString("vd1.material"),false);

    VolumeInfo vd1VInfo(nestTubs( "Vd1",
                                   vd1Params,
                                   vd1Material,
                                   0,                               // rotation,
                                   vd1CenterInWorld,
                                   parentVInfo,
                                   Config.getInt("vd1.copyNumber"), // non-0 for tracking purposes
                                   arodVisible,
                                   G4Color::Blue(),
                                   arodSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// construct and position sensor
//-----------------------------------------------------------------------------
    int    n_photodiodes        = Config.getInt ("nPhotodiodes"  , 0  );
    G4bool sensorVisible        = Config.getBool("sensor.visible",true);
    G4bool sensorSolid          = Config.getBool("sensor.solid"  ,true);

    vector<double> sensorHalfLengths;
    Config.getVectorDouble( "sensor.halfLengths", sensorHalfLengths);

    const G4ThreeVector sensorCenterInWorld(Config.getHep3Vector("sensor.centerInWorld"));

    VolumeInfo sensor(nestBox("GaAsSensor",
                              sensorHalfLengths,
                              gaas,
                              0,                               // no rotation
                              sensorCenterInWorld,
                              parentVInfo,
                              Config.getInt("sensor.copyNumber"), // for tracking purposes
                              sensorVisible,
                              orange,
                              sensorSolid,
                              forceAuxEdgeVisible,
                              placePV,
                              doSurfaceCheck
                              )
                      );
//-----------------------------------------------------------------------------
// construct and position the photodiodes - all 7 of them
// PD offset - offset of the PD center wrt the scintillator center
//-----------------------------------------------------------------------------
    for (int i=0; i<n_photodiodes; i++) {

      G4bool pdVisible        = Config.getBool(Form("pd%i.visible",i),true);
      G4bool pdSolid          = Config.getBool(Form("pd%i.solid"  ,i),true);

      vector<double> pd_half_length;
      Config.getVectorDouble(Form("pd%i.halfLength",i),pd_half_length);

      vector<double> pd_offset;
      Config.getVectorDouble(Form("pd%i.offset",i),pd_offset);

      double xc = sensorCenterInWorld.x()+pd_offset[0];
      double yc = sensorCenterInWorld.y()+pd_offset[1];
      double zc = sensorCenterInWorld.z()+pd_offset[2];

      int pd_copy_number = Config.getInt(Form("pd%i.copyNumber",i));

      G4ThreeVector pd_pos(xc,yc,zc);

      VolumeInfo pd (nestBox( Form("InGaAsPD%i",i),
                              pd_half_length,
                              gaas,
                              0,                     // 0: no rotation
                              pd_pos,
                              parentVInfo,
                              pd_copy_number,        // for volume tracking purposes
                              pdVisible,
                              G4Color::Red(),
                              pdSolid,
                              forceAuxEdgeVisible,
                              placePV,
                              doSurfaceCheck
                              )
                    );
    }
//-----------------------------------------------------------------------------
// now limit the max step
//-----------------------------------------------------------------------------
    double max_g4_step = Config.getDouble("detector.max_g4_step");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    G4UserLimits* stepLimit = reg.add(G4UserLimits(max_g4_step));

    mu2e::Mu2eG4Helper* mgh = art::ServiceHandle<mu2e::Mu2eG4Helper>().get();

    boost::regex expr1("GaAsSensor");
    std::vector<mu2e::VolumeInfo const*> v1 = mgh->locateVolInfo(expr1);

    printf("%s: v1.size() : %i\n",__func__,(int) v1.size());
    for ( auto v : v1 ){
      v->logical->SetUserLimits( stepLimit );
      printf(" Activated step limit for volume %s\n",v->logical->GetName().data());
    }

    for (int i=0; i<n_photodiodes; i++) {
      boost::regex expr2(Form("InGaAsPD%i",i));
      std::vector<mu2e::VolumeInfo const*> v2 = mgh->locateVolInfo(expr2);
    
      printf("%s: regexp:%s: v2.size() : %i\n",__func__,Form("InGaAsPD%i",i),(int) v2.size());

      for ( auto v : v2 ){
        v->logical->SetUserLimits( stepLimit );
        printf(" Activated step limit for volume %s\n",v->logical->GetName().data());
      }
    }

    // sensor.logical->SetUserLimits(stepLimit);

    // for (auto pd : list_of_pd) {
    //   pd.logical->SetUserLimits(stepLimit);
    // }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructGen4Setup)
