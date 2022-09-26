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

  class ConstructGen3Setup: public InitEnvToolBase {
  public:
    ConstructGen3Setup (const fhicl::ParameterSet& PSet);
    ~ConstructGen3Setup();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructGen3Setup::ConstructGen3Setup(const fhicl::ParameterSet& PSet) {
    _name = "Gen3Setup";
  }

//-----------------------------------------------------------------------------
  ConstructGen3Setup::~ConstructGen3Setup() {
  }

//-----------------------------------------------------------------------------
// geometry information passed through the config file
  int ConstructGen3Setup::construct(VolumeInfo const& parentVInfo, SimpleConfig const& Config) {
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

    // polyethilene

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
// the source itself is infinitely thin, coating layer
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
                                   0, // rotation,
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
    G4bool sensorVisible        = Config.getBool("sensor.visible",true);
    G4bool sensorSolid          = Config.getBool("sensor.solid",true);

    vector<double> sensorParams;
    Config.getVectorDouble( "sensor.halfLengths", sensorParams);

    const G4ThreeVector sensorCenterInWorld(Config.getHep3Vector("sensor.centerInWorld"));

    VolumeInfo box(nestBox( "GaasSensor",
			    sensorParams,
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
			    ));
//-----------------------------------------------------------------------------
// construct and position the photodiodes - all 7 of them
//-----------------------------------------------------------------------------
    G4bool pdVisible        = Config.getBool("pd0.visible",true);
    G4bool pdSolid          = Config.getBool("pd0.solid",true);

    vector<double> pd0_params;
    Config.getVectorDouble("pd0.halfLengths",pd0_params);

    double pd_half_thickness = pd0_params[2];
    double pd_z              = sensorCenterInWorld.z()+sensorParams[2]+pd_half_thickness;

    G4ThreeVector pd0_pos;
    pd0_pos.set(sensorCenterInWorld.x()-sensorParams[0]+0.285+pd0_params[0],
		sensorCenterInWorld.y(),pd_z);

    VolumeInfo pd0 (nestBox( "InGaAsPD0",
			     pd0_params,
			     gaas,
			     0,                                 // no rotation
			     pd0_pos,
			     parentVInfo,
			     Config.getInt("pd0.copyNumber"), // for volume tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #1-#6 positioned in Y in 120 um from the lower edge
//-----------------------------------------------------------------------------
    vector<double> pd1_params;
    Config.getVectorDouble("pd1.halfLengths",pd1_params);

    double pd16_ymin = sensorCenterInWorld.y()-sensorParams[1]+0.120;

    G4ThreeVector pd1_pos;
    pd1_pos.set(pd0_pos[0]+pd0_params[0]+0.070+pd1_params[0],pd16_ymin+pd1_params[1],pd_z);

    VolumeInfo pd1 (nestBox( "InGaAsPD1",
			     pd1_params,
			     gaas,
			     0,                                 // no rotation
			     pd1_pos,
			     parentVInfo,
			     Config.getInt("pd1.copyNumber"), // assign copy nuber for volume tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #2
//-----------------------------------------------------------------------------
    vector<double> pd2_params;
    Config.getVectorDouble("pd2.halfLengths",pd2_params);

    G4ThreeVector pd2_pos;
    pd2_pos.set(pd1_pos[0]+pd1_params[0]+0.10+pd2_params[0],pd16_ymin+pd2_params[1],pd_z);

    VolumeInfo pd2 (nestBox( "InGaAsPD2",
			     pd2_params,
			     gaas,
			     0,                                 // no rotation
			     pd2_pos,
			     parentVInfo,
			     Config.getInt("pd2.copyNumber"), // for tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #3
//-----------------------------------------------------------------------------
    vector<double> pd3_params;
    Config.getVectorDouble("pd3.halfLengths",pd3_params);

    G4ThreeVector pd3_pos;
    pd3_pos.set(pd2_pos[0]+pd2_params[0]+0.10+pd3_params[0],pd16_ymin+pd3_params[1],pd_z);

    VolumeInfo pd3 (nestBox( "InGaAsPD3",
			     pd3_params,
			     gaas,
			     0,                                 // no rotation
			     pd3_pos,
			     parentVInfo,
			     Config.getInt("pd3.copyNumber"), // assign copy nuber for volume tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #4
//-----------------------------------------------------------------------------
    vector<double> pd4_params;
    Config.getVectorDouble("pd4.halfLengths",pd4_params);

    G4ThreeVector pd4_pos;
    pd4_pos.set(pd3_pos[0]+pd3_params[0]+0.10+pd4_params[0],pd16_ymin+pd4_params[1],pd_z);

    VolumeInfo pd4 (nestBox( "InGaAsPD4",
			     pd4_params,
			     gaas,
			     0,                                 // no rotation
			     pd4_pos,
			     parentVInfo,
			     Config.getInt("pd4.copyNumber"), // assign copy number for volume tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #5
//-----------------------------------------------------------------------------
    vector<double> pd5_params;
    Config.getVectorDouble("pd5.halfLengths",pd5_params);

    G4ThreeVector pd5_pos;
    pd5_pos.set(pd4_pos[0]+pd4_params[0]+0.100+pd5_params[0],pd16_ymin+pd5_params[1],pd_z);

    VolumeInfo pd5 (nestBox( "InGaAsPD5",
			     pd5_params,
			     gaas,
			     0,                                 // no rotation
			     pd5_pos,
			     parentVInfo,
			     Config.getInt("pd5.copyNumber"), // for tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
//-----------------------------------------------------------------------------
// photodiode #6
//-----------------------------------------------------------------------------
    vector<double> pd6_params;
    Config.getVectorDouble("pd6.halfLengths",pd6_params);

    G4ThreeVector pd6_pos;
    pd6_pos.set(pd5_pos[0]+pd5_params[0]+0.100+pd6_params[0],pd16_ymin+pd6_params[1],pd_z);

    VolumeInfo pd6 (nestBox( "InGaAsPD6",
			     pd6_params,
			     gaas,
			     0,                                 // no rotation
			     pd6_pos,
			     parentVInfo,
			     Config.getInt("pd6.copyNumber"), // assign copy number for volume tracking purposes
			     pdVisible,
			     G4Color::Red(),
			     pdSolid,
			     forceAuxEdgeVisible,
			     placePV,
			     doSurfaceCheck
			     ));
    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructGen3Setup)
