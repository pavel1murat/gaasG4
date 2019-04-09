///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "art/Utilities/ToolMacros.h"

// Mu2e includes.

// #include "Mu2eG4/inc/ConstructMaterials.hh"
// #include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
// #include "Mu2eG4/inc/findMaterialOrThrow.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Helper/inc/G4Helper.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

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

    G4Material* gaas = G4Material::GetMaterial("GaAs",false);

    if (gaas == NULL) {
      gaas = new G4Material("GaAs", 5.3*CLHEP::g/CLHEP::cm3, 2 );

      G4Element* ga = nist->FindOrBuildElement("Ga",true);
      gaas->AddElement( ga, 1);

      G4Element* as = nist->FindOrBuildElement("As",true);
      gaas->AddElement( as, 1);
    }

    // Stainless Steel (Medical Physics, Vol 25, No 10, Oct 1998) based on brachytherapy example

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
//-----------------------------------------------------------------------------
    G4bool tubeVisible        = Config.getBool("tube.visible",true);
    G4bool tubeSolid          = Config.getBool("tube.solid",true);

    TubsParams tubeParams( Config.getDouble("tube.rIn"),
                           Config.getDouble("tube.rOut"),
                           Config.getDouble("tube.halfLength"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    //    MaterialFinder materialFinder(Config);

    //    G4Material* tubeMaterial = materialFinder.get("tube.material");

    const G4ThreeVector tubeCenterInWorld(Config.getHep3Vector("tube.centerInWorld"));

    // CLHEP::HepRotationZ rotZ(_config.getDouble("tube.phiRotZ")*CLHEP::degree);
    // CLHEP::HepRotationY rotY(M_PI);
    // G4RotationMatrix* rotation  = (sgn < 0 )?
    //   reg.add(G4RotationMatrix(rotZ)) :
    //   reg.add(G4RotationMatrix(rotZ*rotY));

    VolumeInfo tubeVInfo(nestTubs( "SourceHolder",
                                   tubeParams,
                                   poly, // material
                                   0, // rotation,
                                   tubeCenterInWorld,
                                   parentVInfo,
                                   Config.getInt("tube.copyNumber"), // non-0 for tracking purposes
                                   tubeVisible,
                                   G4Color::Blue(),
                                   tubeSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// source itself
//-----------------------------------------------------------------------------
//    G4Material* src1Material = materialFinder.get("src1.material");

    TubsParams src1Params( Config.getDouble("src1.rIn"),
                           Config.getDouble("src1.rOut"),
                           Config.getDouble("src1.halfLength"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    const G4ThreeVector src1CenterInWorld(Config.getHep3Vector("src1.centerInWorld"));

    // CLHEP::HepRotationZ rotZ(_config.getDouble("tube.phiRotZ")*CLHEP::degree);
    // CLHEP::HepRotationY rotY(M_PI);
    // G4RotationMatrix* rotation  = (sgn < 0 )?
    //   reg.add(G4RotationMatrix(rotZ)) :
    //   reg.add(G4RotationMatrix(rotZ*rotY));

    VolumeInfo src1VInfo(nestTubs( "Src1",
                                   src1Params,
                                   steel,
                                   0, // rotation,
                                   src1CenterInWorld,
                                   parentVInfo,
                                   Config.getInt("src1.copyNumber"), // non-0 for tracking purposes
                                   tubeVisible,
                                   G4Color::Blue(),
                                   tubeSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// source itself
//-----------------------------------------------------------------------------
//    G4Material* src2Material = materialFinder.get("src2.material");

    TubsParams src2Params( Config.getDouble("src2.rIn"),
                           Config.getDouble("src2.rOut"),
                           Config.getDouble("src2.halfLength"),
                           0.*CLHEP::degree,
                           360.*CLHEP::degree );

    const G4ThreeVector src2CenterInWorld(Config.getHep3Vector("src2.centerInWorld"));

    VolumeInfo src2VInfo(nestTubs( "Src2",
                                   src2Params,
                                   steel,
                                   0, // rotation,
                                   src2CenterInWorld,
                                   parentVInfo,
                                   Config.getInt("src2.copyNumber"), // non-0 for tracking purposes
                                   tubeVisible,
                                   G4Color::Blue(),
                                   tubeSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
//-----------------------------------------------------------------------------
// construct and position sensor
//-----------------------------------------------------------------------------
    G4bool boxVisible        = Config.getBool("box.visible",true);
    G4bool boxSolid          = Config.getBool("box.solid",true);

    vector<double> boxParams;
    Config.getVectorDouble( "box.halfLengths", boxParams);

    const G4ThreeVector boxCenterInWorld(Config.getHep3Vector("box.centerInWorld"));

    VolumeInfo box(nestBox( "GaasSensor",
			    boxParams,
			    gaas,
			    0,                               // no rotation
			    boxCenterInWorld,
			    parentVInfo,
			    Config.getInt("box.copyNumber"), // assign copy nuber for volume tracking purposes
			    boxVisible,
			    orange,
			    boxSolid,
			    forceAuxEdgeVisible,
			    placePV,
			    doSurfaceCheck
			    ));
//-----------------------------------------------------------------------------
// construct and position the photodiode
//-----------------------------------------------------------------------------
    G4bool pdVisible        = Config.getBool("pd.visible",true);
    G4bool pdSolid          = Config.getBool("pd.solid",true);

    // G4Material* pdMaterial = materialFinder.get("pd.material");

    vector<double> pdParams;
    Config.getVectorDouble("pd.halfLengths",pdParams);

    G4ThreeVector pdPos;
    pdPos.set(boxCenterInWorld.x()+boxParams[0]-0.5,
	      boxCenterInWorld.y(),
	      boxCenterInWorld.z()+boxParams[2]+pdParams[2]);

    VolumeInfo pd (nestBox( "InGaAsPD",
			    pdParams,
			    gaas,
			    0,                                 // no rotation
			    pdPos,
			    parentVInfo,
			    Config.getInt("pd.copyNumber"), // assign copy nuber for volume tracking purposes
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
