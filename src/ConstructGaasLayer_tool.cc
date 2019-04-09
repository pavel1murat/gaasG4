//
// Places the test setup inside the box
//

#include "art/Utilities/ToolMacros.h"

// Mu2e includes.

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

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

  class ConstructGaasLayer: public InitEnvToolBase {
  public:
    ConstructGaasLayer(const fhicl::ParameterSet& PSet);
    ~ConstructGaasLayer();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructGaasLayer::ConstructGaasLayer(const fhicl::ParameterSet& PSet) {
    _name = "GaasLayer";
  }

//-----------------------------------------------------------------------------
  ConstructGaasLayer::~ConstructGaasLayer() {
  }

//-----------------------------------------------------------------------------
  int ConstructGaasLayer::construct(VolumeInfo const& parentVInfo, SimpleConfig const& _config) {

    // 2016-11-06 P.Murat: add GaAs

    G4Material* mat = G4Material::GetMaterial("GaAs",false);
    if (mat == NULL) {
      G4Material* GaAs = new G4Material("GaAs", 5.3*CLHEP::g/CLHEP::cm3, 2 );
      G4NistManager* nistMan = G4NistManager::Instance();

      G4Element* ga = nistMan->FindOrBuildElement("Ga",true);
      GaAs->AddElement( ga, 1);

      G4Element* as = nistMan->FindOrBuildElement("As",true);
      GaAs->AddElement( as, 1);
    }

    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract box information from the config file.

    G4bool boxVisible        = _config.getBool("box.visible",true);
    G4bool boxSolid          = _config.getBool("box.solid",true);

    vector<double> boxParams;
    _config.getVectorDouble( "box.halfLengths", boxParams);

    MaterialFinder materialFinder(_config);
    G4Material* boxMaterial = materialFinder.get("box.wallMaterialName");

    const G4ThreeVector boxCenterInWorld(_config.getHep3Vector("box.centerInWorld"));

    G4Colour  orange  (.75, .55, .0);

    VolumeInfo box(nestBox( "GaasLayer",
			    boxParams,
			    boxMaterial,
			    0,                                // no rotation
			    boxCenterInWorld,
			    parentVInfo,
			    _config.getInt("box.copyNumber",2), 
			    // we assign a non 0 copy nuber for
			    // volume tracking purposes
			    boxVisible,
			    orange,
			    boxSolid,
			    forceAuxEdgeVisible,
			    placePV,
			    doSurfaceCheck
			    ));
    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructGaasLayer)
