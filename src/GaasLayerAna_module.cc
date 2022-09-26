///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: GaasLayerAna_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
// $Author: murat $
// $Date: 2014/10/02 17:15:09 $
//
// Contact person:  Pavel Murat
//
// Debug_003: look at the systematics between the StepPointMCs and StrawHitPosition's
// Debug_004: look at various hit-level MC distributions
//
// .fcl file to use: murat/test/trackerMCCheck.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Offline/GeometryService/inc/VirtualDetector.hh"

// #include "Offline/MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "Stntuple/mod/StntupleModule.hh"

// ROOT includes
// #include "TApplication.h"
// #include "TArc.h"
// #include "TArrow.h"
// #include "TCanvas.h"
// #include "TDirectory.h"
// #include "TGraph.h"
#include "TH1F.h"
// #include "TLine.h"
// #include "TBox.h"
// #include "TMarker.h"
// #include "TEllipse.h"
// #include "TText.h"
// #include "TNtuple.h"

// Other includes
// #include "CLHEP/Units/SystemOfUnits.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class GaasLayerAna : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string        fModuleLabel;	             // this module label
    std::string        processName_;
    std::string        fG4ModuleLabel;
    
    std::string        producerName_;
    std::string        fStrawHitMaker;
    std::string        fStrawHitPosMaker;
    std::string        fFlagBgrHitsModuleLabel;
    
    int                fPdgCode;
    int                fGeneratorCode;

    const mu2e::StepPointMCCollection*           fSteps;            //

    struct Hist_t {
      TH1F*   fEDep;                  // radial distance between the StepPointMC and the corresponding hit
    } ;

    Hist_t fHist;


  public:
    explicit GaasLayerAna(fhicl::ParameterSet const& pset);
    virtual ~GaasLayerAna();

    void     getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
    void     Debug_003();   // handles fDr
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     endJob  ();
    virtual bool     filter (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  GaasLayerAna::GaasLayerAna(fhicl::ParameterSet const& pset): 
    StntupleModule            (pset,"GaasLayerAna"),
    fModuleLabel              (pset.get<std::string>("module_label"                 )),
    processName_              (pset.get<std::string>("processName"          ,""     )),
    fG4ModuleLabel            (pset.get<std::string>("g4ModuleLabel"        ,"g4run"))
  {

  }

//-----------------------------------------------------------------------------
  GaasLayerAna::~GaasLayerAna() { 
  }


//-----------------------------------------------------------------------------
  void GaasLayerAna::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }
//-----------------------------------------------------------------------------
  void GaasLayerAna::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

//     art::TFileDirectory tfdir = tfs->mkdir( "CosmicDYB" );
//     _cosmicMultiplicityH = tfdir.make<TH1D>( "MultiplicityH", "Cosmic Multiplicity", 20, -0.5, 19.5);

    fHist.fEDep = tfs->make<TH1F>("edep" ,"Deposited Energy", 1200,0, 6);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void GaasLayerAna::getData(const art::Event* Evt) {
    //    int   rc (0);
    //    const char* oname = "GaasLayerAna::getData";

    art::Handle<StepPointMCCollection> stepsHandle;
    art::Selector sel(art::ProductInstanceNameSelector("stepper") &&
		      art::ProcessNameSelector(processName_) &&
		      art::ModuleLabelSelector(fG4ModuleLabel)  );
    Evt->get(sel, stepsHandle);
    if (stepsHandle.isValid()) fSteps =  (const mu2e::StepPointMCCollection*) &(*stepsHandle);
    else                       fSteps = NULL;
  }

//-----------------------------------------------------------------------------
  void GaasLayerAna::Init(art::Event* Evt) {
  }


  //-----------------------------------------------------------------------------
  bool GaasLayerAna::filter(art::Event& Evt) {
    const char* oname = "GaasLayerAna::filter";

    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

//-----------------------------------------------------------------------------
// get event data and initialize data blocks
//-----------------------------------------------------------------------------
    getData(&Evt);

    Debug_003();
    //    if  (DebugBit(3)) Debug_003();

    return true;
  }

//-----------------------------------------------------------------------------
// plot distribution in radial distance between the StepPointMC and the 
// corresponding hit - fHist.fDr
//-----------------------------------------------------------------------------
  void GaasLayerAna::Debug_003() {
    
    int nsteps;

    const mu2e::StepPointMC        *step;

    nsteps = fSteps->size();
    
    double edep = 0;

    for (int i=0; i<nsteps; i++) {
      step =  &fSteps->at(i);
      if (step->volumeId() == 2) {
	edep += step->totalEDep();
      }
    }

    fHist.fEDep->Fill(edep);
  }
    
}

using mu2e::GaasLayerAna;
DEFINE_ART_MODULE(GaasLayerAna);
