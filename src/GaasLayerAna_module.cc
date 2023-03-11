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

// ROOT includes
// #include "TApplication.h"
// #include "TArc.h"
// #include "TArrow.h"
// #include "TCanvas.h"
// #include "TDirectory.h"
// #include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
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

  class GaasLayerAna : public art::EDAnalyzer {
  private:
    std::string        _spmcCollTag;
    
    const mu2e::StepPointMCCollection*           _spmcColl;            //
    
    struct Hist_t {
      TH1F*   fNSteps;
      TH1F*   fEDepTot;
      TH1F*   fEDepGaas;
      TH1F*   fEDepPD;
      TH2F*   fYVsX;
    } ;

    Hist_t fHist;

  public:
    explicit GaasLayerAna(fhicl::ParameterSet const& pset);
    virtual ~GaasLayerAna();

    void     getData(const art::Event& AnEvent);
    void     Init   (const art::Event& AnEvent);
    void     Debug_003();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob()                       override;
    virtual void     endJob  ()                       override;
    virtual void     analyze (const art::Event& Evt)  override;
  };


//-----------------------------------------------------------------------------
  GaasLayerAna::GaasLayerAna(fhicl::ParameterSet const& pset): 
    EDAnalyzer                (pset),
    _spmcCollTag              (pset.get<std::string>("spmcCollTag"))
  {}

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

    fHist.fNSteps   = tfs->make<TH1F>("nsteps"   ,"N(steps), GaAs+PD"      , 100 ,0, 100);
    fHist.fEDepTot  = tfs->make<TH1F>("edep_tot" ,"Deposited Energy, total", 6000,0, 6);
    fHist.fEDepGaas = tfs->make<TH1F>("edep_gaas","Deposited Energy, GaAs" , 6000,0, 6);
    fHist.fEDepPD   = tfs->make<TH1F>("edep_pd"  ,"Deposited Energy, PD"   , 6000,0, 6);
    fHist.fYVsX     = tfs->make<TH2F>("y_vs_x"   ,"Y vs X, all"            , 200 ,-0.50,0.50,200,-0.50,0.50);
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void GaasLayerAna::getData(const art::Event& AnEvent) {

    art::Handle<StepPointMCCollection> spmccH;

    /* bool ok = */ AnEvent.getByLabel("g4run:stepper",spmccH);

    if (spmccH.isValid()) _spmcColl =  spmccH.product();
    else                  _spmcColl = NULL;
  }

//-----------------------------------------------------------------------------
  void GaasLayerAna::Init(const art::Event& Evt) {
  }


  //-----------------------------------------------------------------------------
  void GaasLayerAna::analyze(const art::Event& Evt) {
    //const char* oname = "GaasLayerAna::analyzer";

    // printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

    getData(Evt);

    Debug_003();
    //    if  (DebugBit(3)) Debug_003();

    //    return true;
  }

//-----------------------------------------------------------------------------
// plot distribution in radial distance between the StepPointMC and the 
// corresponding hit - fHist.fDr
//-----------------------------------------------------------------------------
  void GaasLayerAna::Debug_003() {
    
    int nsteps;

    nsteps = _spmcColl->size();
    
    double edep_tot(0), edep_gaas(0), edep_pd(0);

    for (int i=0; i<nsteps; i++) {
      const mu2e::StepPointMC* step =  &_spmcColl->at(i);

      int vol_id = step->volumeId();

      if ((vol_id > 999) and (vol_id < 2000)) {
        edep_pd += step->totalEDep();
      }
      else if (step->volumeId() == 2000) {
	edep_gaas += step->totalEDep();
      }

      fHist.fYVsX->Fill(step->position().x(),step->position().y());
    }

    edep_tot = edep_gaas+edep_pd;

    fHist.fNSteps->Fill(nsteps);
    fHist.fEDepGaas->Fill(edep_gaas);
    fHist.fEDepPD->Fill(edep_pd);
    fHist.fEDepTot->Fill(edep_tot);

  }
    
}

using mu2e::GaasLayerAna;
DEFINE_ART_MODULE(GaasLayerAna);
