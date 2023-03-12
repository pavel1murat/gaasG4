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
  public:
    enum {
      kNEventHistSets = 100 
    };

  private:
    std::string        _spmcCollTag;
    
    const mu2e::StepPointMCCollection*           _spmcColl;            //
    
    struct EventHist_t {
      TH1F*   fNSteps;
      TH1F*   fEDepTot;
      TH1F*   fEDepGaas;
      TH1F*   fEDepPD;
      TH2F*   fYVsX;
      TH1F*   fDeVsZ;
    } ;

    struct Hist_t {
      EventHist_t* fEvent[kNEventHistSets];
    } _hist;

  public:
    explicit GaasLayerAna(fhicl::ParameterSet const& pset);
    virtual ~GaasLayerAna();

    void     getData(const art::Event& AnEvent);
    void     Init   (const art::Event& AnEvent);
    void     Debug_003();

    void     bookEventHistograms(EventHist_t*    Hist, art::TFileDirectory* Dir);
    void     bookHistograms();
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
  void GaasLayerAna::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {

    art::ServiceHandle<art::TFileService> tfs;

    Hist->fNSteps   = Dir->make<TH1F>("nsteps"   ,"N(steps), GaAs+PD"      , 100 ,0, 100);
    Hist->fEDepTot  = Dir->make<TH1F>("edep_tot" ,"Deposited Energy, total", 6000,0, 6);
    Hist->fEDepGaas = Dir->make<TH1F>("edep_gaas","Deposited Energy, GaAs" , 6000,0, 6);
    Hist->fEDepPD   = Dir->make<TH1F>("edep_pd"  ,"Deposited Energy, PD"   , 6000,0, 6);
    Hist->fYVsX     = Dir->make<TH2F>("y_vs_x"   ,"Y vs X, all"            , 200 ,-0.50,0.50,200,-0.50,0.50);
    Hist->fDeVsZ    = Dir->make<TH1F>("de_vs_z"  ,"dE vs Z"                , 20  , 0,20e-3);
  }

//-----------------------------------------------------------------------------
  void GaasLayerAna::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    char   folder_name[200];

    TH1::AddDirectory(0);
    
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events
    book_event_histset[ 1] = 1;                // events with the total edep 

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],&tfdir);
      }
    }
  }

//-----------------------------------------------------------------------------
  void GaasLayerAna::beginJob() {
    bookHistograms();
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

      float edep = step->totalEDep();
      float z    = step->position().z();

      int vol_id = step->volumeId();

      if ((vol_id > 999) and (vol_id < 2000)) {
        edep_pd   += edep;
      }
      else if (step->volumeId() == 2000) {
//-----------------------------------------------------------------------------
// GaAs sensor
//-----------------------------------------------------------------------------
        _hist.fEvent[0]->fDeVsZ->Fill(-z,edep);
        if (edep > 0) {
          _hist.fEvent[1]->fDeVsZ->Fill(-z,edep);
          edep_gaas += edep;
        }
      }

      _hist.fEvent[0]->fYVsX->Fill(step->position().x(),step->position().y());
      if (edep > 0) {
        _hist.fEvent[1]->fYVsX->Fill(step->position().x(),step->position().y());
      }
    }

    edep_tot = edep_gaas+edep_pd;

    _hist.fEvent[0]->fNSteps->Fill(nsteps);
    _hist.fEvent[0]->fEDepGaas->Fill(edep_gaas);
    _hist.fEvent[0]->fEDepPD->Fill(edep_pd);
    _hist.fEvent[0]->fEDepTot->Fill(edep_tot);

    if (edep_tot > 0) {
      _hist.fEvent[1]->fNSteps->Fill(nsteps);
      _hist.fEvent[1]->fEDepGaas->Fill(edep_gaas);
      _hist.fEvent[1]->fEDepPD->Fill(edep_pd);
      _hist.fEvent[1]->fEDepTot->Fill(edep_tot);
    }
  }
    
}

using mu2e::GaasLayerAna;
DEFINE_ART_MODULE(GaasLayerAna);
