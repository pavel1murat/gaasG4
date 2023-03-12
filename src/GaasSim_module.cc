///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: GaasSim_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
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
#include "art/Framework/Core/EDProducer.h"
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

#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"

#include "gaasG4/inc/GaasSimAlg.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class GaasSim : public art::EDProducer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   spmcCollTag            {Name("spmcCollTag"       )    , Comment("StepPointMC collection tag") };
      fhicl::Atom<int>             debugLevel             {Name("debugLevel"        )    , Comment("debug level"                ) };
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"         )    , Comment("diag level"                 ) };
    };

  public:
    enum {
      kNEventHistSets = 100 
    };

    art::InputTag                       _spmcCollTag;
    const mu2e::StepPointMCCollection*  _spmcColl;            //
    int                                 _debugLevel;
    
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


    GaasSimAlg*  _simAlg;
    int          _nphotons;
    TRandom3     _rn;

  public:
    explicit         GaasSim(const art::EDProducer::Table<Config>& config);    
    virtual         ~GaasSim();

    void             getData(const art::Event& AnEvent);

    void             bookEventHistograms(EventHist_t*    Hist, art::TFileDirectory* Dir);
    void             bookHistograms();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob()                 override;
    virtual void     endJob  ()                 override;
    virtual void     produce (art::Event& Evt)  override;
  };


//-----------------------------------------------------------------------------
  GaasSim::GaasSim(const art::EDProducer::Table<Config>& config) : 
    EDProducer                (config),
    _spmcCollTag              (config().spmcCollTag())
  {
    _simAlg = new GaasSimAlg();
  }

//-----------------------------------------------------------------------------
  GaasSim::~GaasSim() { 
  }


//-----------------------------------------------------------------------------
  void GaasSim::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }


//-----------------------------------------------------------------------------
  void GaasSim::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {

    art::ServiceHandle<art::TFileService> tfs;

    Hist->fNSteps   = Dir->make<TH1F>("nsteps"   ,"N(steps), GaAs+PD"      , 100 ,0, 100);
    Hist->fEDepTot  = Dir->make<TH1F>("edep_tot" ,"Deposited Energy, total", 6000,0, 6);
    Hist->fEDepGaas = Dir->make<TH1F>("edep_gaas","Deposited Energy, GaAs" , 6000,0, 6);
    Hist->fEDepPD   = Dir->make<TH1F>("edep_pd"  ,"Deposited Energy, PD"   , 6000,0, 6);
    Hist->fYVsX     = Dir->make<TH2F>("y_vs_x"   ,"Y vs X, all"            , 200 ,-0.50,0.50,200,-0.50,0.50);
    Hist->fDeVsZ    = Dir->make<TH1F>("de_vs_z"  ,"dE vs Z"                , 20  , 0,20e-3);
  }

//-----------------------------------------------------------------------------
  void GaasSim::bookHistograms() {
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
  void GaasSim::beginJob() {

  // init sensor dimensions : sensor used for absorption measurements : 8.0 x 0.7 x 0.025

  double sensor[3] = { 0.400, 0.0350, 0.00125 };

  // init detectors

  GaasSimAlg::Detector_t det;

  det.fDx      = 0.0020;	 	//  50 um
  det.fDy      = 0.0025;		// 300 um
  det.fDz      = 0.00125;                // just to make the diode visible, thickness is not used so far

  det.fSide    = 0;			// Z-side,top
  
  det.fOffsetX = -0.4000;
  det.fOffsetY =  0.0000;
  det.fOffsetZ =  0.0000;
  det.fType    = 1;       // fiber with the air gap

  int ndet = 1;
  _simAlg->InitGeometry(sensor, &det, ndet);


  _simAlg->fPosMode = 1;

    bookHistograms();
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void GaasSim::getData(const art::Event& AnEvent) {

    art::Handle<StepPointMCCollection> spmccH;

    AnEvent.getByLabel(_spmcCollTag,spmccH);

    if (spmccH.isValid()) _spmcColl =  spmccH.product();
    else                  _spmcColl = NULL;
  }

  //-----------------------------------------------------------------------------
  void GaasSim::produce(art::Event& AnEvent) {
    //const char* oname = "GaasSim::analyzer";
    // printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    getData(AnEvent);

// steppoints are read in, read and convert into number of photons

    float edep_tot = 0;

    if (_simAlg->fPosMode == 0) { 
//-----------------------------------------------------------------------------
// fixed position 
//-----------------------------------------------------------------------------
      // float fX0 = _simAlg->fDiodePos[0]+0.1; // Dist;
      // float fY0 = 0.; // 2*(rn[1]-0.5)*cPos->GetDY();
      // float fZ0 = 0.; // 2*(rn[2]-0.5)*cPos->GetDZ();
    }
    else {
//-----------------------------------------------------------------------------
// read steppoints
//-----------------------------------------------------------------------------
      int nsteps = _spmcColl->size();
    
      for (int i=0; i<nsteps; i++) {
        const mu2e::StepPointMC* step =  &_spmcColl->at(i);
      
        float edep = step->totalEDep();

        edep_tot += edep;

        // float x    = step->position().z();
        // float y    = step->position().z();
        // float z    = step->position().z();
      }
    }

//-----------------------------------------------------------------------------
// need to determine from the energy deposition
//-----------------------------------------------------------------------------
    int fNPhMean = 1*edep_tot;

    if (fNPhMean > 0) _nphotons = _rn.Poisson(fNPhMean);
    else              _nphotons = -fNPhMean;
   
    if (_debugLevel > 0) {
      printf(" >>>> %s : nphotons:  %12i\n",__func__, _nphotons);
    }
   
    _simAlg->simulate(_nphotons);
  }
}

using mu2e::GaasSim;
DEFINE_ART_MODULE(GaasSim);
