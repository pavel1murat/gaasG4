// Pass events with at least one hit satisfying a min momentum cut.
//
#include <string>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  //================================================================
  class SpmcFilter : public art::EDFilter {
    art::InputTag _spmcCollTag;
    double        _minEDep;
                                        // statistics counters
    unsigned      _nInputEvents;
    unsigned      _nPassedEvents;
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> spmcCollTag { Name("spmcCollTag" ), Comment("SPMC coll tag")  };
      fhicl::Atom<double>        minEDep     { Name("minEDep"     ), Comment("min edep"     )  };
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit SpmcFilter(const Parameters& conf);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  SpmcFilter::SpmcFilter(const Parameters& conf) : art::EDFilter{conf}
    , _spmcCollTag  (conf().spmcCollTag())
    , _minEDep      (conf().minEDep    ())
    , _nInputEvents (0)
    , _nPassedEvents(0)
    { }

  //================================================================
  bool SpmcFilter::filter(art::Event& event) {
    bool passed = false;

    auto ih = event.getValidHandle<StepPointMCCollection>(_spmcCollTag);

    double edep = 0;                   // hit is 
    for(const auto& spmc : *ih) {
      if (spmc.volumeId() > 999) { 
        edep += spmc.totalEDep();
      }
    }

    if (edep > _minEDep) passed = true;

    ++_nInputEvents;
    if (passed)  ++_nPassedEvents; 
    return passed;
  }

  //================================================================
  void SpmcFilter::endJob() {
    mf::LogInfo("Summary")
      <<"SpmcFilter_module: passed "
      << _nPassedEvents << " / " << _nInputEvents << " events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SpmcFilter);
