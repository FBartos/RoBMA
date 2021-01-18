#include <module/Module.h>
#include "distributions/DWT1.h"
#include "distributions/DWT2.h"
#include "distributions/DWN1.h"
#include "distributions/DWN2.h"
#include "distributions/DTboost.h"
#include "distributions/DWT1boost.h"
#include "distributions/DWT2boost.h"
#include "distributions/DNboost.h"
#include "distributions/DWN1boost.h"
#include "distributions/DWN2boost.h"
#include "functions/LogWeightedtFun.h"

namespace jags { 
namespace weightedt { // module namespace

  // JAGS module class
  class WTModule : public Module {
    public:
      WTModule();
      ~WTModule();
  };

  // constructor (executed when loading the module)
  WTModule::WTModule() : Module("RoBMA"){
    insert(new DWT1);
    insert(new DWT2);
    insert(new DWN1);
    insert(new DWN2);
    insert(new DTboost);
    insert(new DWT1boost);
    insert(new DWT2boost);
    insert(new DNboost);
    insert(new DWN1boost);
    insert(new DWN2boost);
    //load functions
    insert(new LogWeightedtFun);
  }
  
  // destructor (executed when unloading the module)
  WTModule::~WTModule() {
    std::vector<Function*> const &fvec = functions();
    for (unsigned int i = 0; i < fvec.size(); ++i) {
      delete fvec[i];
    }
    std::vector<Distribution*> const &dvec = distributions();
    for (unsigned int i=0;i<dvec.size();++i) {
      delete dvec[i];
    }
  }

}
}

jags::weightedt::WTModule _weightedt_module;

