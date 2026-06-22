#include <module/Module.h>
#include "distributions/DWN.h"
#include "distributions/DWB.h"
#include "distributions/DWP.h"
#include "distributions/DKNOWNVMNORM.h"

#include "distributions/DSELNORMKERNEL.h"
#include "distributions/DSELNORMSTEP.h"
#include "distributions/DSELNORMSTEPSWITCH.h"

namespace jags {
  namespace RoBMA { // module namespace

    // JAGS module class
    class RoBMAModule : public Module {
      public:
        RoBMAModule();
        ~RoBMAModule();
    };

    // constructor (executed when loading the module)
    RoBMAModule::RoBMAModule() : Module("RoBMA"){

      // distributions
      insert(new DWN);
      insert(new DWB);
      insert(new DWP);
      insert(new DKNOWNVMNORM);

      // mixture distributions
      insert(new DSELNORMSTEP);
      insert(new DSELNORMSTEPSWITCH);
      insert(new DSELNORMKERNEL);
    }

    // destructor (executed when unloading the module)
    RoBMAModule::~RoBMAModule() {
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

jags::RoBMA::RoBMAModule _RoBMA_module;

