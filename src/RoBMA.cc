#include <module/Module.h>
#include <version.h>

static_assert(JAGS_MAJOR >= 5, "RoBMA requires JAGS 5.0.0 or later");

#include "distributions/DWN.h"
#include "distributions/DWB.h"
#include "distributions/DWP.h"
#include "distributions/DKNOWNVMNORM.h"

#include "distributions/DSELNORMKERNEL.h"
#include "distributions/DSELNORMSTEP.h"
#include "distributions/DSELNORMSTEPSWITCH.h"
#include "distributions/DSELNORMMVSTEP.h"
#include "distributions/DSELNORMCLUSTERSTEP.h"
#include "distributions/DSELNORMFACTORSTEP.h"
#include "distributions/DSELNORMSAMPLINGCONDITIONED.h"
#include "samplers/CoarseCorrectedSlice.h"
#include "samplers/CoarseCorrectedBlock.h"
#include "samplers/IndependentBiasJoint.h"
#include "samplers/GaussianJoint.h"
#include "samplers/JointCluster.h"

namespace jags {
  namespace RoBMA { // module namespace

    // JAGS module class
    class RoBMAModule : public Module {
      IndependentBiasJointFactory independent_bias_joint;
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
      insert(new DSELNORMMVSTEP);
      insert(new DSELNORMCLUSTERSTEP);
      insert(new DSELNORMFACTORSTEP);
      insert(new DSELNORMSAMPLINGCONDITIONED);
      // These process-lifetime factories start disabled (the shared corrected
      // sampler configuration) and own no fitted model. JAGS asks the factories
      // of a module loaded later first, so RoBMA's (loaded after BayesTools')
      // claim their blocks before BayesTools' factor slice factory.
      insert(coarse_corrected_slice_factory());
      insert(coarse_corrected_block_factory());
      insert(new DIndependentBiasJointMarker);
      insert(&independent_bias_joint);
      insert(new GaussianJointMarker);
      insert(gaussian_joint_factory());
      insert(new ClusterMarker(false));
      insert(new ClusterMarker(true));
      insert(joint_cluster_factory());
    }

    // destructor (executed when unloading the module)
    RoBMAModule::~RoBMAModule() {
      // Normal process shutdown does not invoke the R namespace unload hook.
      unload();

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

