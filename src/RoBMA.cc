#include <module/Module.h>
#include "distributions/DWT1.h"
#include "distributions/DWT2.h"
#include "distributions/DWN1.h"
#include "distributions/DWN2.h"
#include "distributions/DWN.h"
#include "distributions/DWWN1.h"
#include "distributions/DWWN2.h"
//#include "distributions/DMN.h"
#include "distributions/DWMN1.h"
#include "distributions/DWMN2.h"
#include "distributions/DMNv.h"
#include "distributions/DWMN1v.h"
#include "distributions/DWMN2v.h"
#include "distributions/DWB.h"

#include "transformations/d.h"
#include "transformations/r.h"
#include "transformations/z.h"
#include "transformations/logOR.h"

#include "functions/mnorm.h"
#include "functions/wmnorm.h"

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
      insert(new DWT1);
      insert(new DWT2);
      insert(new DWN1);
      insert(new DWN2);
      insert(new DWN);
      insert(new DWWN1);
      insert(new DWWN2);
      //insert(new DMN);
      insert(new DWMN1);
      insert(new DWMN2);
      insert(new DMNv);
      insert(new DWMN1v);
      insert(new DWMN2v);
      insert(new DWB);

      //effect sizes transformations
      insert(new d2z);
      insert(new d2r);
      insert(new d2logOR);
      insert(new r2d);
      insert(new r2z);
      insert(new r2logOR);
      insert(new z2r);
      insert(new z2d);
      insert(new z2logOR);
      insert(new logOR2d);
      insert(new logOR2z);
      insert(new logOR2r);

      //standard errors transformations
      insert(new se_d2se_z);
      insert(new se_d2se_r);
      insert(new se_d2se_logOR);
      insert(new se_r2se_d);
      insert(new se_r2se_z);
      insert(new se_r2se_logOR);
      insert(new se_z2se_r);
      insert(new se_z2se_d);
      insert(new se_z2se_logOR);
      insert(new se_logOR2se_d);
      insert(new se_logOR2se_z);
      insert(new se_logOR2se_r);

      //prior scaling functions (aproximate linear transformations)
      insert(new scale_d2z);
      insert(new scale_d2logOR);
      insert(new scale_z2d);
      insert(new scale_z2logOR);
      insert(new scale_logOR2d);
      insert(new scale_logOR2z);
      insert(new scale_r2d);
      insert(new scale_r2z);
      insert(new scale_r2logOR);
      insert(new scale_d2r);
      insert(new scale_z2r);
      insert(new scale_logOR2r);

      // likelihood functions
      insert(new mnorm_lpdf);
      insert(new wmnorm_1s_lpdf);
      insert(new wmnorm_2s_lpdf);
      insert(new mnorm_v_lpdf);
      insert(new wmnorm_1s_v_lpdf);
      insert(new wmnorm_2s_v_lpdf);
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

