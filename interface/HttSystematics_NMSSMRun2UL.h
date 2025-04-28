#ifndef SMRun2Legacy_HttSystematics_NMSSMRun2UL_h
#define SMRun2Legacy_HttSystematics_NMSSMRun2UL_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 NMSSM (with SM categories) analysis systematics
// Implemented in src/HttSystematics_MSSMvsSMRun2.cc
   void AddRun2Systematics(CombineHarvester& cb, bool jetfakes, bool embedding, int era);
}

#endif
