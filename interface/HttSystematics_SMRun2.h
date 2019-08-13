#ifndef SMRun2Legacy_HttSystematics_SMRun2_h
#define SMRun2Legacy_HttSystematics_SMRun2_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 MSSM (with SM categories) analysis systematics
// Implemented in src/HttSystematics_MSSMvsSMRun2.cc
void AddSMRun2Systematics(CombineHarvester& cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, int era);
}

#endif
