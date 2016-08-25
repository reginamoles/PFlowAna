// Duplicate checks includes
#include "EventLoopAlgs/DuplicateChecker.h"

int main( ) {
  bool good = EL::DuplicateChecker::processSummary ("hist-mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_JETM3.e3601_s2576_s2132_r7725_r7676_p2666_tid08619356_00.root", "duplicate_info");
  return 0;
}
