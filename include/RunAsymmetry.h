#ifndef RUNASYMMETRY_H
#define RUNASYMMETRY_H

#include "TString.h"
#include "Rtypes.h"

// declaration header
Float_t findCentrality( Float_t HFSum );
void runAsymmetry( TString input, TString output, TString modeFlag = "triggered", Long64_t maxEvents = -1 );

#endif
