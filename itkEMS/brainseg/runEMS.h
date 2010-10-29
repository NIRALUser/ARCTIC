
// Wrapper function that puts everything together (filtering, registration,
// segmentation)

// prastawa@cs.unc.edu 11/2003

#ifndef _runEMS_h
#define _runEMS_h

#include "EMSParameters.h"
#include <vector>
#include <string>

extern void runEMSFull(EMSParameters * emsp, const bool debugflag, const bool writemoreflag, std::string templateVolume, std::vector<std::string> priorsList, std::vector<double> priorsWeightList );
extern void runEMS(EMSParameters * emsp, const bool debugflag, const bool writemoreflag);

#endif
