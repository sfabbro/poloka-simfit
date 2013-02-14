// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMFITPHOT__H
#define SIMFITPHOT__H

#include "simfit.h"
class SimFitPhot {
  
public:
  
  SimFit zeFit;

  SimFitPhot(LightCurveList& Fiducials,bool usegal=true);
  
  void operator() (LightCurve& Lc);
  bool bWriteVignets;
  bool bWriteLC;
  bool bOutputDirectoryFromName;
  
  

};


#endif // SIMFITPHOT__H
