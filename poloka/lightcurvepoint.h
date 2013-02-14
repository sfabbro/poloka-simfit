// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEPOINT__H
#define LIGHTCURVEPOINT__H

#include <iostream>
#include <iomanip>
#include <cmath>

#include <poloka/countedref.h>

//! lightcurve result
class LightCurvePoint: public RefCount {
public:
  LightCurvePoint();
  double modifiedjulianday;
  double flux;
  double eflux;
  double mag;
  double emag_minus;
  double emag_plus;
  double zeropoint;
  void computemag(double zp) {
    zeropoint = zp;
    if(flux<1.e-12) {
      mag=99;
      emag_minus=0;
      emag_plus=0;
    }else{
      mag=-2.5*log10(flux)+zeropoint;
      if((1.-eflux/flux)<0)
	emag_minus=-99;
      else
	emag_minus=2.5*log10(1.-eflux/flux);
      emag_plus=2.5*log10(1.+eflux/flux);
    }
  }

  friend std::ostream& operator << (std::ostream &Stream, const LightCurvePoint &lcp) { 
    Stream << std::setiosflags(std::ios::fixed);
    //Stream << std::setw(14) << std::setprecision(2) << lcp.julianday
    Stream << std::setw(14) << std::setprecision(2) << lcp.modifiedjulianday
	   << std::setw(15) << std::setprecision(3) << lcp.flux
	   << std::setw(15) << std::setprecision(3) << lcp.eflux
      //   << std::setw(15) << std::setprecision(3) << lcp.mag
      //   << std::setw(15) << std::setprecision(3) << lcp.emag_minus
      //   << std::setw(15) << std::setprecision(3) << lcp.emag_plus
	   << std::setw(15) << std::setprecision(3) << lcp.zeropoint;
    return Stream;
  }

};


#endif // LIGHTCURVEPOINT__H
