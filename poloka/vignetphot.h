// This may look like C code, but it is really -*- C++ -*-
#ifndef VIGNETPHOT__H
#define VIGNETPHOT__H

#include <poloka/vignet.h>
#include <poloka/gausspsf.h>

class VignetPhot : public Vignet {

  DPixel *pim, *pw;
  double wf, sumf, sumwf;

  double p, dy2,dist2;

  double rad;
  double radAperMin2, radAperMax2;

  double radSkyMin2,  radSkyMax2;
  double *skyArray;

  double xp, yp, dpdx, dpdy, sumxr, sumyr, res;
  double wx, wy, sumx2, sumy2, sumxy;


public:

  VignetPhot() : pim(0),  pw(0), rad(0), skyArray(0) {}

  ~VignetPhot() { if (skyArray) delete [] skyArray; if (psf) delete psf; }

  GaussPsf* psf;

  void SetRadius(const double& Rad);

  void SetSkyRadius(const double& RadSkyMin, const double& RadSkyMax);


  void Aperture();

  void SkyAnnulus();

  void Recentroid(const int MaxIter=10, const double& Eps=0.001);

  void Photometry();

  void operator () (const PhotStar *Star);

};



#endif // VIGNETPHOT__H
