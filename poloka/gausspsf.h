// This may look like C code, but it is really -*- C++ -*-
#ifndef GAUSSPSF__H
#define GAUSSPSF__H

#include <poloka/reducedimage.h>

class GaussPsf {

  double beta, alphax, alphay, rxy, norm;

public:

  //! empty constructor make sure array are not pointing to anything
  GaussPsf() {}

  //! constructor builds q gaussian psf from seeing properties
  GaussPsf(const ReducedImage &Rim);
    
  //! compute the psf value in a pixel (i,j) from a point (Xc,Yc) and its derivatives
  double Value(const int i, const int j, const double &Xc, const double &Yc, 
	       double &DpDx, double &DpDy) const;

};

#endif // GAUSSPSF__H
