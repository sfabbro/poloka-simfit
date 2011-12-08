// This may look like C code, but it is really -*- C++ -*-
#ifndef VIGNET__H
#define VIGNET__H

#include <dimage.h>
#include <countedref.h>

#include "photstar.h"
#include "fiducial.h"

//
//! \file 
//! \brief A resizeable rectangle, with centered coordinates. 
//
//! Contains a region of an image, with a weight, resizeable and can read from a FITS
//! It can hold local data, weight and associated residual pixels, as well as
//! references to any kind of star and a ReducedImage. If nothing is loaded, it should
//! keep pretty low memory usage, according to the Kernel empty constructor.
//

class Vignet : public Fiducial<Window>, public RefCount {

protected:

  int hx,hy;       // current half size of the vignet (hSizeX and hSizeY are the max sizes).
  void Allocate();

  double exptime; 
  double seeing;
  double mjd;
  
public:

  //! default constructor initialize everything to zero
  Vignet(): hx(0), hy(0), seeing(-1), mjd(-1) {}

  //! associate a ReducedImage
  Vignet(const ReducedImage* Rim) : Fiducial<Window>(Rim), seeing(-1), mjd(-1) {}

  //! extract a squared vignette around a point on a image and a weight of a ReducedImage
  Vignet(const ReducedImage* Rim, const int HMax)
    : Fiducial<Window>(Rim), hx(HMax), hy(HMax), seeing(-1), mjd(-1)   { Resize(hx,hy); }

  //! extract a rectangular vignette around a point on a image and a weight of a ReducedImage
  Vignet(const ReducedImage* Rim, const int HMaxX, const int HMaxY)
    : Fiducial<Window>(Rim), hx(HMaxX), hy(HMaxY), seeing(-1), mjd(-1) { Resize(hx,hy); }

  //! extract a squared vignette around a point on a image and a weight of a ReducedImage
  Vignet(const Fiducial<PhotStar> *Fs, const int HMax)
    : hx(HMax), hy(HMax), seeing(-1), mjd(-1)  { AssignImage(Fs->Image()); Resize(hx,hy); Load(Fs);}

  //! extract a squared vignette around a point on a image and a weight of a ReducedImage
  Vignet(const PhotStar *Star, const ReducedImage* Rim, const int HMax)
    : Fiducial<Window>(Rim), hx(HMax), hy(HMax), seeing(-1), mjd(-1)   { Resize(hx,hy); Load(Star); }

  //! extract a rectangular vignette around a point on a image and a weight of a ReducedImage
  Vignet(const PhotStar *Star, const ReducedImage* Rim, const int HMaxX, const int HMaxY)
    : Fiducial<Window>(Rim), hx(HMaxX), hy(HMaxY), seeing(-1), mjd(-1) { Resize(hx,hy); Load(Star); }

  //! reads a vignet from a vignette FITS file
  Vignet(const string &FitsFileName) : seeing(-1), mjd(-1) {Data.readFits(FitsFileName); }

  // default destructor, copy constructor and assigning operator are OK

  CountedRef<PhotStar> Star;  

  //! the associated data of the fiducial
  Kernel Data;

  //! the associated weight of the image. Follow rules applied to the image itself
  Kernel Weight;

  //! the residual pixels if u need it
  Kernel Resid;
    
  //! initializer from a ReducedImage, loads the calibrated and weight vignets around a star
  virtual void Load(const PhotStar *Star);

  //! resize the Vignet (attempt to modify hx and hy)
  void Resize(const int Hx, const int Hy);

  //! resize the Vignet with a scale factor
  void Resize(const double &ScaleFactor);

  //! shift the center of the vignet
  bool ShiftCenter(const Point& Shift);

  //! function to de-weight high residuals.
  //! inspired from http://nedwww.ipac.caltech.edu/level5/Stetson/frames.html
  void RobustifyWeight(const double& alpha=3., const double& beta=6.);
  
  //! set weights of pixels above or below nsigma*sigma to 0
  //! some kind of robustification
  void KillOutliers(const double& nsigma=5);
  
  //! set to zero resid pixels with null weight
  void ClearResidZeroWeight();

  //! compute and return the chi2 for this vignet
  virtual double Chi2() const;

  //! compute and returns mean residual of the vignet
  double MeanResid() const;
 
   //! compute and returns sigma of residual of the vignet
  double SigmaResid() const;
 
  //! compute and returns the most postivie or negative pixel value of the residual of the vignet
  double MaxPixResid() const;
  
  //! return the number of valid pixels
  int NValidPixels() const;
  
 
  //! insert data into a file
  void writeInImage(const string& FitsFileName) const;

  //! enable "cout << Vignet << endl"
  friend ostream& operator << (ostream & stream, const Vignet& Vig);


  ostream& Info(ostream& stream) const {
    stream << Chi2() << " "
	   << MeanResid() << " "
	   << SigmaResid() << " "
	   << MaxPixResid() << " "
	   << NValidPixels() << " ";
    return stream;
  }

};

//! compute and write the model into a larger FITS image
void writeModel(const Vignet& Vig, const string& FitsFileName);
  
#endif // VIGNET__H
