// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVE__H
#define LIGHTCURVE__H

#include "refstar.h"

//! 
//!  \file lightcurve.h
//!  \brief A set of Fiducial taken at different exposures.
//!
 
//! The same Fiducial monitored in many images in a list
class LightCurve : public list<CountedRef<Fiducial<PhotStar> > > {
public:
  
  //! chi2 of the fit
  double chi2; 

  //! number of degrees of freedom
  int ndf; 
  
  //! chi2 per degree of freedom
  double chi2ndf() const {if(ndf==0) return 0.; return chi2/ndf;};
  
  //! empty constructor does strictly nothing
  LightCurve() {ndf=0;chi2=0.;};

  //! load the Ref with a RefStar, does not contain any measurements
  LightCurve(const RefStar *Star) : Ref(Star) {ndf=0;chi2=0.;};

  // default destructor, copy constructor and assigning operator are OK  

  //! a pointer to the reference star
  CountedRef<RefStar> Ref;
   
  //! push back a star to the LightCurve if image does not exist already
  void push_back(const ReducedImage* Rim, const PhotStar *Star);

  //! push back the reference star with an image if that one does not exist already
  void push_back(const ReducedImage* Rim);

  //! normal push_back of the list reference star with an image if that one does not exist already
  void push_back(const Fiducial<PhotStar> *Star)
  { list< CountedRef<Fiducial<PhotStar> > > ::push_back(Star); }

  //! returns the vector of all julian dates, fluxes, covariance matrix in a C style array
  void ComputeMatVec(double *JulianDates, double *Fluxes, double *Covariance) const;

  //! writes on a stream a simple output with format JD FLUX DFLUX IMAGE with a header.
  void write_lc2fit(ostream &Stream = cout) const;

  //! writes in an xml file the lightcurve result (with LightCurvePoint)
  //  void write_xml(const string &filename) const;
  
  //! header for extended write  
  ostream& write_header(ostream &Stream) const;

  //! compute ELIXIR zero point
  double computeElixirZeroPoint() const;

  //! allows extended write of the full Fiducial<PhotStar>
  friend ostream& operator << (ostream& Stream, const LightCurve& Lc);

  // stuff for simfit
  double totflux, vartotflux;
  double galflux, vargalflux;
  double totsky , vartotsky;
  double resmean, resmed, resrms, resadev;

};


//! a class that contains most information to build light curves of objects. 
//! It is basically a list of LightCurve so it has objects information as well as images
//! The Images and Objects are offered for conviency. The LightCurve should not 
//! copy the data in Images and Objects thanks to the CountedRef nature of the LightCurve and Fiducial
class LightCurveList : public list<LightCurve> {
public:
  
  //! empty constructor does strictly nothing
  LightCurveList() {}

  //! fill Images, Objects and the list of LightCurve from a light file stream
  LightCurveList(istream& LcFileStream);

  // default destructor, copy constructor and assigning operator are OK  

  //! the images to reduce
  ReducedImageList Images;

  //! the objects to monitor
  RefStarList Objects;

  //! a pointer to the photometric reference used in the LightCurve
  ReducedImageRef RefImage;

  //! enables to write in a file and to dump on screen
  friend ostream& operator << (ostream& Stream, const LightCurveList& MyList);

  //! dumps some useful info from the fit
  void write(const string& filename) const;
};


#endif // LIGHTCURVE__H
