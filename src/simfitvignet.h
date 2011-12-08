// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMFITVIGNET__H
#define SIMFITVIGNET__H

#include <countedref.h>
//#include <daophotpsf.h>
#include <reducedimage.h>
#include "photstar.h"
#include "vignet.h"
#include "imagepsf.h"
#include "kernelfit.h"

//
//! \file simfitvignet.h
//! \brief a Vignet with its associated tabulated kernels, psfs, residuals.
//!  Primarly to use for SimultaneousFit.
//

// uncomment this to use one daophot psf per image, kernels are still used for the galaxy and for the photometric ratio
//#define ONEPSFPERIMAGE

#ifdef STORAGE

class TabulatedDaoPsf : public Kernel {

private:
  double mx; // moments to compute second derivative
  double my; 
  double mxx; 
  double myy; 
  double mxy; 
  double det;
  double integral;
public:
  
  //! empty constructor allocate nothing
  TabulatedDaoPsf() {}
  
  //! allocate psf and derivatives of half-sizes Hx and Hy
  TabulatedDaoPsf(const int Hx, const int Hy) : Kernel(Hx,Hy), Dx(Hx,Hy), Dy(Hx,Hy)  {}
  
  //! allocate psf and derivatives of a Radius
  TabulatedDaoPsf(const int Radius) : Kernel(Radius), Dx(Radius), Dy(Radius)  {}
  
  //! allocate psf and derivatives, and fill them with DaoPsf value on that Pt.
  TabulatedDaoPsf(const Point& Pt, const DaoPsf& Psf, const int Radius);

  //! allocate psf and derivatives, and fill them with DaoPsf value on that Pt
  TabulatedDaoPsf(const Point& Pt, const DaoPsf& Dao, const Window& Rect);

  // default destructor, copy constructor and assigning operator are OK  

  //! tabulated derivative of the psf with x
  Kernel Dx;

  //! tabulated derivative of the psf with y
  Kernel Dy;
  
  void Resize(const int Hx, const int Hy);

  void Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius);

  void Tabulate(const Point& Pt, const DaoPsf& Dao, const Window& Rect);

  //! return the current norm of the psf
  double Norm() const { return integral; }

  double Mxx() const {return mxx;};
  double Mxy() const {return mxy;};
  double Myy() const {return myy;};
  

  //! enable "cout << TabulatedDaoPsf << endl"
  //friend ostream& operator << (ostream & stream, const TabulatedDaoPsf& Psf);

  void ComputeMoments();
  
  double Gaus(int i,int j) const;
  double dGausdx2(int i,int j) const;
  double dGausdy2(int i,int j) const;
  double dGausdxdy(int i,int j) const;
  
  void writeGaussian(const string& filename) ;

};
#endif

class TabulatedPsf : public Kernel {

private:
  double mx; // moments to compute second derivative
  double my; 
  double mxx; 
  double myy; 
  double mxy; 
  double det;
  double integral;
public:
  
  //! empty constructor allocate nothing
  TabulatedPsf() {}
  
  //! allocate psf and derivatives of half-sizes Hx and Hy
  TabulatedPsf(const int Hx, const int Hy) : Kernel(Hx,Hy), Dx(Hx,Hy), Dy(Hx,Hy)  {}
  
  //! allocate psf and derivatives of a Radius
  TabulatedPsf(const int Radius) : Kernel(Radius), Dx(Radius), Dy(Radius)  {}
  
  //! allocate psf and derivatives, and fill them with DaoPsf value on that Pt
  TabulatedPsf(const Point& Pt, const ImagePSF& imagepsf, const Window& Rect);

  // default destructor, copy constructor and assigning operator are OK  

  //! tabulated derivative of the psf with x
  Kernel Dx;

  //! tabulated derivative of the psf with y
  Kernel Dy;
  
  void Resize(const int Hx, const int Hy);

  
  void Tabulate(const Point& Pt, const ImagePSF& imagepsf, const Window& Rect);

  //! return the current norm of the psf
  double Norm() const { return integral; }

  double Mxx() const {return mxx;};
  double Mxy() const {return mxy;};
  double Myy() const {return myy;};
  void ComputeMoments();
};


class SimFitRefVignet : public Vignet {

public:
  
  //!
  //CountedRef<DaoPsf> daopsf;
  CountedRef<ImagePSF> imagepsf;
  

  //! build a rough initial galaxy to start the iterative fit
  void makeInitialGalaxy();

  //! empty constructor allocate nothing
  SimFitRefVignet(bool usegal=true) { ImagePSF *p = 0; imagepsf = p; UseGal = usegal; }
  
  //! does not do much
  SimFitRefVignet(const ReducedImage *Rim, bool usegal=true);
  
  // default destructor, copy constructor and assigning operator are OK

  //! tabulated PSF and its derivatives to allow fast computation
  //TabulatedDaoPsf Psf;
  TabulatedPsf Psf;
  
  //! the galaxy underneath at its resolution
  Kernel Galaxy;

  //!
  void UpdatePsfResid();


  void SetStar(const PhotStar *RefStar);
  //!
  void Load(const PhotStar *RefStar);

  void Resize(const int Hx, const int Hy);

  //! enable "cout << SimFitRefVignet << endl"
  //friend ostream& operator << (ostream& stream, const SimFitRefVignet& myVignet);

   bool UseGal;

};

//! A Vignet with tabulated PSF, and Kernel
class SimFitVignet : public Vignet {

private: 
  //! list of boolean to check whether components have been updated when the star is changed
  bool kernel_updated;
  bool psf_updated;
  bool resid_updated;
  bool gaussian_updated;
  
public:

#ifdef ONEPSFPERIMAGE
  CountedRef<DaoPsf> daopsf;
#endif

  KernelFit* kernelFit;

  //! empty constructor allocate nothing
  SimFitVignet();
  
  SimFitVignet(const ReducedImage *Rim,  SimFitRefVignet* Ref);

  //! extract vignet image and weight from a ReducedImage of a max radius of Nfwhm*Fwhm pixels
  //! loads the Kern with the Reference, and also builds the proper Psf.
  SimFitVignet( const PhotStar *Star, const ReducedImage *Rim,  SimFitRefVignet* Ref);
  
  virtual ~SimFitVignet() {delete kernelFit;};

  void ResetFlags();
  void ModifiedResid() {resid_updated = false;};
  
  // default destructor, copy constructor and assigning operator are OK

  bool FitFlux; // do we need to fit the flux
  bool FitPos; //  do we need to fit the position
  bool FitSky; // do we need to fit the sky bg
  bool UseGal; // do we need to use a model for the galaxy (does not mean we necessarly fit it)
  bool CanFitFlux; // do we need to fit the flux
  bool CanFitPos; //  do we need to fit the position
  bool CanFitSky; // do we need to fit the sky bg
  bool CanFitGal; // 
  bool DontConvolve;
  bool forceresize;
  double inverse_gain;
  double ronoise;
  double skysub;

  //! a kernel to convolve a reference to match the current vignet
  Kernel Kern;
  
  //! weight taking into account star flux
  Kernel OptWeight;

  //! tabulated PSF and its derivatives to allow fast computation
  TabulatedPsf Psf;

  //! reference to the simfitvignetref for updating size, star, and psf 
  CountedRef<SimFitRefVignet> VignetRef;
  
  
  void SetStar(const PhotStar *RefStar);
  
  //! resize vignet and call Update()
  void Resize( int Hx_new,  int Hy_new);
  
  // ! this makes a tabulated version of the reference image psf (that of VignetRef)
  void BuildPsf();

  //! this builds the kernel between this image and that of the reference image (that of VignetRef)
  void BuildKernel();

  //! update psf and residuals assuming the model = Star->flux[Kern*RefPsf] + Kern*RefGal + Star->sky
  void UpdateResid_psf_gal();

  //! update psf and residuals assuming the model = Star->flux[Kern*RefPsf] + Star->sky
  void UpdateResid_psf();

  //! update residuals assuming the model = Star->flux*Psf + Kern*RefGal + Star->sky
  void UpdateResid_gal();

  //! update residuals assuming the model = Star->flux*Psf + Star->sky
  void UpdateResid();
  
  //! update kernel, psf, residus if needed
  void Update();
  
  //! auto resize vignet according to the size of vignetref and call Update()
  void AutoResize();
  
  //! do not read images but fills todo list for dimage::readfitsimage()
  void PrepareAutoResize();
  
  //! use Vignet::Chi2 with additionnal debugging info
  double Chi2() const;
  
  void RedoWeight();

  //! 
  double CentralChi2(int &npix) const;
  
  
  //! write a lot of stuff, to use in case of problem 
  void DumpDebug() const;

  //! check whether there are weights>0 on the position of the star
  void CheckWeight();

  //! enable "cout << SimFitVignet << endl"
  //friend ostream& operator << (ostream & stream, const SimFitVignet& myVignet);
};


    
#endif // SIMFITVIGNET__H
