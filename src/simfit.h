// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMFIT__H
#define SIMFIT__H

#include "lightcurve.h"
#include "simfitvignet.h"
#include "matvect.h"

//!  \file simfit.h
//!  \brief A class to fit pixels simultaneously from a set of vignets.
//!
//!  The model on a given vignet to be fitted is 
//!  Model(i,j) = flux [Kernel X PSF(xc,yc)](i,j) + [Kernel X Gal](i,j) + sky
//!  where user chooses what to fit (flux, position, galaxy pixels and constant sky). 
//!  The flux and sky are free to change on each vignets, but the galaxy pixels 
//!  and the point source position are fixed at the best seeing vignet.
//!  The final fluxes are therefore expressed in the best seeing reference image.
//!  PSF and Kernel should be given and will not be re-fitted.

// what to fit
const unsigned int FitFlux = 1;
const unsigned int FitPos  = 2;
const unsigned int FitGal  = 4;
const unsigned int FitSky  = 8;

// what result to write
const unsigned int WriteLightCurve = 1;
const unsigned int WriteGalaxy  = 2;
const unsigned int WriteResid  = 4;
const unsigned int WriteData  = 8;
const unsigned int WritePsf  = 16;
const unsigned int WriteWeight  = 32;
const unsigned int WriteKernel  = 64;
const unsigned int WriteVignetsInfo  = 128;
const unsigned int WriteMatrices  = 256;


typedef ImageList<SimFitVignet>::iterator SimFitVignetIterator;
typedef ImageList<SimFitVignet>::const_iterator SimFitVignetCIterator;

//! Simultaneous fitting of Vignets
class SimFit : public vector< CountedRef<SimFitVignet> > {

private:
  // flags
  bool refill;            // whether we need to refill the gal-gal matrix part
  bool fit_flux;          // whether we fit the point source fluxes
  bool fit_gal;           // whether we fit the galaxy
  bool fit_sky;           // whether we fit the sky
  bool fit_pos;           // whether we fit the point source position
  bool use_gal;           // one can use the galaxy model but not fit it
  bool dont_use_vignets_with_star; // this when you want to fit the galaxy only, see 
  bool fatalerror; // internal bool to quit without core dump
  bool inverted; //check if matrix was inverted

  // vector and matrices for the system Mat*Params=Vec
  Vect Vec;            // vector r.h.s and Params when solved
  Mat PMat;            // matrix l.h.s then covariance matrix when inverted
  Mat MatGal;         // gal-gal matrix part to avoid refilling
  Mat NightMat;      // see fillNightMat

  // indices
  int fluxstart, fluxend; // start and end indices for flux parameters in Mat and Vec
  int xind,yind;          // indices for positional parameters in Mat and Vec
  int galstart, galend;   // start and end indices for galaxy parameters in Mat and Vec
  int skystart, skyend;   // start and end indices for sky parameters in Mat and Vec
  
  // various numbers
  int nfx,nfy,hfx,hfy;    // galaxy sizes and half-sizes
  int nparams, ndata;     // current number of params, data
  double scale,minscale;  // current and minimum scaling factor 
  double chi2;            // current chi2

  // returns the galaxy matrix index given pixel (i,j)
  inline int galind(const int i, const int j) const;

  // Mat and Vec filling routines
  void fillFluxFlux();
  void fillFluxPos();
  void fillFluxGal();
  void fillFluxSky();
  void fillPosPos();
  void fillPosGal();
  void fillPosSky();
  void fillGalGal();
  void fillGalSky();
  void fillSkySky();
  
  // compute the chi2 of the current fit
  double computeChi2() const;

  // perform one Newton-Raphson iteration: fill system and solve, check decreasing of chi2
  double oneNRIteration(double oldchi2);

  // handling of errors
  void FatalError(const char* comment);

public:

  //! simply initialize properly the many private members
  SimFit();

  // default destructor, copy constructor and assigning operator are OK  

  //! reference to the best seeing vignet
  CountedRef<SimFitRefVignet> VignetRef;

  //! set what you want to fit
  void SetWhatToFit(unsigned int ToFit = FitFlux);

  //! get worst seeing value in the list of images
  double GetWorstSeeing();

  //! get the minimum scaling factor to resize the vignets. WorstSeeing is in ReducedImage::Seeing() unit
  void FindMinimumScale(double WorstSeeing);

  //! resize all the vignets of a scale factor, resize matrixes, and compute indices
  void Resize(const double& ScaleFactor);

  //! allow to change full data set to another star
  void Load(LightCurve& Lc, bool keepstar=false, bool only_reserve_images=false);

  //! fill the entire matrix and vectors
  void FillMatAndVec();

  //! iterate on solution and solve the system
  bool IterateAndSolve(int MaxIter=10, double Eps=0.01);

  //! a procedure to fill up the covariance into the Mat and the proper SimFitVignets...
  bool GetCovariance();

  //! update the vignets with the current solution, possibily apply a scale factor to the solution
  bool Update(double Factor=1., bool print=true);

  //! fill, fit and shit: do everything
  bool DoTheFit(int MaxIter=10, double epsilon=0.01);

  //! returns the chi2 of the current fit
  double Chi2() const { return chi2; }

  //! returns the number of degrees of freedom of the current fit
  int Dof() const { return ndata-nparams; }

  //! sigma scale for the variance parameters
  double VarScale() const;

  //! returns the total flux for the galaxy
  double GalFlux() const;

  //! returns the variance for the total flux of the galaxy
  double VarGalFlux() const;

  //! returns the total flux for all the nights
  double TotFlux() const;

  //! returns the variance for TotFlux()
  double VarTotFlux() const;
  
  //! returns the total sky for all the nights
  double TotSky() const;

  //! returns the variance for TotSky()
  double VarTotSky() const;

  //! gets the mean, median, rms and absolute deviation for all the residual pixels
  double MeanMedianRms(double& med, double& rms, double& adev) const;

  //! returns the current used scale factor for vignets
  double Scale() const { return scale; }

  //! write galaxy, covariance matrix and lightcurve on disk
  void write(const string &StarName,const string &DirName=".", unsigned int whattofit = 0);
  
  ostream& DumpMatrices(ostream& Stream=cout) const; 
  
  void UseGalaxyModel(bool useit = true);
  
  //! fit initial galaxy using only vignets without burning star
  void FitInitialGalaxy();
  
  //! fill a boolean matrix
  //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  void fillNightMat(LightCurve& Lc);
  
  //! fill a boolean matrix
  //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  void fillNightMat();

  //! DumpAndAbort
  void DumpAndAbort(const string& message) const ;
  
  //! Check Consistency of vignets
  bool CheckConsistency() const;

  void Chi2PosMap();
  

  //! enable "cout << SimFit << endl;"
  friend ostream& operator << (ostream& Stream, const SimFit& SimFit);
  
};

#endif // SIMFIT__H
