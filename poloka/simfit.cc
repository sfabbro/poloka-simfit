#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <ctime>
#include <algorithm>
#include <poloka/simfitvignet.h>
#include <poloka/simfit.h>
 
#define DEBUG

// #define ONLYPOSITIVEFLUXFORPOSITION
//#define USE_SECOND_DERIVATIVE_OF_POSITION
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  :::::::::::::::::: SimFit stuff   ::::::::::::::::::::::::::
  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

SimFit::SimFit()
{
  fatalerror = false;
  refill = true;
  fit_flux = fit_gal = true; fit_sky = fit_pos = false ;
  use_gal = true;
  fluxstart = galstart = skystart = xind = yind = 0;
  fluxend = galend = skyend = 0;
  scale = 1., minscale = 0.; chi2 = 1e29;
  nfx = nfy = hfx = hfy = nparams = ndata = 0;
  dont_use_vignets_with_star = false;
  inverted = false;
}

void SimFit::UseGalaxyModel(bool useit) {
  use_gal = useit;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    (*itVig)->UseGal = use_gal;
  }
}

void SimFit::SetWhatToFit(unsigned int ToFit)
{
  bool resize = false;
  if (fit_flux != (ToFit & FitFlux)) { fit_flux = (ToFit & FitFlux); resize=true; }
  if (fit_pos  != (ToFit & FitPos))  { fit_pos  = (ToFit & FitPos ); resize=true; }
  if (fit_gal  != (ToFit & FitGal))  { fit_gal  = (ToFit & FitGal ); resize=true; }
  if (fit_sky  != (ToFit & FitSky))  { fit_sky  = (ToFit & FitSky ); resize=true; }
  use_gal |= fit_gal; // use the galaxy if we fit it
  dont_use_vignets_with_star = false;

  int count=0;
  cout << " > SimFit::SetWhatToFit(): fit_flux: " << fit_flux 
       << " fit_pos: " << fit_pos 
       << " fit_gal: " << fit_gal 
       << " fit_sky: " << fit_sky << endl;
  //cout << "What to fit for each vignet (id FitPos UseGal FitSky FitFlux): " << endl;
  bool firstwithsky = true;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    (*itVig)->FitPos = fit_pos && (*itVig)->CanFitPos;
    (*itVig)->UseGal = use_gal && (*itVig)->CanFitGal;
    // micmac pour fixer un niveau de fond de ciel afin de lever la degenerescence avec la galaxie
    
    if( fit_sky && (*itVig)->CanFitSky && firstwithsky && fit_gal) { 
      (*itVig)->FitSky = false;
      firstwithsky = false;
    }else{
      (*itVig)->FitSky =  fit_sky && (*itVig)->CanFitSky;
    }
    (*itVig)->FitFlux = fit_flux && (*itVig)->CanFitFlux;
    //printf("%03d %d %d %d %d\n",count,(*itVig)->FitPos,(*itVig)->UseGal,(*itVig)->FitSky,(*itVig)->FitFlux);
    count++;
  }
}

void SimFit::FindMinimumScale(double WorstSeeing)
{
#ifdef FNAME
  cout << " > SimFit::FindMinimumScale() : WorstSeeing =" << WorstSeeing  << endl;
#endif

  int hmin = max(int(ceil(WorstSeeing*2.3548)), 5);
  int hrefx = VignetRef->Data.HSizeX();
  int hrefy = VignetRef->Data.HSizeY();
  int hkx = 0;
  int hky = 0;

#ifdef DEBUG
  cout << "VignetRef->Data.HSizeX() = " << hrefx << endl;
  cout << "VignetRef->Data.HSizeY() = " << hrefy << endl;
#endif

  for (SimFitVignetIterator it = begin(); it != end(); ++it)
    {
      SimFitVignet *vi = *it;
      if ((vi->Kern.HSizeX() > hkx) || (vi->Kern.HSizeY()> hky))
	{
	  hkx = vi->Kern.HSizeX();
	  hky = vi->Kern.HSizeY();
	  hrefx = max(hmin, hkx) + hkx;
	  hrefy = max(hmin, hky) + hky;
	}
    }

  minscale = double(min(hrefx,hrefy)) / double(min(VignetRef->Data.HSizeX(),VignetRef->Data.HSizeY()));

#ifdef DEBUG
  cout << " SimFit::FindMinimumScale(" << WorstSeeing 
       << ") : Minimum scaling factor = " << minscale << endl;
#endif

}

double SimFit::GetWorstSeeing() {
#ifdef FNAME
  cout << " > SimFit::GetWorstSeeing()" << endl;
#endif 
  double worstSeeing = 0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    double currentseeing = (*itVig)->Seeing();
    
    if(currentseeing>worstSeeing)
      worstSeeing=currentseeing;
  }
  return worstSeeing;
}

void SimFit::Load(LightCurve& Lc, bool keepstar, bool only_reserve_images)
{
#ifdef FNAME
  cout << " > SimFit::Load()" << endl;
#endif
  
  fatalerror = false; // reset

  if (Lc.size() != size()) 
    {
      cerr << " > SimFit::Load() Error : trying to load a LightCurve of different size \n";
      return;
    }

  cout << " > SimFit::Load() : " << Lc.Ref->name << " at x,y = "
       << Lc.Ref->x << ", " << Lc.Ref->y << endl;
  
  // define the size of the reference vignet
   
  // get worst seeing to set the size of reference vignette
  double worst_seeing = GetWorstSeeing();
  
   
  // now build kernels of vignets to get the worst kernel
  int worst_kernel = 0;
  LightCurve::const_iterator itLc = Lc.begin();
  int nim=0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig, ++itLc)
    {      
      SimFitVignet *vi = *itVig;
      cout << flush << " > SimFit::Load() : loading vignets from image #" << nim++ << "\r";
      //vi->Load(*itLc); // this does not modify the kernel so I would better use SetStar which does nothing but set the star
      if(!keepstar) {
	vi->SetStar(*itLc); // just set the star
	vi->BuildKernel(); // rebuild kernel
      } 
      if(vi->Kern.HSizeX()>worst_kernel)
	worst_kernel = vi->Kern.HSizeX();
      if(vi->Kern.HSizeY()>worst_kernel)
	worst_kernel = vi->Kern.HSizeY();
    }
  cout << endl;

  // 2.3548*sigma = full-width at half-maximum [2.3548 = 2.*sqrt(2*log(2.))]
  // seeing (from sextractor SESEEING) is sigma in pixel units
  
  // radius is the size of the reference vignet
  float nseeing = 2.3548; // number of seeings for the vignet radius size
  //float nseeing = 1; // number of seeings for the vignet radius size
  int radius = int(ceil(nseeing*worst_seeing+worst_kernel)); 
  // minscale  = min_radius/radius (min_radius is used for fitting the position)
  minscale = (worst_seeing+worst_kernel)/radius;
  
#ifdef DEBUG
  cout << " > SimFit::Load() : worst_seeing = " << worst_seeing 
       << " worst_kernel = " << worst_kernel
       << " radius = " << radius
       << " minscale = " << minscale << endl;
#endif
  
  // the VignetRef has already been build
  if(!keepstar)
    VignetRef->SetStar(Lc.Ref); // just set the star
  VignetRef->Resize(radius,radius); // now resize, this reloads data, update psf, and makeInitialGalaxy if usegal
#ifdef DEBUG
  cout << " > SimFit::Load() : VignetRef() half size (" 
       << VignetRef->Hx() << ", " << VignetRef->Hy() << ")\n";
#endif 
  
  if(only_reserve_images) {
    nim = 0;
    for (SimFitVignetIterator it = begin(); it != end(); ++it) {
      SimFitVignet *vi = *it;
      cout << flush << " > SimFit::Load() : reserving vignets #" << nim++ << "\r";
      vi->PrepareAutoResize();
    }
    cout << endl;
    return;
  }
  
  // now resize all vignets and initialize residuals
  nim = 0;
  for (SimFitVignetIterator it = begin(); it != end(); ++it)
    {
      SimFitVignet *vi = *it;
      cout << flush << " > SimFit::Load() : checking vignets #" << nim++ << "\r";
      // we will fit the flux according to Lc.Ref and date
      vi->CanFitFlux = Lc.Ref->IsVariable(vi->ModifiedJulianDate());

      vi->CanFitPos  = vi->CanFitFlux;
      //vi->DontConvolve = Lc.Ref->Image()->Name() == vi->Image()->Name(); 
      vi->DontConvolve = false;
      // this resizes the vignet and update everything (kernel, psf, residus) 
      vi->ModifiedResid();
      vi->AutoResize();
      // check whether there are weights>0 on the position of the star,
      // if not, FitFlux, FitSky =false 
      // shunt this (only for simu!!)
      vi->CheckWeight();
    }
  cout << endl;
  
}

void SimFit::Resize(const double& ScaleFactor)
{

#ifdef FNAME
  cout << " > SimFit::Resize() : " << ScaleFactor << endl;
#endif
  
  int hrefx = VignetRef->Hx();
  int hrefy = VignetRef->Hy();
  
  //if(true) {
  if(fabs(ScaleFactor-1)>0.01) {
    if ((!fit_flux) && (!fit_pos) && (!fit_gal) && (!fit_sky) || (size()==0)) 
      {
#ifdef DEBUG
	cerr << " > SimFit::Resize() Warning : nothing to fit, not resizing. " << endl;
#endif
	return;
      }
    
    scale = max(ScaleFactor, minscale);

#ifdef DEBUG    
    cout << " > SimFit::Resize() : scalefactor = " << ScaleFactor << " old scale = " << scale
	 << " new scale = " << scale << endl;
#endif
  
    // resize vignets
    
    VignetRef->Resize(int(ceil(scale*double(VignetRef->Hx()))),
		      int(ceil(scale*double(VignetRef->Hy()))));
  }
  hrefx = VignetRef->Hx();
  hrefy = VignetRef->Hy();
  
  // anyway, resize and update everything for all vignets
  for (SimFitVignetIterator it = begin(); it != end(); ++it) {
    (*it)->ModifiedResid();
    (*it)->AutoResize();
  }
  
  // recompute matrix indices
  hfx = hfy = nfx = nfy = 0;
  fluxstart = fluxend = 0;
  
  // fluxes
  if (fit_flux)
    {
      for (SimFitVignetCIterator it=begin(); it != end(); ++it)
	if ((*it)->FitFlux) fluxend++;
    }  
  fluxend--; // we want fluxend = nfluxes-1
    
  
  // position
  xind = yind = fluxend;
  if (fit_pos) 
    {
      xind += 1;
      yind += 2;
    }
  
  // galaxy pixels
  galstart = galend = yind;
  if (fit_gal)
    {
      if(fluxend>=0)
	galstart += 1;
      hfx = hrefx;
      hfy = hrefy;
      nfx = 2*hfx+1;
      nfy = 2*hfy+1;
      galend = galstart + nfx*nfy-1;
      refill = true;
    }
  
  // skies
  skystart = skyend = galend;
  if (fit_sky)
    {
      skystart += 1;
      for (SimFitVignetCIterator it=begin(); it != end(); ++it)
	if ((*it)->FitSky) skyend++;
    }
  
  
  nparams = (fluxend-fluxstart+1) + (yind-fluxend) + (galend-yind) + (skyend-galend);

  ndata = 0;
  
  for (SimFitVignetCIterator it = begin(); it != end(); ++it) {
    if((*it)->CanFitFlux && dont_use_vignets_with_star)
      continue;
    if(! (*it)->CanFitFlux && ! (*it)->UseGal && ! (*it)->CanFitSky && ! (*it)->CanFitPos  )
      continue;
    ndata += (*it)->NValidPixels();
  }
  //ndata += (2*(*it)->Hx()+1) * (2*(*it)->Hy()+1);
  
#ifdef DEBUG
  cout << " > SimFit::Resize() : "
       << " #images = " <<  size() 
       << " #data = " <<  ndata
       << " #params = " << nparams
       << " #flux = " << (fluxend-fluxstart+1)
       << " #pos = " << (yind-fluxend)
       << " #gal = " << (galend-yind)
       << " #sky = " << (skyend-galend)
       << endl;
  cout << " > SimFit::Resize() : "
       << " flux[" << fluxstart  << ":" << fluxend
       << "] pos[" << xind << ":" << yind
       << "] gal[" << galstart  << ":" << galend
       << "] sky[" << skystart << ":" << skyend 
       << "]\n";
#endif

  if ((fluxend-fluxstart+1 <= 0) && (fit_flux)) {
    FatalError("Want to fit flux but nothing to fit");
    return;
  }

  if ((yind-fluxend <= 0) && (fit_pos)) {
    FatalError("Want to fit pos but nothing to fit");
    return;
  }

  if ((galstart-yind <= 0) && (fit_gal)) {
    FatalError("Want to fit gal but nothing to fit");
    return;
  }

  if ((skyend-galend <= 0) && (fit_sky)) {
    FatalError("Want to fit sky but nothing to fit");
    return;
  }

  Vec.allocate(nparams);
  PMat.allocate(nparams, nparams);
  MatGal.allocate(nfx*nfy,nfx*nfy);
  
#ifdef DEBUG
  cout << "   nparams = " << nparams << endl;
  cout << "   nfx     = " << nfx << endl;
  cout << "   nfy     = " << nfy << endl;
  cout << "   fluxstart,fluxend  = " << fluxstart << "," << fluxend << endl;
  cout << "   xind,yind          = " << xind << "," << yind << endl;
  cout << "   galstart,galend    = " << galstart << "," << galend << endl;
  cout << "   skystart,skyend    = " << skystart << "," << skyend << endl;  
#endif

}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::Matrix filling routines::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IMPORTANT: right now, coordinates are all centered on (0,0), center of 
  the vignet. That makes it easier to read. If you change the origin for 
  the data, weights or psf kernels, you should recode everything. Do the 
  merguez barbecue instead. 
*/


void SimFit::FillMatAndVec()
{

#ifdef FNAME
  cout << " > SimFit::FillMatAndVec() : Initialize " << endl;  
#endif
  inverted = false;
  Vec.Zero();
  PMat.Zero();
  
#ifdef DEBUG
  cout << " > SimFit::FillMatAndVec() : Compute matrix and vectors " << endl;  
#endif

  if (fit_flux)            fillFluxFlux();
  if (fit_flux && fit_pos) fillFluxPos();
  if (fit_flux && fit_gal) fillFluxGal();
  if (fit_flux && fit_sky) fillFluxSky();

  if (fit_pos)             fillPosPos();
  if (fit_pos && fit_gal)  fillPosGal();
  if (fit_pos && fit_sky)  fillPosSky();

  if (fit_gal)             fillGalGal();
  if (fit_gal && fit_sky)  fillGalSky();

  if (fit_sky)             fillSkySky();

  // no need to symmetrize the matrix

#ifdef DEBUG  
  if(PMat.SizeX()<20 && PMat.SizeY()<20) {
    cout << "===== PMat =====" << endl;
    cout << PMat << endl;
    cout << "===== PMat =====" << endl;
  }
#endif
}

int SimFit::galind(const int i, const int j) const
{
  return galstart + (i+hfx)*nfy + (j+hfy);
}

void SimFit::fillFluxFlux()
{
  //*********************************************
  // flux-flux matrix terms and flux vector terms
  //*********************************************

#ifdef FNAME
  cout << " > SimFit::fillFluxFlux()" << endl;
#endif

  // loop over vignets
  int fluxind = 0;
  DPixel *pw, *pres, *ppsf;
  double sumvec, summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;
      
      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      sumvec = 0.;
      summat = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
#ifdef DEBUG_TOTO
	  cout << "pw,pres,ppsf " << *pw << "," << *pres << "," << *ppsf << endl;
#endif

	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvec += (*pres) * (*ppsf) * (*pw);
	      summat += (*ppsf) * (*ppsf) * (*pw);
	      ++ppsf; ++pres; ++pw;
	    }
	}

      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;

      Vec(ind) = sumvec;
      PMat(ind,ind) = summat;
      
      ++fluxind;
    }
}

void SimFit::fillFluxPos()
{
  //**********************
  // flux-pos matrix terms
  //**********************

#ifdef FNAME
  cout << " > SimFit::fillFluxPos()" << endl;
#endif

  // loop over vignets
  int fluxind = 0;
  DPixel *pw, *ppsf, *ppdx, *ppdy;
  double summatx, summaty;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;

      if (!vi->FitFlux || !vi->FitPos) continue;

#ifdef ONLYPOSITIVEFLUXFORPOSITION
      if(vi->Star->flux<=0)  {
	fluxind++;
	continue;
      }
#endif


      // sum over the smallest area between psf and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      summatx = 0.;
      summaty = 0.;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
	  ppdx = &(vi->Psf.Dx)(-hx,j);
	  ppdy = &(vi->Psf.Dy)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillFluxPos pw,ppsf,ppdx,ppdy " 
	       << *pw << "," 
	       << *ppsf << ","
	       << *ppdx << ","
	       << *ppdy << endl;
	    
#endif
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summatx += (*ppsf) * (*ppdx) * (*pw);
	      summaty += (*ppsf) * (*ppdy) * (*pw);
	      ++pw; ++ppsf; ++ppdx; ++ppdy;
	    }
	}
      
      summatx *= vi->Star->flux; //JG
      summaty *= vi->Star->flux; //JG
      
      // now fill in the matrix and vector
      int ind = fluxstart+fluxind;


      PMat(xind,ind) = summatx; // ok cause xind>ind
      PMat(yind,ind) = summaty; // ok cause ymin>ind

      fluxind++;
    }  
}


#define KERNIND(HK,H,IND,A,B)\
int At = (HK-IND < H) ? IND-HK : -H;\
int B = (HK+IND < H) ? IND+HK : H;\
int A = At<B ? At : B;\
B = At>B ? At : B;\

/*
#define KERNIND(HK,H,IND,A,B)\
int A = (HK-IND < H) ? IND-HK : -H;\
int B = (HK+IND < H) ? IND+HK : H;\
A = A<B ? A : B;\
B = A>B ? A : B;\
*/

void SimFit::fillFluxGal()
{ 
  //*********************
  // flux-gal matrix part
  //*********************

#ifdef FNAME
  cout << " > SimFit::fillFluxGal()" << endl;
#endif


  int fluxind = 0;
  DPixel *pw, *ppsf, *pkern;
  double summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitFlux) continue;
      if (!vi->UseGal) {
	++fluxind;
	continue;
      }
      
      
      int hx = vi->Hx();
      int hy = vi->Hy();
  
      int ind = fluxstart+fluxind;

      // the dirac case : sum[ psf(x) * dirac(y) ] = psf(y)	  
      if (vi->DontConvolve)
	{
	  // loop over fitting coordinates 
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      ppsf = &(vi->Psf)   (-hx,j);
	      pw   = &(vi->OptWeight)(-hx,j);
	      for (int i=-hx; i<=hx; ++i)
		{
		  PMat(galind(i,j),ind) = (*ppsf) * (*pw); // ok cause galind > ind
		  ++ppsf; ++pw;
		}
	    }
	  ++fluxind;
	  continue;
	}
      
      // the normal cases
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      // loop over fitting coordinates
      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      summat = 0.;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  ppsf  = &(vi->Psf)   (ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summat += (*ppsf) * (*pkern) * (*pw);
		      ++pkern; ++ppsf; ++pw;
		    }
		}
	      PMat(galind(is,js),ind) = summat;  // ok cause galind > ind
	    }
	}
      ++fluxind;

    }
}

void SimFit::fillFluxSky()
{
  //***********************
  // flux-sky matrix terms
  //***********************
  
#ifdef FNAME
  cout << " > SimFit::fillFluxSky()" << endl;
#endif

  int skyind  = 0;
  int fluxind = 0;
  double summat;
  DPixel *pw, *ppsf;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (vi->FitFlux && vi->FitSky) {

      // sum over smallest between psf and data coordinates
      int hx = vi->Hx();
      int hy = vi->Hy();
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  ppsf = &(vi->Psf)   (-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      summat += (*ppsf) * (*pw);
	      ++ppsf; ++pw;
	    }  
	}
      // now fill in matrix part

      //PMat(fluxstart+fluxind,skystart+skyind) = summat; // that's wrong 
      PMat(skystart+skyind,fluxstart+fluxind) = summat; // correct cause skystart+skyind > fluxstart+fluxind
      
      }
      if(vi->FitFlux)
	++fluxind; 
      if(vi->FitSky)
	++skyind;
      
    }
  //cout << "in fillFluxSky() skyind = " << skyind << endl;
}

void SimFit::fillPosPos()
{
  //*********************************************
  // pos-pos matrix terms and pos vector terms
  //*********************************************

#ifdef FNAME
  cout << " > SimFit::fillPosPos()" << endl;
#endif

  double sumvecx = 0.;
  double sumvecy = 0.;
  double summatx = 0.;
  double summaty = 0.;
  double summatxy = 0.;
  DPixel *pw, *pres, *ppdx, *ppdy;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitPos) continue;

      // sum over the smallest area between psf derivative and vignet
      int hx = vi->Hx();
      int hy = vi->Hy();

      double flux = vi->Star->flux;
#ifdef ONLYPOSITIVEFLUXFORPOSITION
      if(flux<=0) continue;
#endif
      double flux2 = flux*flux;

      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  ppdx = &(vi->Psf.Dx)(-hx,j);
	  ppdy = &(vi->Psf.Dy)(-hx,j);
	  
	  
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillPosPos pw,pres,ppdx,ppdy " 
	       << *pw << "," 
	       << *pres << ","
	       << *ppdx << ","
	       << *ppdy << endl;
#endif
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvecx  += (*pres) * (*ppdx) * (*pw) * flux;
	      sumvecy  += (*pres) * (*ppdy) * (*pw) * flux;
	      summatx  += (*ppdx) * (*ppdx) * (*pw) * flux2;
	      summaty  += (*ppdy) * (*ppdy) * (*pw) * flux2;
	      summatxy += (*ppdx) * (*ppdy) * (*pw) * flux2;
#ifdef USE_SECOND_DERIVATIVE_OF_POSITION
	      // new : use second derivative of pos 
	      summatx  += vi->Psf.dGausdx2(i,j) * (*pres) * (*pw) * flux;
	      summaty  += vi->Psf.dGausdy2(i,j) * (*pres) * (*pw) * flux;
	      summatxy += vi->Psf.dGausdxdy(i,j) * (*pres) * (*pw) * flux;
#endif	      
	      ++pw; ++pres; ++ppdx; ++ppdy;
	    }
	}    
    }
  
  Vec(xind) = sumvecx;
  Vec(yind) = sumvecy;
  PMat(xind,xind) = summatx;
  PMat(yind,yind) = summaty;
  PMat(yind,xind) = summatxy;
}

void SimFit::fillPosGal()
{
  //*********************
  // pos-gal matrix part
  //*********************


#ifdef FNAME
  cout << " > SimFit::fillPosGal()" << endl;
#endif

  DPixel *ppdx, *ppdy, *pw, *pkern;
  double summatx, summaty;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (!vi->FitPos) continue;

#ifdef ONLYPOSITIVEFLUXFORPOSITION
      if(vi->Star->flux<=0) continue;
#endif

      if (!vi->UseGal) continue;
      
      int hx = vi->Hx();
      int hy = vi->Hy();
  
      // the dirac case
      if (vi->DontConvolve)
	{
	  // loop over psf deriv:sum[ dpdx(x) * dirac(y) ] = dpdx(y)	  
	  for (int j=-hy;  j<=hy; ++j)
	    {
	      ppdx = &(vi->Psf.Dx)(-hx,j);
	      ppdy = &(vi->Psf.Dy)(-hx,j);
	      pw   = &(vi->OptWeight)(-hx,j);
#ifdef DEBUG_FILLMAT
	  cout << "   SimFit::fillPosGal pw,ppdx,ppdy " 
	       << *pw << "," 
	       << *ppdx << ","
	       << *ppdy << endl;
#endif
	      for (int i=-hx; i<=hx; ++i)
		{
		  PMat(galind(i,j),xind) += (*ppdx) * (*pw); // ok cause galind > xind
		  PMat(galind(i,j),yind) += (*ppdy) * (*pw); // ok cause galind > yind
		  ++ppdx; ++ppdy; ++pw;
		}
	    }
	  continue;
	}

      // loop over psf deriv in fitting coordinates
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      KERNIND(hky,hy,js,jkstart,jkend);
	      summatx = 0.;
	      summaty = 0.;
	      // sum over kernel centered on (i,j)
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  ppdx = &(vi->Psf.Dx)(ikstart,jk);
		  ppdy = &(vi->Psf.Dy)(ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      summatx += (*ppdx) * (*pkern) * (*pw);
		      summaty += (*ppdy) * (*pkern) * (*pw);
		      ++pkern; ++ppdx; ++ppdy; ++pw;
		    }
		}
	      
	      summatx *= vi->Star->flux;
	      summaty *= vi->Star->flux;
	      
	      PMat(galind(is,js),xind) += summatx; // ok cause galind > xind
	      PMat(galind(is,js),yind) += summaty; // ok cause galind > yind    

	    }
	} // end of loop on pixels

    } //end of loop on vignets
}

void SimFit::fillPosSky()
{
  //***********************
  // pos-sky matrix terms
  //***********************

#ifdef FNAME
  cout << " > SimFit::fillPosSky()" << endl;
#endif

  int skyind = 0;
  DPixel *pw, *ppdx, *ppdy;
  double summatx, summaty;

  // loop over vignets
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if (vi->FitSky && vi->FitPos) {

#ifdef ONLYPOSITIVEFLUXFORPOSITION
	if(vi->Star->flux<=0) {
	  ++skyind;
	  continue;
	}
#endif

	// sum over the smallest area between psf and vignet
	int hx = vi->Hx();
	int hy = vi->Hy();
	summatx = 0.;
	summaty = 0.;
	
	for (int j=-hy; j<=hy; ++j)
	  {
	    pw   = &(vi->OptWeight)(-hx,j);
	    ppdx = &(vi->Psf.Dx)(-hx,j);
	    ppdy = &(vi->Psf.Dy)(-hx,j);
	    for (int i=-hx; i<=hx; ++i) 
	      {
		summatx += (*ppdx) * (*pw);
		summaty += (*ppdy) * (*pw);
		++ppdx; ++ppdy; ++pw;
	      }
	  }
	
	summatx *= vi->Star->flux;
	summaty *= vi->Star->flux;
	
	// now fill matrix part
	PMat(skystart+skyind,xind) = summatx; // ok cause skystart+skyind > xind
	PMat(skystart+skyind,yind) = summaty; // ok cause skystart+skyind > yind
	
	++skyind;
      }
    }
}

void SimFit::fillGalGal()
{

  //******************************************
  // gal-gal matrix terms and gal vector terms
  //******************************************

#ifdef FNAME
  cout << " > SimFit::fillGalGal()" << endl;
#endif

  double sumvec;
  DPixel *pkern, *pres, *pw;

  // loop over vignets
#ifdef DEBUG_FILLMAT
  cout << "  Loop over vignets ..." << endl;
#endif
  int count = 0;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      if((*it)->CanFitFlux && dont_use_vignets_with_star)
	continue;

      if( ! (*it)->UseGal) {
	continue;
      }
      count ++;

      const SimFitVignet *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();

      // the dirac case
      if (vi->DontConvolve)
	{
	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		Vec(galind(i,j)) += (vi->Resid)(i,j) * (vi->OptWeight)(i,j);
	      }
	  continue;
	}


      /* galaxy contribution to vector terms
	 basically a convolution. Can't just simply use the Convolve routine: 
	 gotta choose some rules for the borders, depending on sizes user has chosen */

      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
      int hsy = (hy + hky) > hfy ? hfy : (hy + hky);

      for (int is=-hsx;  is<=hsx; ++is)
	{
	  KERNIND(hkx,hx,is,ikstart,ikend);
	  int ikstartis = ikstart-is;
	  for (int js=-hsy;  js<=hsy; ++js)
	    {
	      sumvec = 0.;
	      KERNIND(hky,hy,js,jkstart,jkend);
	      
	      // sum over kernel, stay in fitting coordinates
	      for (int jk=jkstart; jk<=jkend; ++jk)
		{
		  pkern = &(vi->Kern)  (ikstartis,jk-js);
		  pres  = &(vi->Resid) (ikstart,jk);
		  pw    = &(vi->OptWeight)(ikstart,jk);

#ifdef DEBUG_FILLMAT
		  //cout << "   fillGalGal pkern,pres,pw " << *pkern << "," << *pres << "," << *pw << endl;
#endif
		  for (int ik=ikstart; ik<=ikend; ++ik)
		    {
		      sumvec += (*pkern) * (*pres) * (*pw);
		      ++pkern; ++pres; ++pw;
		    }
		}
	      Vec(galind(is,js)) += sumvec;
	    }
	}
    }
#ifdef DEBUG
  cout << " > SimFit::fillGalGal() : nvignets in vector (galgal) = " << count << endl; 
#endif
  // now if weights do not change (that is not robustify), 
  // we do not need to refill this part at each iteration
  if (false && !refill) 
    {

      int ngal = galend-galstart+1; 
#ifdef DEBUG
      cout << "     case notrefile " << endl;
#endif
      for (int j=0; j<ngal; ++j) 
	for (int i=j+1; i<ngal; ++i) {
	  PMat(galstart+i,galstart+j) = MatGal(i,j);
	}
      return;
    }
  
  /* galaxy-galaxy matrix terms: longest loop.
     use the simplified relation, which in one dimension can be written as:
     matrix(m,n) = sum_ik ker[ik-(in-im)] * ker(ik) * weight(ik+im). 
     This is also a convolution, but again, borders are messy. */
  
#ifdef DEBUG
  cout << " > SimFit::fillGalGal() : Computing galaxy-galaxy matrix terms (longest loop) ... " << endl;
#endif
  count = 0;
  
  double summat;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {

   
       if((*it)->CanFitFlux && dont_use_vignets_with_star)
	continue;
       if( ! (*it)->UseGal) {
	 continue;
       }
       

       count++;

      const SimFitVignet *vi = *it;
      int hx = vi->Hx();
      int hy = vi->Hy();
      
      // the dirac case
      if (vi->DontConvolve)
	{

	  for (int j=-hy; j<=hy; ++j)
	    for (int i=-hx; i<=hx; ++i)
	      {
		int ind = galind(i,j);
		PMat(ind,ind) += (vi->OptWeight)(i,j);
	      }
	  continue;
	}

   
      int hkx = vi->Kern.HSizeX();
      int hky = vi->Kern.HSizeY();
      int min_m = -2*(hkx*nfy + hky);
      int npix = nfx*nfy;
      
      int m,n; // index of the galgal matrix
      int im,jm; // index of pixels in the galaxy image, corresponding to matrix index m, 
      int in,jn; // index of pixels in the galaxy image, corresponding to matrix index n
      // m = (im+hfx)*nfy + (jm+hfy)  same as galind
      int imn,jmn;
      int ik,jk;
      
      int ikmin,ikmax,jkmin,jkmax;
      DPixel *pkern1,*pkern2,*pw;
      
      // loop over fitting coordinates
      for ( m=0; m<npix; ++m)
	{
	  im = (m / nfy) - hfx;
	  jm = (m % nfy) - hfy;
	  for ( n=((min_m+m > 0) ? min_m+m : 0); n<=m; ++n) // ???
	    {
	      in = (n / nfy) - hfx;
	      jn = (n % nfy) - hfy;
	      imn = im-in;
	      jmn = jm-jn;
	      summat = 0.;
	      jkmin = max(max(-hky,-hy-jm),-hky-jmn);
	      jkmax = min(min(hky,hy-jm), hky-jmn);
	      
	      for (jk=jkmin; jk<=jkmax; ++jk) {
		ikmin = max(max(-hkx,-hx-im),-hkx-imn);
		ikmax = min(min(hkx,hx-im),hkx-imn);
		pkern1 = &( (vi->Kern)  (ikmin,jk));
		pkern2 = &( (vi->Kern)  (ikmin+imn,jk+jmn));
		pw     = &( (vi->OptWeight) (ikmin+im,jk+jm));
		for (ik=ikmin; ik<=ikmax; ++ik, ++pkern1, ++pkern2, ++pw)  
		  summat += (*pkern1)*(*pkern2)*(*pw);
	      }
	      PMat(galstart+m,galstart+n) += summat; 	
	    }
	}
    }
  
#ifdef DEBUG
  cout << " SimFit::fillGalGal is ending ... " << endl;
  cout << " galstart =" << galstart << endl;
  cout << " galend   =" << galend << endl;
  
#endif
  int ngal = galend-galstart+1;
  for (int j=0; j<ngal; ++j) 
    for (int i=j+1; i<ngal; ++i) { 
      MatGal(i,j) = PMat(galstart+i,galstart+j);
    }
  refill = false;
#ifdef DEBUG
  cout << " SimFit::fillGalGal() done " << endl;
  cout << " nvignets in matrix (galgal) = " << count << endl; 
#endif
}



void SimFit::fillGalSky()
{
  
  //**********************
  // gal-sky matrix terms
  //*********************
  
#ifdef FNAME
  cout << " > SimFit::fillGalSky()" << endl;
#endif

  // loop over vignets
  DPixel *pw, *pkern;
  int skyind = 0;
  double summat;

  for (SimFitVignetCIterator it = begin(); it != end(); ++it) {
    
    const SimFitVignet *vi = *it;
    if(! vi->FitSky )
      continue;
    if(! vi->UseGal ) {
      ++skyind;
      continue;
    }

    
    int hx = vi->Hx();
    int hy = vi->Hy();
    
    int ind = skystart+skyind;
    
    // dirac case
    if (vi->DontConvolve) {
      for (int j=-hy;  j<=hy; ++j) {
	pw = &(vi->OptWeight)(-hx,j);
	for (int i=-hx; i<=hx; ++i) {
	  // PMat(galind(i,j),ind) = *pw++; // that's wrong
	  PMat(ind,galind(i,j)) = *pw++; // ok cause ind > galind
	}
      }
      ++skyind;
      continue;
    }
	
    // loop over smallest area among fitting size and weight
    int hkx = vi->Kern.HSizeX();
    int hky = vi->Kern.HSizeY();
    int hsx = (hx + hkx) > hfx ? hfx : (hx + hkx);
    int hsy = (hy + hky) > hfy ? hfy : (hy + hky);
    
    for (int is=-hsx;  is<=hsx; ++is) {
      KERNIND(hkx,hx,is,ikstart,ikend);
      int ikstartis = ikstart-is;
      for (int js=-hsy;  js<=hsy; ++js) {
	KERNIND(hky,hy,js,jkstart,jkend);
	summat = 0.;
	// sum over kernel
	for (int jk=jkstart; jk<=jkend; ++jk) {
	  pkern = &(vi->Kern)  (ikstartis,jk-js);
	  pw    = &(vi->OptWeight)(ikstart,jk);
	  for (int ik=ikstart; ik<=ikend; ++ik) {
	    summat += (*pkern) * (*pw);
	    ++pkern; ++pw;
	  }
	}
	// now fill matrix elements
	// PMat(galind(is,js),ind) = summat; // that's wrong
	PMat(ind,galind(is,js)) = summat;  // ok cause ind > galind
      }
    }
    ++skyind;
  }     
}

void SimFit::fillSkySky()
{
  //*********************************************
  // sky-sky matrix terms and flux vector terms
  //*********************************************

#ifdef FNAME
  cout << " > SimFit::fillSkySky()" << endl;
#endif

  // loop over vignets
  int skyind = 0;
  double sumvec, summat;
  DPixel *pres, *pw;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if(vi->FitSky) {

      // sum over vignet
      int hx = vi->Hx();
      int hy = vi->Hy();
      sumvec = 0.;
      summat = 0.;
      for (int j=-hy; j<=hy; ++j)
	{
	  pw   = &(vi->OptWeight)(-hx,j);
	  pres = &(vi->Resid) (-hx,j);
	  for (int i=-hx; i<=hx; ++i) 
	    {
	      sumvec += (*pw) * (*pres);
	      summat += (*pw);
	      ++pres; ++pw;
	    }
	}

      // now fill out the matrix and vector
      int ind = skystart+skyind;
      Vec(ind) = sumvec;
      PMat(ind,ind) = summat;
      ++skyind;
      }
    }
  //cout << "in fillSkySky() skyind = " << skyind << endl;
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ::::::::::::::::::    Solving routines   ::::::::::::::::::::::::::::
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/

double SimFit::computeChi2() const
{
#ifdef FNAME
  cout << " > SimFit::computeChi2()" << endl;
#endif
  double c = 0.;
  int count = 0;
  double ic = 0;
  for (SimFitVignetCIterator it=begin(); it != end(); ++it) {
    if((*it)->CanFitFlux && dont_use_vignets_with_star)
      continue;
    if(! (*it)->CanFitFlux && ! (*it)->UseGal && ! (*it)->CanFitSky && ! (*it)->CanFitPos  )
      continue;
    
    ic = (*it)->Chi2();
    c += ic;
    count ++;
  }
  //cout << "   in SimFit::computeChi2() number of vignets = " << count << endl;
  return c;
}


bool SimFit::CheckConsistency() const {
  
  for (SimFitVignetCIterator it=begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      if( vi->Star->x != VignetRef->Star->x || vi->Star->y != VignetRef->Star->y ) {
	DumpAndAbort("Inconsistent coordinates of star between vignets");
	return false;
      }
    }
  return true;
}

void SimFit::DumpAndAbort(const string& message) const {
  cout << "DumpAndAbort cause: " << message << endl;
  cout << "===============================================" << endl;
  cout << *VignetRef << endl;
  for (SimFitVignetCIterator it=begin(); it != end(); ++it) {
    (*it)->DumpDebug();
  }
}

bool SimFit::Update(double Factor, bool print)
{
  if(fit_pos) {
    float maxoffset=1;  // we don't want to get offsets larger than one pixel
    if(fabs(Vec(xind)*Factor)>maxoffset) {
      return Update(0.9*maxoffset/fabs(Vec(xind)));
    }
    if(fabs(Vec(yind)*Factor)>maxoffset) {
      return Update(0.9*maxoffset/fabs(Vec(yind)));
    }
  }

  ios::fmtflags oldflags = cout.flags();
  if (print) cout.setf(ios::fixed);

  // update the galaxy
  if (fit_gal) {
    double gf=0.;
    if (print) 
      gf = GalFlux();
    for (int j=-hfy; j<=hfy; ++j) 
      for (int i=-hfx; i<=hfx; ++i)
	VignetRef->Galaxy(i,j) += Vec(galind(i,j))*Factor;
    if (print) {
      cout << " [gal "  << setprecision(3) << setw(6) 
	   << (GalFlux() - gf)/gf*100   << "%]";
    }
  }

  // update the reference position and psf
  if (fit_pos) {
    if (!(VignetRef->ShiftCenter(Point(Vec(xind)*Factor, Vec(yind)*Factor)))) {
      FatalError("in Update, ShiftCenter failure");
      return false;
    }
    VignetRef->UpdatePsfResid();
    if (print) {
      double dx = Vec(xind) * Factor;
      double dy = Vec(yind) * Factor;
      double dpos = sqrt(dx*dx + dy*dy);
      cout << " [pos "  << setprecision(4) << setw(6) 
	   << dpos << "px]";
    }
  }

  // update all vignets
  int fluxind = fluxstart;
  int skyind  = skystart;

  double tf=0., ts=0.;
  if (print) {
    if (fit_flux) tf = TotFlux();
    if (fit_sky)  ts = TotSky();
  }

  for (SimFitVignetIterator it=begin(); it != end(); ++it) {
    SimFitVignet *vi = *it;
    
    // update flux
    if ((fit_flux) && (vi->FitFlux))
      vi->Star->flux += Vec(fluxind++)*Factor;
    
    // update pos (not really required if not vignetref)
    if ((fit_pos) && !(vi->ShiftCenter(Point(Vec(xind)*Factor, Vec(yind)*Factor)))) {
      FatalError("in Update, ShiftCenter failure");
      return false;
    }
    
    // update sky
    if (fit_sky && (vi->FitSky))
      vi->Star->sky += Vec(skyind++)*Factor;
    
    // update residuals and convolved psf
    vi->ModifiedResid();
    vi->Update();
  }

  if (print) {
    if (fit_flux) 
      cout << " [flux " << setprecision(3) << setw(6) 
	   << (TotFlux()-tf)/tf*100 << "%]";      
    if (fit_sky)
      cout << " [sky "  << setprecision(3) << setw(6) 
	   << (TotSky()-ts)/ts*100   << "%]";
    
    cout.flags(oldflags);
    cout << endl;
  }

  return true;
}

void SimFit::FatalError(const char* comment) 
{
  cerr << " > SimFit::FatalError() : " << comment << endl;
  fatalerror=true;
  // dump vignet info
  cerr << " > SimFit::FatalError() : ref star "
       << " [x=" << VignetRef->Star->x
       << "  y=" << VignetRef->Star->y
       << "  flux= " << VignetRef->Star->flux
       << "]\n ";

  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    {
      const SimFitVignet *vi = *it;
      cerr << " > SimFit::FatalError() : " << fixed << right
	   << setw(10) << vi->Image()->Name() 
	   << " mjd " << setprecision(2) << setw(7) << vi->ModifiedJulianDate()
	   << " [flux=" << setprecision(1) << setw(7) << vi->Star->flux << " fit=" << boolalpha << vi->CanFitFlux
	   << "] [pos=(" << setprecision(1) << setw(6) << vi->Star->x 
	   << ", " << setprecision(1) << setw(6) << vi->Star->y << ") fit=" << boolalpha << vi->CanFitPos
	   << "] [sky= " << setprecision(1) << setw(7) << vi->Star->sky << " fit=" << boolalpha << vi->CanFitSky
	   << "]\n";
    }

  // set 0 to all fluxes
  for (SimFitVignetIterator it=begin(); it != end(); ++it) {
    SimFitVignet *vi = *it;
    vi->Star->sky = 0;
    vi->Star->flux = 0;
  }
}


double SimFit::oneNRIteration(double OldChi2)
{
#ifdef FNAME
  cout << " > SimFit::oneNRIteration()" << endl;
#endif
  cout << " > SimFit::oneNRIteration() : Filling  ...";
  FillMatAndVec();
  if(PMat.SizeX()==0) {
    FatalError(" > SimFit::oneNRIteration() Error : NULL matrix");
    return -12;
    // try to exit without abort
  }

  
  /*
  cout << "je verifie que la matrice est triangulaire" << endl;
  for(int i=0;i<PMat.SizeX();i++) {
    if( PMat(i,i)==0 ) {
      cerr << "arggg null diagonal element at i = " << i << endl;
      abort();
    }
    for(int j=0;j<i;j++) {
      if ( PMat(i,j)*PMat(i,j)>PMat(i,i)*PMat(j,j) ) {
	cerr << "arggg pb at i,j = " << i << " " << j << " " << endl;
      }
      if ( PMat(j,i) != 0 ) {
	cerr << "not null val at i,j = " << i << " " << j << " de l'autre cote" << PMat(i,j) << endl;
      }
    }
  }
  */
  
  Mat PMatcopy = PMat;
  Vect Vectcopy = Vec;
  cout << "\r" << flush << " > SimFit::oneNRIteration() : Solving  ...";
  if (cholesky_solve(PMat,Vec,"L") != 0) {
    cout << flush << endl;
    FatalError("\n > SimFit::oneNRIteration() Error : cholesky_solve failure");
    if(false) {
      FatalError("in oneNRIteration, cholesky_solve failure");
      cout << "writing DEBUG_pmat.{fits,mat} and weight vignets before exit ... " << endl;
      PMatcopy.writeFits("DEBUG_pmat.fits");
      PMatcopy.writeASCII("DEBUG_pmat.dat");
      write("sn","./", WriteWeight);
      return -12;
    }else{
      float scaling = 0.995;
      cout << "Error with cholesky, assuming it is unlikely to be a bug," << endl;
      cout << "decrease by " << scaling << " no diagonal values" << endl;
      PMat = PMatcopy;
      Vec = Vectcopy;
      for(unsigned int i=0;i<PMat.SizeX();i++)
	for(unsigned int j=0;j<i;j++)
	  PMat(i,j)*=scaling;
      if(cholesky_solve(PMat,Vec,"L")!=0) {
	FatalError("in oneNRIteration, cholesky_solve failure (after a try to fix matrix)");
	cout << "writing DEBUG_pmat.{fits,mat} and weight vignets before exit ... " << endl;
	PMatcopy.writeFits("DEBUG_pmat.fits");
	PMatcopy.writeASCII("DEBUG_pmat.dat");
	write("sn","./", WriteWeight);
	return -12;
      }
    }
  }
  cout << "\r" << flush << " > SimFit::oneNRIteration() : Updating ";
  Update();
  
  if (fatalerror) return -12;

  double curChi2 = computeChi2();
 
  if (curChi2-OldChi2>0.01) {
#ifdef DEBUG
    cout << " SimFit::oneNRIteration(" 
	 << OldChi2 << "):  chi2=" << curChi2 
	 << " increased (diff=" << curChi2-OldChi2 <<"), reducing corrections" << endl;
#endif     
     
    double dof  = ndata-nparams;
    double fact = 1.;
    double step = 0.1;
    while ((curChi2 > OldChi2) && (fact > step*1.5)) {
      Update(-fact, false); if(fatalerror) return -12;
      fact -= step;
      Update(fact, false);  if(fatalerror) return -12;
      curChi2 = computeChi2();
#ifdef DEBUG
      cout << " > SimFit::oneNRIteration():  curChi2=" 
	   << curChi2/dof << " diff/dof = " << (curChi2-OldChi2)/dof
	   << " fact = " << fact << endl;
#endif     
    }
    
    // reducing corrections had no effect 
    if (curChi2 > OldChi2) {
      Update(-fact, false); if(fatalerror) return -12;
      curChi2 = computeChi2();
      cout << " > SimFit::oneNRIteration(): diff/dof = " 
	   << setprecision(5) << (curChi2-OldChi2)/dof << " fact = " << fact
	   << " return to beginning" << endl;
    } else 
      cout << " > SimFit::oneNRIteration():  curChi2/dof =" 
	   << setprecision(5) << curChi2/dof << " fact = " << fact << endl;      
  }

  return curChi2;
}

bool SimFit::IterateAndSolve(const int MaxIter,  double Eps)
{
  double oldchi2;
  chi2 = computeChi2();
  if(fatalerror) {
    FatalError("IterateAndSolve after computeChi2");
    return false;
  }
  int iter = 0;
  double diff = 1000;
  cout << " > SimFit::IterateAndSolve() : Initialize  chi2/dof = " << chi2/(ndata-nparams) 
       << " [MaxIter=" << MaxIter << " Eps="
       << Eps << "]" << endl; 
  bool redoweight = false;
  do
    {
      oldchi2 = chi2;
      chi2 = oneNRIteration(oldchi2);
      if(fatalerror) {
	FatalError("IterateAndSolve after oneNRIteration"); 
	return false;
      }
      diff = fabs(chi2-oldchi2);
      cout << " > SimFit::IterateAndSolve(): Iteration " << setw(2) << iter << " chi2/dof = " 
	   << setprecision(5) << chi2/(ndata - nparams)
	   << " dchi2 = " << chi2-oldchi2
	   << endl;
      
      if (redoweight && diff < Eps*5.) {
	chi2 += Eps*10.;
	diff += Eps*10.;
	for (SimFitVignetIterator it=begin(); it != end(); ++it)
	  (*it)->RedoWeight();
	redoweight = false;
      }
    } while ((iter++ < MaxIter) && (diff>Eps));
  
  cout << " > SimFit::IterateAndSolve(): Done" << endl;
  return true;
} 


bool SimFit::GetCovariance()
{
  //#ifdef FNAME
  cout << " > SimFit::GetCovariance()" << endl;
  //#endif
  if (!inverted) {
    int status = cholesky_invert(PMat,"L");
    if (status != 0) {
      cerr << " SimFit::GetCovariance() : Error: inverting failed. Lapack status=: " 
	   << status << endl;
      FatalError(" in GetCovariance, inverting failed");
      return false;
    }
    inverted = true;
  }

  // rescale covariance matrix with estimated global sigma scale factor
  // it corrects for initially under-estimated (ex: correlated) weights if chi2/dof < 1 
  // or for error in our model if chi2/dof > 1

  double sigscale = chi2 / (ndata - nparams);
  //double sigscale = 1;
  if(sigscale<1) sigscale = 1;
  
  int fluxind = 0;
  int skyind  = skystart;
  for (SimFitVignetIterator it=begin(); it != end(); ++it)
    {
      if ((fit_flux) && (*it)->FitFlux) {
	// (*it)->Star->varflux = sigscale * PMat(fluxind,fluxind++); // DO NOT USE THIS
	(*it)->Star->eflux = sqrt(sigscale * PMat(fluxind,fluxind));
	(*it)->Star->sigscale_varflux = sigscale ;
	fluxind++;
      }
      if (fit_sky && (*it)->FitSky) {  
	// (*it)->Star->varsky  = sigscale * PMat(skyind,skyind++); // DO NOT USE THIS
	(*it)->Star->varsky  = sigscale * PMat(skyind,skyind);
	skyind++;
      }
    }

  if (fit_pos)
    {
      VignetRef->Star->vx  = sigscale * PMat(xind,xind);
      VignetRef->Star->vy  = sigscale * PMat(yind,yind);
      VignetRef->Star->vxy = sigscale * PMat(xind,yind);
      for (SimFitVignetIterator it=begin(); it != end(); ++it) {
	(*it)->Star->vx  = sigscale * PMat(xind,xind);
	(*it)->Star->vy  = sigscale * PMat(yind,yind);
	(*it)->Star->vxy = sigscale * PMat(xind,yind);
      }
    }
  
  return true;
}

void SimFit::FitInitialGalaxy() {
#ifdef FNAME
  cout << " > SimFit::FitInitialGalaxy()" << endl;
#endif
  int nzero=0;
  for (SimFitVignetIterator it=begin(); it != end(); ++it) {
    if( (!((*it)->CanFitFlux)) && (*it)->CanFitGal )
      nzero ++;
  }
  cout << " > SimFit::FitInitialGalaxy() with " << nzero << " exposures" << endl;  
  if(nzero==0) {
    FatalError(" in FitInitialGalaxy, nzero=0");
    abort();
  }
  SetWhatToFit(FitGal);
  dont_use_vignets_with_star = true;
  Resize(1);
  IterateAndSolve(2,0.1);
}


bool SimFit::DoTheFit(int MaxIter, double epsilon)
{
#ifdef FNAME
  cout << " > SimFit::DoTheFit() : Starting" << endl;  
#endif
  if (!(fit_pos || fit_sky || fit_flux || fit_gal)) 
    {
      cerr << " > SimFit::DoTheFit() Error : nothing to fit" << endl;
      return false;
    }

#ifdef DEBUG
  clock_t tstart = clock(); 
#endif 

  Resize(1);
  if (!IterateAndSolve(MaxIter,epsilon) || !GetCovariance()) {
    FatalError("DoTheFit set all fluxes to zero before quitting because of failure");
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      (*it)->Star->flux = 0;
      (*it)->Star->eflux = 0;
    }    
    return false;
  }
  
#ifdef DEBUG
  clock_t tend = clock();
  cout << " > SimFit::DoTheFit() : CPU comsumed : " 
       <<  float(tend-tstart) / float(CLOCKS_PER_SEC) << endl;
#endif
  return true;
}

void SimFit::write(const string& StarName,const string &DirName, unsigned int whattowrite) 
{
#ifdef FNAME
  cout << " > SimFit::write(" << StarName << ")" << endl;
#endif
  cout << " > SimFit::write()  ";
  // write lc
  if(whattowrite & WriteLightCurve) {
    ofstream lstream(string(DirName+"/lightcurve_"+StarName+".dat").c_str());
    lstream << setiosflags(ios::fixed);
    front()->Star->WriteHeader(lstream);
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      (*it)->Star->writen(lstream);
      lstream << endl;
    }
    lstream.close();
  }
  
  // write vignets info
  if(whattowrite & WriteVignetsInfo) {
    ofstream vignetstream(string(DirName+"/vignets_"+StarName+".dat").c_str());
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      
      if((*it)->CanFitFlux && dont_use_vignets_with_star)
	continue;
      if(! (*it)->CanFitFlux && ! (*it)->UseGal && ! (*it)->CanFitSky && ! (*it)->CanFitPos  )
	continue;
      SimFitVignet *vi = *it;
      vignetstream << *vi << endl;
      int npix=0;
      float chi2 = vi->CentralChi2(npix);
      vignetstream << "centralchi2 npix chi2pdf " << chi2 << " " << npix << " " << chi2/(npix-1) << endl;
    }
    vignetstream.close();
  }
  
  // write galaxy fits
  if(whattowrite & WriteGalaxy)
    VignetRef->Galaxy.writeFits(DirName+"/galaxy_"+StarName+".fits");
  

  // write residuals
  if(whattowrite & WriteResid){
    cout << " resids ";
    for (SimFitVignetIterator it=begin(); it != end() ; ++it){      
      (*it)->ClearResidZeroWeight();
      (*it)->Resid.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_resid.fits");
    }
  }
  
  // write Data
  if(whattowrite & WriteData){
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->Data.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_data.fits");
  }
  
  // write Psf
  if(whattowrite & WritePsf){
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      (*it)->Psf.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_psf.fits");
      //(*it)->Psf.writeGaussian(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_gaus.fits");
    }
  }
  // write Kernel
  if(whattowrite & WriteKernel){
    for (SimFitVignetIterator it=begin(); it != end() ; ++it) {
      (*it)->Kern.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_kern.fits");
    }
  }
  // write weights
  if(whattowrite & WriteWeight){
    for (SimFitVignetIterator it=begin(); it != end() ; ++it)
      (*it)->OptWeight.writeFits(DirName+"/"+(*it)->Image()->Name()+"_"+StarName+"_weight.fits");
  }
  //write matrices
  if(whattowrite & WriteMatrices) {
    PMat.writeFits(DirName+"/pmat_"+StarName+".fits");
  }  
  int i=0;
  for (SimFitVignetIterator it=begin(); it != end() ; ++it)
    if((*it)->FitFlux)
      i++;
  Mat vm(1,i);
  i=0;
  for (SimFitVignetIterator it=begin(); it != end() ; ++it)
    if((*it)->FitFlux) {
      vm(0,i)= (*it)->Star->flux;
      i++;
    }
  vm.writeFits(DirName+"/vec_"+StarName+".fits");
    
  // create matrix of Nights and Images
  fillNightMat();
  NightMat.writeFits(DirName+"/nightmat_"+StarName +".fits");

  cout << endl;
}

double SimFit::VarScale() const 
{
  double varscale = Chi2() / Dof();
  return varscale < 1 ? 1. : varscale;
}

double SimFit::VarTotFlux() const
{
  if (!fit_flux) return 0.;
  if (!inverted) {
    cerr << " ERROR : matrix not inverted " << endl;
    return -1;
  }

  double vartotflux = 0.;
  int nflux = fluxend - fluxstart + 1;

  for (int j=0; j<nflux; ++j) {
    vartotflux += PMat(fluxstart+j, fluxstart+j);
    for (int i=j+1; i<nflux; ++i)
      vartotflux += 2. * PMat(fluxstart+i, fluxstart+j);
  }

  return VarScale() * vartotflux;
} 

  
double SimFit::GalFlux() const
{
  if (!fit_gal && use_gal) return 0.;
  return VignetRef->Galaxy.sum();
}

double SimFit::VarGalFlux() const
{
  if (!fit_gal && use_gal) return 0.;
  if (!inverted) {
    cerr << " ERROR : matrix not inverted " << endl;
    return -1;
  }
  double vargalflux = 0.;
  int ngal = galend - galstart + 1;
  for (int j=0; j<ngal; ++j) {
    vargalflux += PMat(galstart+j, galstart+j);
    for (int i=j+1; i<ngal; ++i) { 
      vargalflux += 2. * PMat(galstart+i, galstart+j);
    }
  }
 
  return VarScale() * vargalflux;
}

double SimFit::TotFlux() const 
{
  if (!fit_flux) return 0.;
  double totflux = 0.;
  for (SimFitVignetCIterator it = begin(); it != end(); ++it)
    if ((*it)->FitFlux)
      totflux += (*it)->Star->flux;
  return totflux;
}

double SimFit::TotSky() const 
{
  if (!fit_sky) return 0.;
  double totsky = 0.;
  for (SimFitVignetCIterator itVig = begin(); itVig != end(); ++itVig)
    if ((*itVig)->FitSky) 
      //totsky  += (*itVig)->Star->sky + (*itVig)->Image()->OriginalSkyLevel();
      totsky  += (*itVig)->Star->sky;
  return totsky;
}

double SimFit::VarTotSky() const
{
  if (!fit_sky) return 0.;

  if (!inverted) {
    cerr << " ERROR : matrix not inverted " << endl;
    return -1;
  }

  double vartotsky = 0.;
  int nsky = skyend - skystart + 1;
  for (int j=0; j<nsky; ++j) {
    vartotsky += PMat(skystart+j, skystart+j);
    for (int i=j+1; i<nsky; ++i)
      vartotsky += 2. * PMat(skystart+i, skystart+j);
  }
   
  return VarScale() * vartotsky;

}

double mymedian(vector<double>& array)
{
  sort(array.begin(), array.end());
  size_t n = array.size();
  return n&1? array[n/2] : (array[n/2-1] + array[n/2])*0.5;
}


double SimFit::MeanMedianRms(double& med, double& rms, double& adev) const 
{
  double mean = 0.;
  size_t  nim = 0;
  for (SimFitVignetCIterator it=begin(); it != end() ; ++it) {
    mean += (*it)->Resid.sum();
    nim  += (*it)->Resid.Nx() * (*it)->Resid.Ny();
  }
  mean /= nim;

  double var = 0.;
  adev = 0.;
  vector<double> resid;
  resid.resize(nim);
  vector<double>::iterator itr = resid.begin();

  for (SimFitVignetCIterator it=begin(); it != end() ; ++it) {
    DPixel *pres = (*it)->Resid.begin();
    size_t n  = (*it)->Resid.Nx() * (*it)->Resid.Ny();
    for (size_t i=0; i<n ; ++i, ++itr, ++pres) {
      double delta = *pres - mean;
      var += delta*delta;
      adev += fabs(delta);
      *itr = *pres;
    }
  }

  var /= (nim-1);
  rms = sqrt(var);
  adev /= nim;
  med = mymedian(resid);

  return mean;
}

void SimFit::fillNightMat(LightCurve& Lc) {
#ifdef FNAME
  cout << " > SimFit::fillNightMat" << endl;
#endif
  
  
      
  LightCurve::const_iterator itLc = Lc.begin();
  LightCurve::const_iterator endLc = Lc.end();
  double jd; // julian day
  vector<double> nightdates;
  bool isinnight;
  double timediff = 10./24.; // 10 hours
  int nimages = 0;
  for (;itLc != endLc; ++itLc) {
    jd = (*itLc)->ModifiedJulianDate();
    if(Lc.Ref->IsVariable(jd)) {
      nimages++;
      isinnight=false;
      for(size_t day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  isinnight = true;
	  break;
	}
      }
      if(!isinnight) {
	nightdates.push_back(jd);
      }
    }
  } 
  int nnights = nightdates.size(); // number of nights for this lightcurve
  NightMat.allocate(nnights,nimages); // szie of matrix
  // now we fill this matrix
    //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  itLc = Lc.begin();
  int im = 0;
  for (;itLc != endLc; ++itLc) {
    jd = (*itLc)->ModifiedJulianDate();
    if(Lc.Ref->IsVariable(jd)) {
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  NightMat(day,im)=1;
	  break;
	}
      }
      im++;
    }
  }
  
  cout << "========= NightMat ==========" << endl;
  cout << NightMat << endl;
  return;
}


void SimFit::fillNightMat() { // using simfitvignet
#ifdef FNAME
  cout << " > SimFit::fillNightMat" << endl;
#endif
  
  double jd; // julian day
  vector<double> nightdates;
  bool isinnight;
  double timediff = 10./24.; // 10 hours
  int nimages = 0;
  
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    if((*itVig)->FitFlux) {
      jd = (*itVig)->ModifiedJulianDate();
      nimages++;
      isinnight=false;
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  isinnight = true;
	  break;
	}
      }
      if(!isinnight) {
	nightdates.push_back(jd);
      }
    }
  } 
  int nnights = nightdates.size(); // number of nights for this lightcurve
  NightMat.allocate(nnights,nimages); // szie of matrix
  // now we fill this matrix
    //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  
  int im = 0;
  for (SimFitVignetIterator itVig = begin(); itVig != end(); ++itVig) {
    if((*itVig)->FitFlux) {
      jd = (*itVig)->ModifiedJulianDate();
      for(unsigned int day=0;day<nightdates.size(); ++day) {
	if(fabs(jd-nightdates[day])<timediff) {
	  NightMat(day,im)=1;
	  break;
	}
      }
      im++;
    }
  }
  //cout << "========= NightMat ==========" << endl;
  //cout << NightMat << endl;
  return;
}


ostream& operator << (ostream& Stream, const CountedRef<SimFitVignet> &Vig)
{
  Stream << *Vig;
  return Stream;
}


ostream& operator << (ostream& Stream, const SimFit &MyFit)
{
  return Stream;
}


ostream& SimFit::DumpMatrices(ostream& Stream) const {
  Stream << "====  SimFit::DumpMatrices ====" << endl;
  Stream  << endl;
  Stream  << "Vec" << endl;
  Stream  << Vec << endl;
  Stream  << "PMat" << endl;
  Stream  << PMat << endl;
  return Stream;

}

void SimFit::Chi2PosMap() {
  cout << "Chi2PosMap Start" << endl;
  int hx = 15;
  int hy = 15;
  float pixstep = 0.02;
  float fittedposx = VignetRef->Star->x;
  float fittedposy = VignetRef->Star->y;
  float posx,posy;
  Mat PosMap(2*hx+1,2*hy+1);
  for(int j=-hy;j<=hy;j++) {     
    for(int i=-hx;i<=hx;i++) { 
      cout << "Chi2PosMap " << i << " " << j << endl;
      posx = fittedposx+i*pixstep;
      posy = fittedposy+j*pixstep;
      VignetRef->Star->x = posx;
      VignetRef->Star->y = posy;
      VignetRef->UpdatePsfResid();
      for (SimFitVignetIterator it=begin(); it != end(); ++it) {
	SimFitVignet *vi = *it;
	vi->Star->x = posx;
	vi->Star->y = posy;
	vi->ModifiedResid();
	vi->Update();
      }
      PosMap(i + hx,j + hy) = computeChi2();
    }
  }
  PosMap.writeFits("Chi2PosMap.fits");
}

// could not make it work. was designed to replace the ugly test in the fillGalGal loop above
#if 0
DPixel *pkerm, *pkern;
#define MAX3(A,B,C) (A>B ? (A>C ? A : C) : (B>C ? B : C))
#define MIN3(A,B,C) (A<B ? (A<C ? A : C) : (B<C ? B : C))

	  int ikstart = MAX3(-hkx, -hkx-imn, -hsx-im);
	  int jkstart = MAX3(-hky, -hky-jmn, -hsy-jm);
	  int ikend = MIN3(hkx, hkx-imn, hsx-im);
	  int jkend = MIN3(hky, hky-jmn, hsy-jm);
          ikstart = (ikstart < ikend) ? ikstart : ikend;  
          jkstart = (jkstart < jkend) ? jkstart : jkend;  
	  summat = 0.;
	  for (int jk=jkstart; jk<=jkend; ++jk)
	    {
	      pkerm = &(vi->Kern)  (ikstart,jk);
	      pkern = &(vi->Kern)  (ikstart+imn,jk+jmn);
	      pw    = &(vi->OptWeight)(ikstart+im,jk+jm);
	      for (int ik=ikstart; ik<=ikend; ++ik)
		{
		  summat += (*pkerm) * (*pkern) * (*pw); 
		  ++pkerm; ++pkern; ++pw;
		}
	    }
	  Mat[(galstart+m)*nparams+galstart+n] = summat; 	
	}
    }
}
#endif

