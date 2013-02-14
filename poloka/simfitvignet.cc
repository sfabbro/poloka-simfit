#include <fstream> 
#include <poloka/simfitvignet.h>
#include <poloka/vignetserver.h>
#include <poloka/kernelfitter.h>


#define DEBUG_KERNEL

// uncomment this to include model in weigths
//#define VALCUTOFF 0. 

// uncomment this to normalize  psf, warning: need to modify those functions to do this:
// SimFitRefVignet::UpdatePsfResid
// SimFitVignet::UpdateResid_psf_gal
// SimFitVignet::UpdateResid_psf
// 
#define NORMALIZE_PSF 

static double sq(double x) {return x*x;};

////////////////////////////////////////////////////////////////////////////////////
//  TabulatedPsf
////////////////////////////////////////////////////////////////////////////////////

TabulatedPsf::TabulatedPsf(const Point& Pt, const ImagePSF& imagepsf, const Window& Rect)
  : Kernel(Rect.Hx(), Rect.Hy()), Dx(hSizeX, hSizeY), Dy(hSizeX, hSizeY)
{ 
  Tabulate(Pt, imagepsf, Rect);
}

void TabulatedPsf::Resize(const int Hx, const int Hy)
{
#ifdef FNAME
  cout << " > TabulatedDaoPsf::Resize(const int Hx, const int Hy)" << endl;
#endif

#ifdef DEBUG
 cout << "before Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;
#endif
  if (HSizeX() != Dx.HSizeX() ||
      HSizeY() != Dx.HSizeY() ||
      HSizeX() != Dy.HSizeX() ||
      HSizeY() != Dy.HSizeY() ||
      HSizeX() != Hx ||
      HSizeX() != Hy)
    {
      Allocate(2*Hx+1, 2*Hy+1);
      Dx.Allocate(Nx(), Ny());
      Dy.Allocate(Nx(), Ny());
    }
#ifdef DEBUG
  cout << "after Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;
#endif
}

void TabulatedPsf::Tabulate(const Point& Pt, const ImagePSF& imagepsf, const Window& Rect)
{
  Resize(Rect.Hx(), Rect.Hy());
  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();
  Vect der(2);
  integral = 0;
  for (int j=Rect.ystart; j<Rect.yend; ++j)
    for (int i=Rect.xstart; i<Rect.xend; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = imagepsf.PSFValue(Pt.x,Pt.y,i,j,&der);
	*ppdx = der(0);
	*ppdy = der(1);
	integral += *ppsf;
      }
  
  // normalization
#ifdef NORMALIZE_PSF 
  double norme = 1./integral;
  ppsf = begin();
  ppdx = Dx.begin();
  ppdy = Dy.begin();
  
  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf *= norme;
	*ppdx *= norme;
	*ppdy *= norme;
      }
  integral=1;
#endif

  

}

void TabulatedPsf::ComputeMoments() {
  integral = 0;
  DPixel *ppsf = begin();
  for (int j=-hSizeY; j<=hSizeY; ++j) {
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf) {
      integral += (*ppsf);
    }
  }
  

  double wxx=0.5;
  double wyy=0.5;
  double wxy=0;
  int iter=0;
  mx = 0;
  my = 0;
  det=0;
  while (iter < 10)
    {
      iter++;
      
      double sumx = 0;
      double sumy = 0;
      double sumxx = 0;
      double sumyy = 0;
      double sumxy = 0;
      double sumw = 0;
      
      DPixel *ppsf = begin();
      for (int j=-hSizeY; j<=hSizeY; ++j) {
	for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf) {
	  double dx = i-mx;
	  double dy = j-my;
	  double wg = wxx*dx*dx + wyy*dy*dy + 2.*wxy*dx*dy;
	  if (wg > 16) continue; // 4 sigmas, and avoids overflows
	  wg = exp(-0.5*wg);
	  wg *= (*ppsf);
	  sumx += wg*dx;
	  sumy += wg*dy;
	  sumxx += wg*dx*dx;
	  sumyy += wg*dy*dy;
	  sumxy += wg*dx*dy;
	  sumw += wg;
	}
      }
      
      sumx /= sumw;
      sumy /= sumw;
      sumxx /= sumw;
      sumyy /= sumw;
      sumxy /= sumw;
      
      mx += sumx;
      my += sumy;
      sumxx -= sq(sumx);
      sumyy -= sq(sumy);
      sumxy -= sumx*sumy;
      
      det = sumxx*sumyy - sq(sumxy);
      wxx = 0.5*sumyy/det;
      wyy = 0.5*sumxx/det;
      wxy = -0.5*sumxy/det;
    } // end of iteration.
  
  det = wxx*wyy-sq(wxy);
  mxx = wyy/det;
  myy = wxx/det;
  mxy = -wxy/det;
  det = 1./det;
  
  
}

////////////////////////////////////////////////////////////////////////////////////
//  TabulatedDaoPsf
////////////////////////////////////////////////////////////////////////////////////
#ifdef STORAGE
TabulatedDaoPsf::TabulatedDaoPsf(const Point& Pt, const DaoPsf& Dao, const int Radius)
  : Kernel(Radius), Dx(Radius), Dy(Radius)
{ 
  Tabulate(Pt, Dao, Radius);
}

TabulatedDaoPsf::TabulatedDaoPsf(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
  : Kernel(Rect.Hx(), Rect.Hy()), Dx(hSizeX, hSizeY), Dy(hSizeX, hSizeY)
{ 
  Tabulate(Pt, Dao, Rect);
}

void TabulatedDaoPsf::Resize(const int Hx, const int Hy)
{
#ifdef FNAME
  cout << " > TabulatedDaoPsf::Resize(const int Hx, const int Hy)" << endl;
#endif

#ifdef DEBUG
 cout << "before Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;
#endif
  if (HSizeX() != Dx.HSizeX() ||
      HSizeY() != Dx.HSizeY() ||
      HSizeX() != Dy.HSizeX() ||
      HSizeY() != Dy.HSizeY() ||
      HSizeX() != Hx ||
      HSizeX() != Hy)
    {
      Allocate(2*Hx+1, 2*Hy+1);
      Dx.Allocate(Nx(), Ny());
      Dy.Allocate(Nx(), Ny());
    }
#ifdef DEBUG
  cout << "after Hx Hy HSizeX() HSizeY() Nx() Ny() hSizeX hSizeY " 
       << Hx << " "
       << Hy << " "
       << HSizeX() << " "
       << HSizeY() << " "
       << Nx() << " "
       << Ny() << " "
       << hSizeX << " "
       << hSizeY << endl;
#endif
}

void TabulatedDaoPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius)
{
#ifdef FNAME
  cout << " > TabulatedDaoPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius)" << endl;
#endif

  Resize(Radius, Radius);

  const int ic = int(Pt.x);
  const int jc = int(Pt.y);

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();
  

  integral = 0;

  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(ic+i,jc+j, Pt, *ppdx, *ppdy);
	integral += *ppsf;
      }
  
  cout << " PSF_INTEGRAL " << Pt.x << " " << Pt.y << " " << integral << endl;

  if(!(integral>0.5)) {
    cout << "Very strange psf integral in TabulatedDaoPsf::Tabulate (1) : " << integral << endl;
    cout << "hSizeX,hSizeX=" << hSizeX << "," << hSizeY << endl;
    cout << "better quit (first write psf as fits)" << endl;
    writeFits("psf_DEBUG.fits");
    abort();
  }


#ifdef NORMALIZE_PSF
  // normalization
  double norme = 1./integral;
  
  ppsf = begin();
  ppdx = Dx.begin();
  ppdy = Dy.begin();
  
  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf *= norme;
	*ppdx *= norme;
	*ppdy *= norme;
      }

  integral = 1;
#endif
  
  ComputeMoments();
}

void TabulatedDaoPsf::Tabulate(const Point& Pt, const DaoPsf& Dao, const Window& Rect)
{

  Resize(Rect.Hx(), Rect.Hy());

  DPixel *ppsf = begin();
  DPixel *ppdx = Dx.begin();
  DPixel *ppdy = Dy.begin();

  integral = 0;
#ifdef DEBUG
  cout << "in  TabulatedDaoPsf::Tabulate DUMP , Pt.x, Pt.y = " << Pt.x << " " << Pt.y << endl;
#endif
  for (int j=Rect.ystart; j<Rect.yend; ++j)
    for (int i=Rect.xstart; i<Rect.xend; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf = Dao.Value(i,j, Pt, *ppdx, *ppdy);
	integral += *ppsf;
      }
  cout << " PSF_INTEGRAL " << Pt.x << " " << Pt.y << " " << integral << endl;
  if(!(integral>0.5)) {
    cout << "Very strange psf integral in TabulatedDaoPsf::Tabulate (2) : " << integral << endl;
    cout << "Rect=" << Rect << endl;
    cout << "better quit (first write psf as fits)" << endl;
    writeFits("psf_DEBUG.fits");
    abort();
  }

// normalization
#ifdef NORMALIZE_PSF 
  double norme = 1./integral;
  ppsf = begin();
  ppdx = Dx.begin();
  ppdy = Dy.begin();
  
  for (int j=-hSizeY; j<=hSizeY; ++j) 
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf, ++ppdx, ++ppdy) 
      {
	*ppsf *= norme;
	*ppdx *= norme;
	*ppdy *= norme;
      }
  integral=1;
#endif
  
  ComputeMoments();
}


//#define DEBUG_GAUS
void TabulatedDaoPsf::ComputeMoments() {
#ifdef FNAME
  cout << " > TabulatedDaoPsf::ComputeMoments()" << endl;
#endif
#ifdef DEBUG_GAUS
  cout << " > TabulatedDaoPsf::ComputeMoments()" << endl;
#endif

  integral = 0;
  DPixel *ppsf = begin();
  for (int j=-hSizeY; j<=hSizeY; ++j) {
    for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf) {
      integral += (*ppsf);
    }
  }
  

  double wxx=0.5;
  double wyy=0.5;
  double wxy=0;
  int iter=0;
  mx = 0;
  my = 0;
  det=0;
  while (iter < 10)
    {
      iter++;
      
      double sumx = 0;
      double sumy = 0;
      double sumxx = 0;
      double sumyy = 0;
      double sumxy = 0;
      double sumw = 0;
      
      DPixel *ppsf = begin();
      for (int j=-hSizeY; j<=hSizeY; ++j) {
	for (int i=-hSizeX; i<=hSizeX; ++i, ++ppsf) {
	  double dx = i-mx;
	  double dy = j-my;
	  double wg = wxx*dx*dx + wyy*dy*dy + 2.*wxy*dx*dy;
	  if (wg > 16) continue; // 4 sigmas, and avoids overflows
	  wg = exp(-0.5*wg);
	  wg *= (*ppsf);
	  sumx += wg*dx;
	  sumy += wg*dy;
	  sumxx += wg*dx*dx;
	  sumyy += wg*dy*dy;
	  sumxy += wg*dx*dy;
	  sumw += wg;
	}
      }
      
      sumx /= sumw;
      sumy /= sumw;
      sumxx /= sumw;
      sumyy /= sumw;
      sumxy /= sumw;
      
      mx += sumx;
      my += sumy;
      sumxx -= sq(sumx);
      sumyy -= sq(sumy);
      sumxy -= sumx*sumy;
      
      det = sumxx*sumyy - sq(sumxy);
      wxx = 0.5*sumyy/det;
      wyy = 0.5*sumxx/det;
      wxy = -0.5*sumxy/det;
    } // end of iteration.
  
  det = wxx*wyy-sq(wxy);
  mxx = wyy/det;
  myy = wxx/det;
  mxy = -wxy/det;
  det = 1./det;

#ifdef DEBUG_GAUS
  cout << "TabulatedDaoPsf::ComputeMoments: mx,my,mxx,myy,mxy = " 
       << mx << ", "
       << my << ", "
       << mxx << ", "
       << myy << ", "
       << mxy << endl;
  
  cout << "TabulatedDaoPsf::ComputeMoments: det = " << det << endl;
  cout << "TabulatedDaoPsf::ComputeMoments: integral = " << integral << endl;
  cout << "TabulatedDaoPsf::ComputeMoments Gaus,dGausdx2,dGausdy2,dGausdxdx(0,0) = " 
       << Gaus(0,0) << ", "
       << dGausdx2(0,0) << ", "
       << dGausdy2(0,0) << ", "
       << dGausdxdy(0,0)
       << endl;
  cout << "TabulatedDaoPsf::ComputeMoments Gaus,dGausdx2,dGausdy2,dGausdxdx(2,2) = " 
       << Gaus(2,2) << ", "
       << dGausdx2(2,2) << ", "
       << dGausdy2(2,2) << ", "
       << dGausdxdy(2,2)
       << endl;
  
  
    
#endif
}


double TabulatedDaoPsf::Gaus(int i,int j) const {
  //cout << "TabulatedDaoPsf::Gaus " << det << " " << integral << " " << integral/(6.2831853072*sqrt(det)) << endl;

  return integral/(6.2831853072*sqrt(det))*exp(-1./2./det* ( myy*sq(i-mx)+mxx*sq(j-my)-2*mxy*(i-mx)*(j-my) ) ); 
}


double  TabulatedDaoPsf::dGausdx2(int i,int j) const {
  double dGausdx_over_Gaus = ( mxy*(j-my) - myy*(i-mx) )/det;
  return Gaus(i,j)*( -myy/det + sq(dGausdx_over_Gaus) );
}

double  TabulatedDaoPsf::dGausdy2(int i,int j) const {
  double dGausdy_over_Gaus = ( mxy*(i-mx) - mxx*(j-my) )/det;
  return Gaus(i,j)*( -mxx/det + sq(dGausdy_over_Gaus) );
}

double  TabulatedDaoPsf::dGausdxdy(int i,int j) const {
  double dGausdx_over_Gaus = ( mxy*(j-my) - myy*(i-mx) )/det;
  double dGausdy_over_Gaus = ( mxy*(i-mx) - mxx*(j-my) )/det;
  return Gaus(i,j)*( mxy/det + dGausdx_over_Gaus*dGausdy_over_Gaus);
}

void  TabulatedDaoPsf::writeGaussian(const string& filename)  {
  ComputeMoments();
  Kernel tabulatedgaussian = (*this);
  for(int j=-hSizeY;j<=hSizeY;++j)
    for(int i=-hSizeX;i<=hSizeX;++i)
      tabulatedgaussian(i,j) = Gaus(i,j);
  tabulatedgaussian.writeFits(filename);
}
#endif

////////////////////////////////////////////////////////////////////////////////////
//  SimFitRefVignet
////////////////////////////////////////////////////////////////////////////////////
//=========================================================================================
SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim,bool usegal)
  : Vignet(Rim), imagepsf(new ImagePSF(*Rim,false))
{
#ifdef FNAME
  cout << " > SimFitRefVignet::SimFitRefVignet(const ReducedImage *Rim)" << endl;
#endif
  UseGal=usegal;
}

// make sure that psf, vignet and data have proper sizes.
// allocate if needed
void SimFitRefVignet::Resize(const int Hx, const int Hy)
{
#ifdef FNAME
  cout << " > SimFitRefVignet::Resize(const int Hx, const int Hy)" << endl;
#endif
  
  if  (Hx < 0 || Hy < 0) 
    {
      cerr << " SimFitRefVignet::Resize(" << Hx << "," << Hy << ") : impossible \n";
      return;
    } 
  
  // resize Data, Weight, Resid if necessary
  Vignet::Resize(Hx,Hy);
  
  // resize Psf, Psf.Dx, Psf.Dy 
  Psf.Tabulate(*Star,*imagepsf,*this);
  
  // resize Galaxy 
  if(UseGal) {
    if(Galaxy.Nx()!=Data.Nx() || Galaxy.Ny()!=Data.Ny())
      makeInitialGalaxy();
  }
}

void SimFitRefVignet::Load(const PhotStar *Star)
{
#ifdef FNAME
  cout << " > SimFitRefVignet::Load(const PhotStar *Star)" << endl;
#endif
#ifdef VALCUTOFF
  cout << "VALCUTOFF is defined" << endl;
#endif
#ifdef NORMALIZE_PSF 
  cout << "NORMALIZE_PSF is defined" << endl;
#endif



  if (!Star) return;
  Vignet::Load(Star);

  //printf(" in SimFitRefVignet::Load x,y = %10.10g,%10.10g\n",Star->x,Star->y);
  
  if(!imagepsf) imagepsf = new ImagePSF(*rim,false);
  Psf.Tabulate(*Star, *imagepsf, *this);
  if(UseGal) {
    makeInitialGalaxy(); 
  }
  UpdatePsfResid();
}

void SimFitRefVignet::makeInitialGalaxy()
{  
#ifdef FNAME
  cout << " > SimFitRefVignet::makeInitialGalaxy()" << endl;
#endif

  // allocate galaxy as biggest vignet possible, so resizing it is easy.
  Galaxy.Allocate(Data.Nx(), Data.Ny(), 1);

  // initial galaxy: best resolution image - sky - [flux*psf]  
  DPixel *pdat = Data.begin(), *pgal = Galaxy.begin(), *ppsf = Psf.begin();

  for (int i=Nx()*Ny(); i ; --i)
    {
      //*pgal = *pdat - Star->sky - Star->flux * *ppsf;
      *pgal = 0.;
      ++pgal; ++pdat; ++ppsf;
    }
}


void SimFitRefVignet::UpdatePsfResid()
{
#ifdef FNAME
  if(UseGal)
    cout << " > SimFitVignet::UpdatePsfResid() : updating residuals with galaxy" << endl;
  else
    cout << " > SimFitVignet::UpdatePsfResid() : updating residuals without galaxy" << endl;
#endif
  // re-allocate Resid if too small

  Psf.Tabulate(*Star, *imagepsf, *this);
  
  DPixel *pdat = Data.begin(), *pres = Resid.begin();
  DPixel *ppsf = Psf.begin(), *ppdx = Psf.Dx.begin(), *ppdy = Psf.Dy.begin();


  if(UseGal) {
    DPixel  *pgal = Galaxy.begin();
    for (int j=ystart; j<yend; ++j)
      for (int i=xstart; i<xend; ++i)
	{
	  *pres = *pdat - Star->flux * *ppsf - *pgal - Star->sky; // JG 
	  ++ppsf; ++ppdx; ++ppdy; ++pres; ++pdat; ++pgal; // JG
	}
    
  }else{
    for (int j=ystart; j<yend; ++j)
      for (int i=xstart; i<xend; ++i)
	{
	  *pres = *pdat - Star->flux * *ppsf - Star->sky;
	  ++ppsf; ++ppdx; ++ppdy; ++pres; ++pdat;
	}
  }
}
  

void SimFitRefVignet::SetStar(const PhotStar *RefStar) {
  Star = RefStar;
  //cout << " > SimfitRefVignet::SetStar x,y=" << Star->x << "," << Star->y << endl;
}

//=========================================================================================

void SimFitVignet::ResetFlags() {
  // state of components
  kernel_updated = false;
  psf_updated = false;
  resid_updated = false;
  gaussian_updated = false;
  // things that are to be fitted
  FitFlux = false;
  FitPos = false;
  UseGal = false;
  forceresize = true;
  // things that can be fitted
  CanFitFlux=true;
  CanFitSky=true;
  CanFitPos=true;
  CanFitGal=true;
}

void SimFitVignet::SetStar(const PhotStar *RefStar) {
  Star = RefStar;
  ResetFlags();
  //cout << " in SimfitVignet::SetStar x,y=" << Star->x << "," << Star->y << endl;
} 

SimFitVignet::SimFitVignet() {
  ResetFlags();
  kernelFit = 0;
}

SimFitVignet::SimFitVignet(const ReducedImage *Rim,  SimFitRefVignet* Ref)
  : Vignet(Rim)
{
  VignetRef = Ref;
  ResetFlags();

  double gain = Rim->Gain();
  inverse_gain = gain>0 ? 1./gain : 1.;
  ronoise = Rim->ReadoutNoise();
  skysub = Rim->OriginalSkyLevel();
  kernelFit = 0;
}

SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim,   SimFitRefVignet* Ref)
  : Vignet(Star, Rim, 0), FitFlux(false)
{  
#ifdef FNAME
  cout << " > SimFitVignet::SimFitVignet(const PhotStar *Star, const ReducedImage *Rim,  const SimFitRefVignet& Ref)" << endl;
#endif
  VignetRef = Ref;
  ResetFlags();
  inverse_gain = 1./Rim->Gain();
  kernelFit = 0;
}


void SimFitVignet::PrepareAutoResize() {
#ifdef FNAME
  cout << " > SimFitVignet::PrepareAutoResize()" << endl;
#endif
  if(!Star) {
    cerr << "SimFitVignet::PrepareAutoResize ERROR you need a star to update this vignet, use SetStar for this" << endl;
  }
  int hx_ref = VignetRef->Hx();
  int hy_ref = VignetRef->Hy();
  if(!kernel_updated)
    BuildKernel();
  int hx_kernel = Kern.HSizeX();
  int hy_kernel = Kern.HSizeY();
  //  Resize(hx_ref-hx_kernel,hy_ref-hy_kernel); // resize vignet, this will actually read the data
  int hx = hx_ref-hx_kernel;
  int hy = hy_ref-hy_kernel;
  int xc = int(Star->x);
  int yc = int(Star->y);
  
  // xstart,ystart,xend,yend
  Window window(xc-hx,yc-hy,xc+hx+1,yc+hy+1);
  
  reserve_vignet_in_server(FitsName(),window);
  if (HasWeight()) {
    reserve_vignet_in_server(FitsWeightName(),window); 
    if (HasSatur()) reserve_vignet_in_server(FitsSaturName(),window); 
  }
}

void SimFitVignet::AutoResize() {
#ifdef FNAME
  cout << " > SimFitVignet::AutoResize()" << endl;
#endif
  if(!Star) {
    cerr << "SimFitVignet::AutoResize ERROR you need a star to update this vignet, use SetStar for this" << endl;
  }
  int hx_ref = VignetRef->Hx();
  int hy_ref = VignetRef->Hy();
  if(!kernel_updated)
    BuildKernel();
  int hx_kernel = Kern.HSizeX();
  int hy_kernel = Kern.HSizeY();

#ifdef DEBUG
  SimFitRefVignet *toto = VignetRef;
  cout << "   in SimFitVignet::AutoResize VignetRef = " << toto << endl;
  cout << "   in SimFitVignet::AutoResize hx_ref    = " << hx_ref << endl;
  cout << "   in SimFitVignet::AutoResize hy_ref    = " << hy_ref << endl;
  cout << "   in SimFitVignet::AutoResize hx_kernel = " << hx_kernel << endl;
  cout << "   in SimFitVignet::AutoResize hy_kernel = " << hy_kernel << endl;  
#endif
  Resize(hx_ref-hx_kernel,hy_ref-hy_kernel); // resize vignet, this will actually read the data
}

void SimFitVignet::Update()
{
#ifdef FNAME
  cout << " > SimFitVignet::Update()" << endl;
#endif
  if(!Star) {
    cerr << "SimFitVignet::Update ERROR you need a star to update this vignet, use SetStar for this" << endl;
  }
  if(!kernel_updated)
    BuildKernel();
  
  if(!psf_updated)
    BuildPsf();
  
#ifdef DEBUG
  // check sizes
  //int hx_kernel = Kern.HSizeX();
  //int hy_kernel = Kern.HSizeY();
  int hx_data = Data.HSizeX();
  int hy_data = Data.HSizeY();
  int hx_weight = Weight.HSizeX();
  int hy_weight = Weight.HSizeY();
  int hx_resid = Resid.HSizeX();
  int hy_resid = Resid.HSizeY();
  int hx_psf = Psf.HSizeX();
  int hy_psf = Psf.HSizeY();
#ifndef ONEPSFPERIMAGE
  if (hx_data !=  hx_weight || hx_data != hx_resid || hx_data != hx_psf
      || hy_data !=  hy_weight || hy_data != hy_resid || hy_data != hy_psf ) {
#else
  if (hx_data !=  hx_weight || hx_data != hx_resid
      || hy_data !=  hy_weight || hy_data != hy_resid ) {
#endif
      cout << "   SimFitVignet::Update ERROR wrong size " << endl;
      cout << Kern << endl;
      cout << Data << endl;
      cout << Weight << endl;
      cout << Resid << endl;
      cout << Psf << endl;
      abort();
  }
#endif 

  if(!resid_updated) { 
#ifdef ONEPSFPERIMAGE
    BuildPsf(); // this loads a psf
    if(UseGal)
      UpdateResid_gal();
    else
      UpdateResid();
#else
    if(UseGal) 
      UpdateResid_psf_gal(); // compute psf as convolution of reference psf*kernel
    else
      UpdateResid_psf();
#endif
  }

  if(!gaussian_updated) {
    Psf.ComputeMoments(); // compute psf approximation with a gaussian
    gaussian_updated=true;
  }

}
  

  // check whether there are weights>0 on the position of the star,
  // if not, FitFlux=false 
void SimFitVignet::CheckWeight() {
  
  DPixel sumw=0;
  int center_hy = min(2,hy);
  int center_hx = min(2,hx);
  

  for (int j=-center_hy; j<=center_hy; ++j)
    for (int i=-center_hx; i<=center_hx; ++i)
      sumw += Weight(i,j);
  
  if(sumw<1.e-20) {
    cout << "WARNING SimFitVignet::CheckWeight : null weight at center, set FitFlux,FitSky to false" << endl;
    CanFitFlux=false;
    CanFitSky=false;
    CanFitPos=false;
    CanFitGal=false;
  }
}

// make sure that psf, vignet and data have proper sizes.
// allocate if needed, do only with a Star
 void SimFitVignet::Resize( int Hx_new,  int Hy_new)
{
#ifdef FNAME
  cout << " > SimFitVignet::Resize( int Hx_new,  int Hy_new)" << endl;
#endif
  if (!Star) return;
  
  // resize data, weight, resid if necessary
  if(Hx()==Hx_new && Hy()==Hy_new && !forceresize) {
#ifdef DEBUG
    cout << "    actually do not resize" << endl;
#endif
  }else{
    Vignet::Resize(Hx_new,Hy_new);
    OptWeight.Allocate(Nx(),Ny());
#ifndef ONEPSFPERIMAGE
    Psf.Allocate(Nx(),Ny());  // if use kernel, Psf is computed with resid
#else
    psf_updated    = false;
#endif
    kernel_updated = false;
    resid_updated  = false;
    forceresize = false;
    gaussian_updated = false;
  }
  Update();
}

void  SimFitVignet::BuildPsf() {
#ifdef FNAME
  cout << " > SimFitVignet::BuildPsf()" << endl;
#endif

#ifndef ONEPSFPERIMAGE
  
  if(!VignetRef->imagepsf) {
    cout << "SimFitVignet::BuildPsf ERROR VignetRef has no psf" <<endl;
    abort();
  }
  // we just want to allocate size for this Psf
  Psf.Resize(Hx(),Hy());

  //Psf.Tabulate(*Star,*(VignetRef->psf),*this);
  // now scale this according to kernel integral
  //if(!kernel_updated)
  //BuildKernel();
  //double photomratio = Kern.sum();
  //Psf*=photomratio;
  
#else
  
  if(!psf) {
    psf = new DaoPsf(*rim);
  }
  Psf.Tabulate(*Star,*psf,*this);
  // now scale this according to kernel integral
  if(!kernel_updated)
    BuildKernel();
  double photomratio = Kern.sum();
  Psf*=photomratio;
  
  gaussian_updated=true;
#endif
  
  psf_updated = true;
}

void SimFitVignet::BuildKernel()
{

#ifdef FNAME
  cout << " > SimFitVignet::BuildKernel(const ReducedImage* Ref)" << endl;
#endif
  if (!Star) return;
  
  if(! kernelFit ) {

#ifdef DEBUG_KERNEL
    cout << "KERNEL: no psfmatch in memory" << endl;
#endif
    const string kernelpath = rim->Dir()+"kernel_from_"+VignetRef->Image()->Name()+".dat";

    if(FileExists(kernelpath)) {
      
#ifdef DEBUG_KERNEL
      cout << "KERNEL: load from file" << endl;
#endif

      cout << "   Reading kernel " << kernelpath << " ..." << endl;
      kernelFit = new KernelFit();
      kernelFit->read(kernelpath);
      cout << "   done" << endl;
      
    } else
      {
	
#ifdef DEBUG_KERNEL
    cout << "KERNEL: compute kernel" << endl;
#endif
	
	cout << " SimFitVignet::BuildKernelPsf() : cannot find kernel " 
	     << kernelpath << ", so we do it" << endl;
	KernelFitter fitter(VignetRef->Image(), rim,true);
	fitter.DoTheFit();	  
	kernelFit = new KernelFit(fitter);
      }
  }
  
  kernelFit->KernAllocateAndCompute(Kern, Star->x, Star->y);
  Star->photomratio = Kern.sum();
  kernel_updated = true;

  double photom_ratio_threshold = 0.1;
  if(Kern.sum()<photom_ratio_threshold) {
    cout << "WARNING SimFitVignet::Buildkernel : photom_ratio too low :" << Kern.sum() 
	 << "; Set CanFitFlux,CanFitSky,CanFitPos,CanFitGal to 0." << endl;
    CanFitFlux=false;
    CanFitSky=false;
    CanFitPos=false;
    CanFitGal=false;
  }
}

void SimFitVignet::UpdateResid_psf_gal()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_psf_gal() : convolving Psf, Galaxy and updating residuals " << endl;
  cout << "Ref.Galaxy(0,0) = " << VignetRef->Galaxy(0,0) << endl;
#endif
  
  SimFitRefVignet& Ref = *VignetRef;

  // convolve all of them at same time  
  double sump, sumx, sumy, sumg;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy, *prgal;
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy, ++pw, ++pow)
	{
	  sumg = sump = sumx = sumy = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prpsf = &Ref.Psf   (i+hkx, j-jk);
	      prpdx = &Ref.Psf.Dx(i+hkx, j-jk);
	      prpdy = &Ref.Psf.Dy(i+hkx, j-jk);
	      prgal = &Ref.Galaxy(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik, ++pkern)
		{
		  sump += (*pkern) * (*prpsf); 
		  sumx += (*pkern) * (*prpdx);
		  sumy += (*pkern) * (*prpdy);
		  sumg += (*pkern) * (*prgal);
		  --prpsf; --prpdx; --prpdy; --prgal;
		}
	    }
	  
	  if( (!(sump>0)) && (!(sump<=0))) {
	    cout << "ERROR nan with sump in UpdateResid_psf_gal" << sump << endl;
	    DumpDebug();
	    abort();
	  }
	  val = Star->flux*sump+sumg+Star->sky;
	  *pres = *pdat - val;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
#ifdef VALCUTOFF
	  if(*pw == 0)
	    *pow = 0;
	  else {
	    if(val>VALCUTOFF)
	      *pow = 1./(1./(*pw)+val*inverse_gain);
	    else
	      *pow = *pw;
	  }
#else
	  *pow = *pw;
#endif
	  ++pdat;
	}
    }
  resid_updated = true;
}
  
// update psf and residuals with a convolved psf
void SimFitVignet::UpdateResid_psf()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_psf() : convolving Psf and updating residuals " << endl;
#endif
  // convolve all of them at same time  
  double sump, sumx, sumy;
  DPixel *pdat, *pres, *ppsf, *ppdx, *ppdy;
  DPixel *pkern, *prpsf, *prpdx, *prpdy;
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;
  
  const TabulatedPsf& RefPsf = VignetRef->Psf;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data   (-hx,j);
      pres = &Resid  (-hx,j);
      ppsf = &Psf    (-hx,j);
      ppdx = &Psf.Dx (-hx,j);
      ppdy = &Psf.Dy (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++ppsf, ++ppdx, ++ppdy, ++pw, ++pow)
	{
	  sump = sumx = sumy = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prpsf = &RefPsf   (i+hkx, j-jk);
	      prpdx = &RefPsf.Dx(i+hkx, j-jk);
	      prpdy = &RefPsf.Dy(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik, ++pkern) 
		{
		  sump += (*pkern) * (*prpsf); 
		  sumx += (*pkern) * (*prpdx); 
		  sumy += (*pkern) * (*prpdy);
		  --prpsf; --prpdx ; --prpdy;
		}
	    }
	  val =  Star->flux * sump + Star->sky;
	  *pres = *pdat - val;
	  *ppsf = sump;
	  *ppdx = sumx;
	  *ppdy = sumy;
#ifdef VALCUTOFF
	   if(*pw == 0)
	    *pow = 0;
	   else {
	     if(val>VALCUTOFF)
	       *pow  = 1./(1./(*pw)+val*inverse_gain);
	     else
	       *pow = *pw;
	   }
#else
	   *pow = *pw;
#endif
	   ++pdat;
	}
    }
   resid_updated = true;
}


void SimFitVignet::UpdateResid_gal()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid_gal() : convolving galaxy and updating residuals " << endl;
#endif
#ifdef DEBUG
  cout << " in SimFitVignet::UpdateResid Kern  = " << Kern  << endl;
  cout << " in SimFitVignet::UpdateResid Data  = " << Data  << endl;
  cout << " in SimFitVignet::UpdateResid Resid = " << Resid << endl;
  cout << " in SimFitVignet::UpdateResid Psf   = " << Psf   << endl;
  cout << " in SimFitVignet::UpdateResid Star->flux = " << Star->flux << endl;
#endif

  double sumg;
  DPixel *pdat, *pres, *ppsf, *pkern, *prgal;
  DPixel *pw,*pow;
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  double val;
  const Kernel& RefGal = VignetRef->Galaxy;

  for (int j=-hy; j<=hy; ++j)
    {
      pdat = &Data  (-hx,j);
      pres = &Resid (-hx,j);
      ppsf = &Psf   (-hx,j);
      pw   = &Weight (-hx,j);
      pow  = &OptWeight (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pres, ++pw, ++pow)
	{
	  sumg = 0.;
	  pkern = Kern.begin();
	  for (int jk =-hky; jk <= hky; ++jk)
	    {
	      prgal = &RefGal(i+hkx, j-jk);
	      for (int ik = -hkx; ik <= hkx; ++ik)
		{
		  sumg += (*pkern) * (*prgal);
		  ++pkern; --prgal;
		}
	    }
	  val = Star->flux* *ppsf + sumg + Star->sky;
	  *pres = *pdat - val;
#ifdef VALCUTOFF
	  if(*pw == 0)
	    *pow = 0;
	  else {
	    if(val>VALCUTOFF)
	      *pow = 1./(1./(*pw)+val*inverse_gain);
	    else
	      *pow = *pw;
	  }
#else
	  *pow = *pw;
#endif
	  ++pdat;  ++ppsf; 
	}
    }
  resid_updated = true;
}


// update residuals with existing psf, no galaxy
void SimFitVignet::UpdateResid()
{
#ifdef FNAME
  cout << " > SimFitVignet::UpdateResid() : updating residuals, no galaxy" << endl;
#endif
  
  DPixel *pdat = Data.begin(), *pres = Resid.begin(), *ppsf = Psf.begin();
  DPixel *pw = Weight.begin();
  DPixel *pow = OptWeight.begin();
  double val;
  for (int i=Nx()*Ny(); i; --i)
    {
      val = Star->flux * *ppsf + Star->sky;
      *pres = *pdat - val;
#ifdef VALCUTOFF
       if(*pw == 0)
	 *pow = 0;
       else {
	 if(val>VALCUTOFF)
	   *pow  = 1./(1./(*pw)+val*inverse_gain);
	 else
	   *pow=*pw;
       }
#else
       *pow=*pw;
#endif
       ++pres; ++pdat; ++ppsf; ++pw; ++pow;
    }
   resid_updated = true;
}

void SimFitVignet::RedoWeight()
{  
  //  if (Image()->HasWeight()) 
  //    Weight.readFromImage(Image()->FitsWeightName(), *this, 0);

  DPixel *pw   = Weight.begin();
  DPixel *pow  = OptWeight.begin();
  DPixel *pdat = Data.begin();
  DPixel *pres = Resid.begin();
  
  // avoid dividing by zero while keeping the zeros
  for (int i=Nx()*Ny(); i; --i) {
    double count = fabs((*pres - *pdat) * inverse_gain);
    double invweight = (*pw >0) ? 1./(*pw) : 0.;
    *pow =  (invweight>0) ? 1./(invweight+count) : 0.;
    ++pres;
    //double var = (*pdat + skysub)*inverse_gain + ronoise*ronoise;
    //*pow = 1./var;
    ++pow; ++pdat; ++pw;
  }
}



void SimFitVignet::DumpDebug() const {  
  cout << "############# SimFitVignet::DumpDebug #############" << endl;
  string name =  Image()->Name();
  cout << " image = " << name << endl; 
  cout << " star = " << *Star << endl; 
  cout << " galaxy(0,0) = " << VignetRef->Galaxy(0,0) << endl;
  cout << " Data(0,0) = " << Data(0,0) << endl;
  cout << " Psf(0,0) = " << Psf(0,0) << endl;
  cout << " writing *_DEBUG.fits" << endl;
  Weight.writeFits(name+"_weight_DEBUG.fits");
  Resid.writeFits(name+"_resid_DEBUG.fits");
  Data.writeFits(name+"_data_DEBUG.fits");
  Psf.writeFits(name+"_psf_DEBUG.fits");
  VignetRef->Psf.writeFits(name+"_refpsf_DEBUG.fits");
  Kern.writeFits(name+"_kern_DEBUG.fits");
}

static double sqr(const double& x) { return x*x; }

double SimFitVignet::Chi2() const {
  
  double chi2 = 0;
  DPixel *pow=OptWeight.begin(), *pres=Resid.begin();//, *pdat = Data.begin();
  for (int i=Nx()*Ny(); i; --i, ++pres, ++pow) {
    chi2 += *pow * sqr(*pres);

    // do this when redoing weight with Poisson noise from star & galaxy
    // Gauss MLE, NIMPA A 457, p.394, eq (28) 
    //double invw = (*pow)>0 ? 1./(*pow) : 1;
    //double ci = skysub + *pdat++;
    //double cip = sqrt(0.25 + sqr(ci)) - 0.5;
    //chi2 += log(invw/cip) - sqr(ci - cip)/cip;
  }

  if (chi2 >= 0) return chi2;

  // now debug 
  cout << "############# SimFitVignet::Chi2 ERROR chi2=" << chi2 << " #############" << endl;
  DumpDebug();
  cout << " now I quit ..." << endl;
  abort();
}
 
double SimFitVignet::CentralChi2(int &npix) const {
  
  double chi2 = 0;
  
  int chx = int(ceil(sqrt(Psf.Mxx())));
  int chy = int(ceil(sqrt(Psf.Myy())));
  if(chx>hx)
    chx=hx;
  if(chy>hy)
    chy=hy;
  
  
  for(int j=-chy;j<=chy;j++)
    for(int i=-chx;i<=chx;i++)
      chi2 += OptWeight(i,j)*Resid(i,j)*Resid(i,j);
  
  npix = (2*chx+1)*(2*chy+1);
  return chi2;
}
 




//=============================================================
