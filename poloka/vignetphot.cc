#include <poloka/reducedimage.h>
#include <poloka/vutils.h>

#include <poloka/vignetphot.h>

static double sqr(const double &x) { return x*x; }

void VignetPhot::SetRadius(const double& Rad)
{
  rad = Rad;
  radAperMin2 = sqr(rad-0.5);
  radAperMax2 = sqr(rad+0.5);
}


void VignetPhot::SetSkyRadius(const double& RadSkyMin, const double& RadSkyMax)
{
  radSkyMin2  = sqr(RadSkyMin);
  radSkyMax2  = sqr(RadSkyMax);
  if (skyArray) delete skyArray;
  skyArray = new double[int(4*radSkyMax2)];
}


void VignetPhot::Aperture()
{

  pim = Data.begin();
  pw  = Weight.begin();
  sumf = sumwf = 0.;

  for(int j=ystart; j<yend; ++j)
    {
      dy2 = sqr(j - Star->y);
      for(int i=xstart; i<xend; ++i, ++pw, ++pim)
	{
	  if (*pw == 0.) continue;
	  dist2 = sqr(i-Star->x) + dy2;
	  if (dist2 > radAperMax2) continue;
	  wf = *pw;
	  if (dist2 > radAperMin2) wf *= 0.5 - sqrt(dist2) + rad;
	  sumf  += wf * (*pim - Star->sky);
	  sumwf += wf;
	} 
    }

  if (sumwf <= 0.) return;
  Star->flux  = sumf / sumwf;
  Star->eflux = sqrt(1. / sumwf);

}

void VignetPhot::SkyAnnulus()
{
  pim = Data.begin();
  pw  = Weight.begin();

  sumf = sumwf = sumx2 = 0.;

  double *psky = &skyArray[0];
  int nsky = 0;

  for(int j=ystart; j<yend; ++j)
    {
      dy2 = sqr(j - Star->y);
      for(int i=xstart; i<xend; ++i, ++pw, ++pim)
	{
	  if (*pw == 0.) continue;
	  dist2 = sqr(i-Star->x) + dy2;
	  if (dist2 > radSkyMax2  || dist2 < radSkyMin2) continue;
	  *psky++ = *pim;
	  nsky++;
	  sumf  += *pim;
	  sumwf += *pw;
	  sumx2  = sqr(*pim);
	}
    }

  if (sumwf <= 0.) return;
  double skyMean   = sumf / sumwf;
  double skyMedian = DArrayMedian(skyArray, int(nsky));

  if (skyMedian < skyMean) Star->sky = 3.*skyMedian - 2.*skyMean;
  else Star->sky = skyMean;      
  Star->varsky = sumx2/nsky - sqr(skyMean);
}


void VignetPhot::Photometry()
{
  pim = Data.begin();
  pw  = Weight.begin();
  sumf = sumwf = sumx2 = sumy2 = sumxy = 0.;
  
  for(int j=ystart; j<yend; ++j)
    for(int i=xstart; i<xend; ++i, ++pw, ++pim)
      {
	p  = psf->Value(i, j, Star->x, Star->y, dpdx, dpdy);
	wf = p    * (*pw);
	wx = dpdx * (*pw);
	wy = dpdy * (*pw);
	sumf  += wf * (*pim - Star->sky);
	sumwf += wf * p;
	sumx2 += wx * dpdx;
	sumy2 += wy * dpdy;
	sumxy += wx * dpdy;
      }
  
  if (sumwf <= 0.) return;
  Star->flux  = sumf / sumwf;
  Star->eflux = sqrt(1. / sumwf);
  
  // compute variance position only if significative flux
  if (fabs(Star->flux) > 0.)
    {
      double det = sqr(Star->flux) * (sumx2*sumy2 - sqr(sumxy)); 
      Star->vx  = sumx2 / det;
      Star->vy  = sumy2 / det;
      Star->vxy = sumxy / det;
    }
}

void VignetPhot::Recentroid(const int MaxIter, const double& Eps)
{  

  double det;
  int iter=0;
  do
    {
      xp = Star->x;
      yp = Star->y;
      
      pim = Data.begin(); 
      pw  = Weight.begin();
      sumxy = sumx2 = sumy2 = sumxr = sumyr = 0.;
      
      for(int j=ystart; j<yend; ++j)
	for(int i=xstart; i<xend; ++i, ++pw, ++pim)
	  {
	    res = *pim - Star->flux * psf->Value(i, j, Star->x, Star->y, dpdx, dpdy) - Star->sky;
	    wx  = dpdx * (*pw);
	    wy  = dpdy * (*pw);
	    sumxr += wx * res;
	    sumyr += wy * res;
	    sumx2 += wx * dpdx;
	    sumy2 += wy * dpdy;
	    sumxy += wx * dpdy;
	  }
      det = Star->flux * (sumx2*sumy2 - sqr(sumxy));
      Star->x += (sumxr*sumx2 - sumyr*sumxy) / det;
      Star->y += (sumyr*sumy2 - sumxr*sumxy) / det;
    }
  while ( ((fabs(xp-Star->x) > Eps) || (fabs(yp-Star->y) > Eps)) && (iter++ < MaxIter) );

}


void VignetPhot::operator () (const PhotStar *Star)
{
  Load(Star);
  SkyAnnulus();
  Recentroid();
  Photometry();
  Aperture();
}
