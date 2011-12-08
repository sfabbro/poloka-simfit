#include <iomanip>

#include <sestar.h>

#include "photstar.h"

PhotStar::PhotStar()
  : BaseStar(0,0,0), sky(0.), varflux(0.), varx(0.), vary(0.), 
    covxy(0.), varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false;n_saturated_pixels = 0 ; image_seeing=0;}

PhotStar::PhotStar(const BaseStar &BStar)
  : BaseStar(BStar), sky(0.), varflux(0.), varx(0.), vary(0.), 
    covxy(0.), varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false; n_saturated_pixels = 0 ;image_seeing=0;}

PhotStar::PhotStar(const SEStar &SStar)
  : BaseStar(SStar), sky(SStar.Fond()), varflux(SStar.EFlux()*SStar.EFlux()),
    varx(0.), vary(0.), covxy(0.), varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false; n_saturated_pixels = 0 ;image_seeing=0;}

void PhotStar::writen(ostream &Stream) const
{
  int old = Stream.precision();
  BaseStar::writen(Stream); 
  Stream << ' ' 
	 << sky     << ' '
	 << varflux << ' '
	 << varx    << ' '
	 << vary    << ' '
	 << covxy   << ' '
	 << varsky; 
  Stream << ' ' 
	 << image_seeing     << ' '
	 << photomratio << ' '
	 << sigscale_varflux << ' ';

  Stream << setprecision(old);
}


/*
PhotStar* PhotStar::read(istream & Stream, const char *Format)
{
  PhotStar *pstar = new PhotStar();
  pstar->BaseStar::read(Stream, Format);

  //  if (Stream >> pstar->sky 
      >> pstar->varflux
      >> pstar->varx 
      >> pstar->vary 
      >> pstar->covxy 
      >> pstar->varsky) return pstar;
 if (Stream >> pstar->sky 
      >> pstar->varflux
      >> pstar->varx 
      >> pstar->vary 
      >> pstar->covxy 
      >> pstar->varsky
      >> pstar->image_seeing
      >> pstar->photomratio
      >> pstar->sigscale_varflux) return pstar;

  return 0;    
}
*/


void PhotStar::dumpn(ostream &Stream) const
{
  BaseStar::dumpn(Stream);
  Stream << " sky " << sky << " varflux " << varflux << " varx " << varx << " vary " << vary 
	 << " covxy " << covxy << " varsky " << varsky;
}

string PhotStar::WriteHeader_(ostream & Stream, const char* Num) const
{
  if (!Num) Num = "";
  const string baseStarFormat = BaseStar::WriteHeader_(Stream,Num);
  Stream << "# sky"     << Num << " : mean sky value per pixel" << endl 
	 << "# varflux" << Num << " : variance of measured flux" << endl
	 << "# varx"    << Num << " : variance in x position (pixels)" << endl 
	 << "# vary"    << Num << " : variance in y position (pixels)" << endl 
	 << "# covxy"   << Num << " : covariance in x and y position (pixels)" << endl 
	 << "# varsky"  << Num << " : variance in sky measurement" << endl;
  
  Stream << "# seeing"     << Num << " : image seeing " << endl
	 << "# phratio" << Num << " : initial image photometric ratio " << endl
	 << "# sigscale"    << Num << " : flux variance rescaling factor" << endl ;

  
  // return baseStarFormat + "PhotStar 0 ";
  return baseStarFormat + "PhotStar 1 "; // on ajoute le write de photomratio et de sigscale et de seing
}

// Force instanciation
//#include <starlist.cc>
//template class StarList<PhotStar>; 

