#include <iomanip>

#include <poloka/sestar.h>

#include <poloka/photstar.h>

PhotStar::PhotStar()
  : BaseStar(0,0,0), sky(0.), varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false;n_saturated_pixels = 0 ; image_seeing=0;}

PhotStar::PhotStar(const BaseStar &BStar)
  : BaseStar(BStar), sky(0.), varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false; n_saturated_pixels = 0 ;image_seeing=0;}

PhotStar::PhotStar(const SEStar &SStar)
  : BaseStar(SStar), sky(SStar.Fond()), 
    varsky(0.), photomratio(1.), sigscale_varflux(-1)
{has_saturated_pixels = false; n_saturated_pixels = 0 ;image_seeing=0;}

void PhotStar::writen(ostream &Stream) const
{
  int old = Stream.precision();
  BaseStar::writen(Stream); 
  Stream << ' ' 
	 << sky     << ' '
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
  Stream << " sky " << sky << " varsky " << varsky;
}

string PhotStar::WriteHeader_(ostream & Stream, const char* Num) const
{
  if (!Num) Num = "";
  const string baseStarFormat = BaseStar::WriteHeader_(Stream,Num);
  Stream << "# sky"     << Num << " : mean sky value per pixel" << endl 
	 << "# seeing"     << Num << " : image seeing " << endl
	 << "# phratio" << Num << " : initial image photometric ratio " << endl
	 << "# sigscale"    << Num << " : flux variance rescaling factor" << endl ;

  
  // 1: add photomratio, sigscale, seeing
  // 2: remove all varflux and var pos, now in basestar
  
  return baseStarFormat + "PhotStar 2 ";
}

// Force instanciation
//#include <starlist.cc>
//template class StarList<PhotStar>; 

