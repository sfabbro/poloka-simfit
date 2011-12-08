#include "refstar.h"


ostream& operator << (ostream &Stream, const RefStar &Star)
{
  Star.writen(Stream);
  return Stream;
}

/*
RefStar* RefStar::read(istream & Stream, const char *Format)
{
  RefStar *rstar = new RefStar;
  Stream >> rstar->name 
	 >> rstar->type
	 >> rstar->jdmin
	 >> rstar->jdmax
	 >> rstar->ra
	 >> rstar->dec 
	 >> rstar->varra 
	 >> rstar->vardec
	 >> rstar->covradec
	 >> rstar->band;
  rstar->PhotStar::read(Stream, Format);
  return rstar;
}
*/


void RefStar::writen(ostream &Stream) const
{
  Stream << name << ' '
	 << type << ' '
	 << jdmin << ' '
	 << jdmax << ' '
	 << ra << ' '
	 << dec << ' '
	 << varra << ' '
	 << vardec << ' '
	 << covradec << ' '
  	 << band << ' ';
  PhotStar::writen(Stream);
}

string RefStar::WriteHeader_(ostream & Stream, const char* Num) const
{
  if (!Num) Num = "";
  Stream << "# name : star name \n"
	 << "# type : star type (0:sn+gal, 1:star, 2:gal) \n"
	 << "# jdmin : julian date of birth \n"
	 << "# jdmax : julian date of death \n"
	 << "# ra :  right ascension (J2000) \n"
	 << "# dec : declination (J2000) \n"
	 << "# varra : variance of ra \n"
	 << "# vardec : variance of dec \n"
	 << "# covradec: covariance term ra,dec \n"
	 << "# band: band (char) \n";
  
  return PhotStar::WriteHeader_(Stream, Num) + "RefStar 0 ";
}

// instantiate
//#include <starlist.cc>
//template class StarList<RefStar>;

