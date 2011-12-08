// This may look like C code, but it is really -*- C++ -*-
#ifndef PHOTSTAR__H
#define PHOTSTAR__H

#include <starlist.h>
#include <basestar.h>

class SEStar;

class PhotStar : public BaseStar {

public:

  PhotStar();

  PhotStar(const BaseStar &BStar);

  PhotStar(const SEStar &SStar);

  double sky, varflux, varx, vary, covxy, varsky, photomratio, sigscale_varflux;
  bool has_saturated_pixels;
  int n_saturated_pixels;
  double image_seeing;
  double MJD; 
  
  // photomratio needed to keep the information from simfitvignet -> lightcurve (we have to move this somewhere else)

  friend ostream& operator << (ostream &Stream, const PhotStar &PStar)
   { PStar.dumpn(Stream); return Stream; }


  // all the shitty I/O routines for StarList

  static PhotStar* read(istream& Stream, const char *Format);

  virtual void  dumpn(ostream& Stream=cout) const;

  virtual void writen(ostream& Stream=cout) const;

  string WriteHeader_(ostream& Stream, const char* Num=0) const;

  ostream& write_header(ostream& Stream) const { WriteHeader_(Stream); return Stream; }
  
  static const char* TypeName() {return "PhotStar";}
};

typedef CountedRef<PhotStar>         PhotStarRef;
typedef StarList<PhotStar>           PhotStarList;
typedef PhotStarList::iterator       PhotStarIterator;
typedef PhotStarList::const_iterator PhotStarCIterator;


#endif // PHOTSTAR__H

