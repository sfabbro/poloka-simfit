// This may look like C code, but it is really -*- C++ -*-
#ifndef REFSTAR__H
#define REFSTAR__H

#include <starlist.h>

#include "fiducial.h"
#include "photstar.h"

class RefStar : public Fiducial<PhotStar> {

public:

  RefStar() : name(""), type(0), jdmin(0), jdmax(0), varra(0), vardec(0), covradec(0), band('?')  {}

  RefStar(const ReducedImage *Rim) : Fiducial<PhotStar>(Rim) {}
  
  string name;
  int type;
  double jdmin,jdmax;
  double ra,dec;
  double varra, vardec, covradec;
  char band;
  
  bool IsVariable(const double& Jd) const { return (Jd <= jdmax && Jd >= jdmin); }

  friend ostream& operator << (ostream &Stream, const RefStar &Star);

  // all the shitty I/O routines for StarList

  static RefStar* read(istream& Stream, const char *Format);

  virtual void writen(ostream& Stream=cout) const;

  string WriteHeader_(ostream& Stream, const char* Num=0) const;
  
  static const char* TypeName() {return "RefStar";}
};

typedef CountedRef<RefStar>         RefStarRef;
typedef StarList<RefStar>           RefStarList;
typedef RefStarList::iterator       RefStarIterator;
typedef RefStarList::const_iterator RefStarCIterator;


#endif // REFSTAR__H
