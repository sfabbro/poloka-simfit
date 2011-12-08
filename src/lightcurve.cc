#include <algorithm>  // min_element, copy, for_each
#include <functional> // bind2nd
#include <iterator>   // ostream_iterator
#include <iomanip>    // setw, fixed ...
#include <fstream>
#include <fitsimage.h>
#include <reducedimage.h>
#include "polokaexception.h"

#include "lightcurve.h"
#include "lcio.h"

// for io
#include "lightcurvepoint.h"
//#include "lightcurvepoint_dict.h"
//#include "objio.h"
//#include "typemgr.h"


// instantiate
template class Fiducial<PhotStar>;

void LightCurve::push_back(const ReducedImage* Rim, const PhotStar *Star)
{
  // look if image is already there
  LightCurve::iterator it = find_if(begin(), end(), bind2nd(mem_fun(&Fiducial<PhotStar>::HasImage), Rim));

  if (it != end())
    {
      cerr << " LightCurve::push_back() : Warning: image "
	   << Rim->Name() << " is already in the LightCurve " << endl;
      return;
    }

  Fiducial<PhotStar> *fidStar = new Fiducial<PhotStar>(Star, Rim);

  push_back(fidStar);
}

void LightCurve::push_back(const ReducedImage* Rim)
{
  // look if image is already there
  LightCurve::iterator it = find_if(begin(), end(), bind2nd(mem_fun(&Fiducial<PhotStar>::HasImage), Rim));

  if (it != end())
    {
      cerr << " LightCurve::push_back() : Warning: image "
	   << Rim->Name() << " is already in the LightCurve " << endl;
      return;
    }

  Fiducial<PhotStar> *fidStar = new Fiducial<PhotStar>(Rim);

  push_back(fidStar);
}

void LightCurve::write_lc2fit(ostream& Stream) const
{
 
  double elixir_zp = computeElixirZeroPoint();

  ios::fmtflags oldflags = Stream.flags();

  if (front()->Image()) Stream << "#Date : (MJD)\n";
  Stream << "#Flux : in units of ADU in reference image\n"  
         << "#Fluxerr : \n"
	 << "#ZP : elixir zp\n"
	 << "#seeing: SEseeing\n"
	 << "#exptime : exposure time\n"
	 << "#phratio : photom ratio\n"
	 << "#gseeing : GFseeing\n"
	 << "#sesky : SEsky\n"
	 << "#sigsky : SIGsky\n"
	 << "#sigscale : sigma scale factor\n";
  
  if (front()->Image()) Stream << "#Image : \n";
  Stream << "@INSTRUMENT MEGACAM\n";
  Stream << "@BAND " << Ref->band << "\n";
  Stream << "@MAGSYS AB\n";
  Stream.setf(ios::fixed);

  LightCurvePoint lcp;
  for (LightCurve::const_iterator it = begin(); it != end(); ++it)
    {      
      const Fiducial<PhotStar> *fs = *it;

      // can generate problems if fitted flux is by chance in this range
      if(fabs(fs->flux)<1.e-30) // do not print unfitted fluxes
	continue;
      lcp.modifiedjulianday = fs->ModifiedJulianDate();
      lcp.flux = fs->flux;
      lcp.eflux = sqrt(fs->varflux);
      //lcp.computemag(elixir_zp);
      lcp.zeropoint = elixir_zp;
      Stream << lcp;
      Stream << " " << fs->Seeing();
      Stream << " " << fs->ExposureTime();
      Stream << " " << fs->photomratio;
      try {
	Stream << " " << fs->GFSeeing();
      } catch (PolokaException p) {
	p.PrintMessage(cout);
	Stream << " " << 0.;
      }
      Stream << " " << fs->SESky();
      Stream << " " << fs->SIGSky();
      Stream << " " << fs->sigscale_varflux; 
      if (fs->Image()) 
	{
	  string aligned_name = fs->Name() ; 
	  size_t align_pos = aligned_name.find("enlarged");
	  string dbim_name;
	  if(align_pos!=string::npos) dbim_name = aligned_name.erase(0,align_pos+8); 
	  size_t pos = dbim_name.find("p");
	  if(pos!=string::npos) dbim_name.replace(pos,1,"");
	  Stream << "  " << dbim_name;
	}
      else 
	if(front()->Image())
	  Stream << " none ";
      Stream << endl;
    }
  Stream.flags(oldflags);
}


ostream& operator << (ostream& Stream, const CountedRef<Fiducial<PhotStar> > &Star)
{
  Stream << *Star;
  return Stream;
}

ostream& operator << (ostream& Stream, const LightCurve& Lc)
{
  copy(Lc.begin(), Lc.end(), ostream_iterator<CountedRef<Fiducial<PhotStar> > >(Stream, "\n"));
  return Stream;
}


LightCurveList::LightCurveList(istream& LcFileStream) 
{
#ifdef FNAME
  cout << " > LightCurveList::LightCurveList(istream& LcFileStream)" << endl;
#endif
  lc_read(LcFileStream, Objects, Images);
  
  RefImage = Objects.front()->Image();
  cout << " > LightCurveList::LightCurveList() : read " 
       << Objects.size() << " objects, " 
       << Images.size() << " images. Reference is " 
       << RefImage->Name() << endl;

  // fill up the light curve
  int nobj=0;
  for (RefStarCIterator it = Objects.begin(); it != Objects.end(); ++it)
    {
      LightCurve lc(*it);
      // foreach object, link the list of images
      cout << flush << " > LightCurveList::LightCurveList() : filling object " << nobj++ << "\r";
      for (ReducedImageCIterator im=Images.begin(); im != Images.end(); ++im) {
	PhotStar *fidPhot = new PhotStar(BaseStar((*it)->x, (*it)->y, 0.));
	lc.push_back(*im, fidPhot); // add one image and one PhotStar
      }
      push_back(lc);
    }  
  
  cout << endl;
}

ostream& operator << (ostream& Stream, const LightCurveList& Fiducials)
{
  copy(Fiducials.begin(), Fiducials.end(), ostream_iterator<LightCurve>(Stream));
  return Stream;
}




double LightCurve::computeElixirZeroPoint() const {
  // > COMMENT   Formula for Photometry, based on keywords given in this header:
  // > COMMENT   m = -2.5*log(DN) + 2.5*log(EXPTIME)
  // > COMMENT   M = m + PHOT_C + PHOT_K*(AIRMASS - 1) + PHOT_X*(PHOT_C1 - PHOT_C2)
  
  
  string photometric_image_fitsname = Ref->Image()->FitsName(); // default
  double photomratio = 1; // flux(photometric_image)/flux(reference_image=Ref->Image())
  
  double expo,PHOT_C,PHOT_K,AIRMASS;
  expo=1;
  PHOT_C=0;
  PHOT_K=0;
  AIRMASS=1;
  
  FitsHeader refhead(photometric_image_fitsname);
  if(refhead.HasKey("PHOT_C")) {
    expo =  refhead.KeyVal("TOADEXPO");
    PHOT_C =  refhead.KeyVal("PHOT_C");
    PHOT_K =  refhead.KeyVal("PHOT_K");
    AIRMASS =  refhead.KeyVal("AIRMASS");
  }else{
    cout << "LightCurve::computeElixirZeroPoint WARNING " << photometric_image_fitsname << " has no PHOT_C, try another one ..." << endl;
    for (LightCurve::const_iterator it = begin(); it != end(); ++it) {
      const Fiducial<PhotStar> *fs = *it;
      FitsHeader head(fs->Image()->FitsName());
      if(head.HasKey("PHOT_C")) {
	expo =  head.KeyVal("TOADEXPO");
	PHOT_C =  head.KeyVal("PHOT_C");
	PHOT_K =  head.KeyVal("PHOT_K");
	AIRMASS =  head.KeyVal("AIRMASS");
	photometric_image_fitsname = fs->Image()->FitsName();
	photomratio = fs->photomratio;
	cout << "LightCurve::computeElixirZeroPoint using " << photometric_image_fitsname << endl;
	break;
      }
    }
  }
  if(PHOT_C==0) {
    cout << "LightCurve::computeElixirZeroPoint NO ELIXIR ZP AT ALL !!!! " << endl;
    return 0;
  }
  return 2.5*log10(expo) + PHOT_C + PHOT_K*(AIRMASS-1.) + 2.5*log10(photomratio); // if photomratio>1, flux in ref < flux photometric => mag in ref > mag  photometric
}

void LightCurveList::write(const string& filename) const
{
  ofstream out(filename.c_str());
  out << "# jdmin : \n"
      << "# jdmax : \n"
      << "# x :\n"
      << "# y :\n"
      << "# flux :\n"
      << "# sky :\n"
      << "# varx :\n"
      << "# vary :\n"
      << "# covxy :\n"
      << "# varflux :\n"
      << "# varsky :\n"
      << "# totflux : \n"
      << "# totsky :\n"
      << "# galflux :\n"
      << "# vargalflux :\n"
      << "# vartotflux :\n"
      << "# vartotsky :\n"
      << "# chi2 :\n"
      << "# ndf : \n"
      << "# resmean :\n"
      << "# end \n";

  for (LightCurveList::const_iterator it = begin(); it != end(); ++it) {
    if (fabs(it->Ref->jdmin)>1e11 && fabs(it->Ref->jdmax)>1e11) {
      out.setf(ios::scientific);
      out.unsetf(ios::fixed);
      out << setprecision(1) << it->Ref->jdmin << ' '
	  << setprecision(1) << it->Ref->jdmax << ' ';
    } else {
      out.unsetf(ios::scientific);
      out.setf(ios::fixed); 
      out << setprecision(4) << setw(11) << it->Ref->jdmin << ' '
	  << setprecision(4) << setw(11) << it->Ref->jdmax << ' ';
    }
    out.unsetf(ios::scientific);
    out.setf(ios::fixed); 
    out << setprecision(4) << setw(11) << it->Ref->x << ' '
	<< setprecision(4) << setw(11) << it->Ref->y << ' '
	<< setprecision(4) << setw(11) << it->Ref->flux << ' '
	<< setprecision(4) << setw(11) << it->Ref->sky << ' '
	<< setprecision(7) << setw(9) << it->Ref->varx << ' '
	<< setprecision(7) << setw(9) << it->Ref->vary << ' '
	<< setprecision(7) << setw(9) << it->Ref->covxy << ' '
	<< setprecision(6) << setw(13) << it->Ref->varflux << ' '
	<< setprecision(6) << setw(13) << it->Ref->varsky << ' '
	<< setprecision(2) << setw(12) << it->totflux << ' '
	<< setprecision(2) << setw(12) << it->totsky << ' '
	<< setprecision(2) << setw(12) << it->galflux << ' '
	<< setprecision(2) << setw(12) << it->vargalflux << ' '
	<< setprecision(2) << setw(12) << it->vartotflux << ' '
	<< setprecision(2) << setw(12) << it->vartotsky << ' '
	<< setprecision(4) << setw(11) << it->chi2 << ' '
	<< setw(11) << it->ndf << ' '
	<< setprecision(4) << setw(11) << it->resmean
	<< endl;
  }
    
}
