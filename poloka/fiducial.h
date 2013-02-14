// This may look like C code, but it is really -*- C++ -*-
#ifndef FIDUCIAL__H
#define FIDUCIAL__H

#include <iomanip>

#include <poloka/countedref.h>
#include <poloka/reducedimage.h>

//! a template to use when an pointer element belongs to an image
template<class S> 
class Fiducial : public S {

protected:

  CountedRef<ReducedImage> rim;
  double seeing,mjd;
  double exptime;
  double sesky, sigsky;
  int has_weight,has_satur;
  string fits_name,fits_weight_name,fits_satur_name;
public:

  Fiducial() : sesky(-1),sigsky(-1),seeing(-1),mjd(-1),has_weight(-1),has_satur(-1), exptime(-1) {}

  Fiducial(const S *Fid) : S(*Fid) {}

  Fiducial(const ReducedImage *Rim) : rim(Rim), sesky(-1),sigsky(-1),seeing(-1), mjd(-1), exptime(-1), has_weight(-1),has_satur(-1) {}

  Fiducial(const S *Fid, const ReducedImage *Rim) : S(*Fid), rim(Rim), sesky(-1),sigsky(-1),seeing(-1), mjd(-1), exptime(-1), has_weight(-1),has_satur(-1) {}

  const ReducedImage* Image() const { return rim; }

  bool HasImage(const ReducedImage *Rim) const
  {
    return (rim->Name() == Rim->Name());
  }

  string Name() const { return rim->Name();}

  void AssignImage(const ReducedImage *Rim) 
  { 
    if (rim) { cerr << " Fiducial::AssignImage() : Error : " 
		    << rim->Name() << " already assigned \n";  return; }
    rim = Rim;
    seeing = -1;
    mjd = -1;
    exptime= -1;
    //cout << "in AssignImage rim=" << rim->Name() << endl;
  }

  ostream& write_header(ostream &Stream=cout) const
  {
    if (rim) 
      {
	Stream << "# name : DbImage name \n";
	Stream << "# band : integer for the band (ascii code) \n";
      }
    S::write_header(Stream);
    Stream << "# end \n";

    return Stream;
  }

  friend ostream& operator << (ostream &Stream, const Fiducial<S>& Fs)
  {
    ios::fmtflags oldflags = Stream.flags();
    Stream.setf(ios::fixed);
    if (Fs.rim)
      Stream << setw(12) << setprecision(2) << Fs.ModifiedJulianDate()
	     << int(Fs.rim->Band()[0]);
    Stream << (S) (Fs);
    Stream.flags(oldflags);    
    return Stream;
  }

  friend bool IncreasingJulianDate(const Fiducial<S> *one, const Fiducial<S> *two)
   { return (one->ModifiedJulianDate() < two->ModifiedJulianDate()); }

  friend bool IncreasingSeeing(const Fiducial<S> *one, const Fiducial<S> *two)
   { return (one->Seeing() < two->Seeing()); }

  double Seeing() const {
    if(seeing<0) {
      const_cast<double&>(seeing) = rim->Seeing();  
    }
    return seeing;
  }
  
  double ModifiedJulianDate() const {
    if(mjd<0) {
      const_cast<double&>(mjd) = rim->ModifiedJulianDate();  
    }
    return mjd;
  }
  
  double ExposureTime() const {
    if(exptime<0) {
      const_cast<double&>(exptime) = rim->Exposure();  
    }
    return exptime;
  }

  double SESky() const {
    if(sesky<0) {
      const_cast<double&>(sesky) = rim->BackLevelNoSub();  
    }
    return sesky;
  }

  double SIGSky() const {
    if(sigsky<0) {
      const_cast<double&>(sigsky) = rim->SigmaBack();  
    }
    return sigsky;
  }

  bool HasWeight() const { 
    if(has_weight<0) 
      const_cast<int&>(has_weight) = (rim->HasWeight())?1:0;  
    return has_weight==1;
  }

  bool HasSatur() const { 
    if(has_satur<0) 
      const_cast<int&>(has_satur) = (rim->HasSatur())?1:0;  
    return has_satur==1;
  }
  string FitsName() const { 
    if(fits_name=="")
      const_cast<string&>(fits_name) = rim->FitsName();
    return fits_name;
  }
  string FitsWeightName() const { 
    if(fits_weight_name=="")
      const_cast<string&>(fits_weight_name) = rim->FitsWeightName();
    return fits_weight_name;
  }
  string FitsSaturName() const { 
    if(fits_satur_name=="")
      const_cast<string&>(fits_satur_name) = rim->FitsSaturName();
    return fits_satur_name;
  }
  

};


#endif // FIDUCIAL__H
