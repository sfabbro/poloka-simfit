#include <string>

using namespace std;
#include <algorithm>


#include <astroutils.h>   // JulianDay
#include "lcio.h"

// routine to read an object in a lightcurve file. 
// so far kept local to this program
static RefStar* read_object(const string& line)
{
  vector<string> words;
  DecomposeString(words, line);
  if (words.size() < 2)
    {
      cerr << " read_object() Error: wrong format in line \n'"
	   << line << "'\n";
      return false;
    }

  RefStar* star = new RefStar;
  star->type=0; // defalut is star+galaxy
  star->x = atof(words[0].c_str());
  star->y = atof(words[1].c_str());
  vector<string>::const_iterator it=words.begin();
  ++it; ++it;
  for ( ; it != words.end(); ++it)
    {
      vector<string> option;
      DecomposeString(option, *it, "=");
      if      (option[0]=="DATE_MIN") star->jdmin = atof(option[1].c_str());
      else if (option[0]=="DATE_MAX") star->jdmax = atof(option[1].c_str());
      else if (option[0]=="NAME")  star->name  = option[1];
      else if (option[0]=="TYPE")  star->type  = atoi(option[1].c_str());
      else if (option[0]=="BAND")  star->band  = option[1][0];      
      else cerr << " read_object() Error: unknown argument :'" << option[0] 
		<< "' in line \n   '" << line << "'\n";
    }

  if ((star->jdmin > 0) && (star->jdmax > 0) && (star->jdmin >= star->jdmax))
    {
      cerr << " read_object() Error: \n";
      cerr << "    dates are not correct " << endl;
      return false;
    }
  //cout << "jdmin jdmax = " << star->jdmin << " " << star->jdmax << endl;
  return star;
}

// routine to read a DbImage in a lightcurve file. 
// so far kept local to this program
static ReducedImage* read_image(const string& line)
{
  vector<string> words;
  DecomposeString(words, line);
  RemovePattern(words.front()," ");
  cout << "new img " << words.front() << endl;
  ReducedImage *red = new ReducedImage(words.front());
  //cout << "red = " << red << endl;

  if (!red->ActuallyReduced())
    {
      cerr << " read_image() Error: dbimage '" 
	   << words.front() << "' was not reduced " << endl;
      delete red;
      return 0;
    }
  return red;
}

// driver routine to read the objects and dbimages in a lightcurve file. 
// so far kept local to this program
static istream& read_lcfile(istream& stream, 
			    RefStarList &objects, 
			    ReducedImageList  &images,
			    ReducedImage* &phoref)
{

  bool objin=false, rimin=false, phoin=false;
  string line;
 
  while (getline(stream, line))
    {
      if (line[0] == '#' || line.length() == 0) continue; 
      
      if (line.find("OBJECTS") == 0)
	{ rimin=false; objin=true; phoin=false;  continue; }
      
      if (line.find("IMAGES") == 0)   
	{ 
	  //cout << "Read key IMAGES" << endl;
	  rimin=true;  objin=false; phoin=false; continue; }

      if (line.find("PHOREF") == 0)   
	{ rimin=false;  objin=false; phoin=true; continue; }

      if (objin) 
	{ RefStar *star = read_object(line); if (star) objects.push_back(star); continue; }

      if (rimin) 
	{ ReducedImage *im = read_image(line); if (im) images.push_back(im); continue; }

      if (phoin) 
	{ ReducedImage *im = read_image(line); if (im) phoref = im; continue; }
    }

  return stream;
}

istream& lc_read(istream& Stream, RefStarList& Objects, ReducedImageList& Images)
{
  if (!Stream) return Stream;

  ReducedImage *ref = 0;
  read_lcfile(Stream, Objects, Images, ref);
  cout << " Number of images = " << Images.size() << endl;


  if (ref) for_each(Objects.begin(), Objects.end(),
  		    bind2nd(mem_fun(&Fiducial<PhotStar>::AssignImage), ref));
  
  if (Objects.size() == 0) 
    {
      cerr << " lc_read() : Error: no valid objects " << endl;
      Stream.setstate(ios::badbit);
    }

  if (Images.size() < 2)
    {
      cerr << " lc_read() : Error: need at least 2 valid images " << endl;
      Stream.setstate(ios::badbit);
    }
  
  if (Stream) cout << " lc_read() : Read in successfully " 
		   << Objects.size() << " object(s) and " 
		   << Images.size()  << " images\n";

#ifdef DEBUG
  cout << " lcread() : debug : loaded the following lightfile: \n";
  lc_write(cout, Objects, Images);
#endif


  return Stream;
}

ostream& lc_write(ostream& Stream, const RefStarList& Objects, const ReducedImageList& Images)
{
  Stream << "OBJECTS\n";
  
  for(RefStarCIterator it = Objects.begin(); it != Objects.end(); ++it)
    {
      Stream << (*it)->x << " " 
	     << (*it)->y << " "
	     << " DATE_MIN=" << (*it)->jdmin
	     << " DATE_MAX=" << (*it)->jdmax
	     << " NAME="     << (*it)->name << endl;      
    }

  Stream << "IMAGES\n";
  for(ReducedImageCIterator it = Images.begin(); it != Images.end(); ++it) 
    Stream << (*it)->Name() << endl;
  
  Stream << "PHOREF\n"  << Objects.front()->Image()->Name() << endl;

  return Stream;
}

