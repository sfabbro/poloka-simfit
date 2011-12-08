#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dictfile.h>
#include <basestar.h>
#include <frame.h>
#include <wcsutils.h>
#include <reducedimage.h>
#include <fitsimage.h>
#include <gtransfo.h>
#include <photstar.h>
#include <lightcurve.h>
#include <simfitphot.h>
#include <vutils.h>
#include <imagematch.h>
#include <gtransfo.h>

#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <dbimage1> <dbimage2> <dbimage3> ... -r <referenceimage> -p <catalog4positions> <referenceimageforpositions>" << endl ;
  cerr << "options:"<< endl;
  cerr << "     -o <catalog> : output catalog name (default is calibration.list)" << endl;
  cerr << "     -n # : max number of images (default is unlimited)" << endl;
  exit(1);
}

static string getkey(const string& prefix, const DictFile& catalog) {
  if(!catalog.HasKey(prefix)) {
    cerr << "catalog does not not have info about " << prefix << endl;
    exit(2);
    return "absent";
  }
  return prefix;
}

class CalibratedStar : public BaseStar {
 public:
  CalibratedStar() {};
  CalibratedStar(BaseStar& toto) : BaseStar(toto) {};
  double ra,dec;
  double u,g,r,i,z,x,y;
  double ue,ge,re,ie,ze;
  int id;
};

int main(int argc, char **argv)
{
  string referencedbimage = "";
  string matchedcatalogname = "calibration_fixpos.list";
  vector<string> dbimages;
  size_t maxnimages = 0;
  string referencedbimage_for_positions = "";
  string catalognamefor_positions = "";
  
  if (argc < 9)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  dbimages.push_back(arg);continue;
	}
      switch (arg[1])
	{
	case 'r' : referencedbimage = argv[++i]; break;
	case 'o' : matchedcatalogname = argv[++i]; break;
	case 'n' : maxnimages = atoi(argv[++i]); break;
	case 'p' : catalognamefor_positions = argv[++i]; referencedbimage_for_positions = argv[++i]; break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  
  cout << "catalog          = " << catalognamefor_positions << endl;
  cout << "referencedbimage = " << referencedbimage << endl;
  cout << "referencedbimage_for_positions = " << referencedbimage_for_positions << endl;  
  cout << "n. dbimages      = " << dbimages.size();
  if( maxnimages > 0 && dbimages.size() > maxnimages) cout << " limited to " << maxnimages;
  cout << endl;
  
  CountedRef<Gtransfo> direct,reverse;
  // find transfo from referencedbimage_for_positions to referencedbimage
  {
    ReducedImage image_from(referencedbimage_for_positions);
    ReducedImage image_to(referencedbimage); 
    ImageListMatch(image_from,image_to,direct, reverse);
  }
  
  // put all of this info in a LightCurveList which is the food of the photometric fitter
  LightCurveList lclist;
  lclist.RefImage = new ReducedImage(referencedbimage);
  for (size_t im=0;im<dbimages.size();++im) {
    if ( maxnimages > 0 && im >=  maxnimages ) break;
    lclist.Images.push_back(new ReducedImage(dbimages[im]));
  }

  FitsHeader header(lclist.RefImage->FitsName());
  // frame of reference image
  Frame W = lclist.RefImage->UsablePart();

  // we know want to put new objects in the list
  lclist.Objects.clear();
  lclist.clear();

  DictFile catalog(catalognamefor_positions);
  
  // get keys for mag
  string band = header.KeyVal("TOADBAND");
  string mag_key=getkey(band,catalog);
  
  BaseStar star;
  int count_total=0;
  int count_total_stars=0;
  int count_ok=0;
  
  double mag;
  char name[100];
  //double mag_med;
  //double mag_rms;
  map<RefStar*,CalibratedStar> assocs;
  int img;
  for(DictFileCIterator entry=catalog.begin();entry!=catalog.end();++entry) {

    img=entry->Value("img");
    if(img>1) continue;
    
    count_total++;
    
    mag=entry->Value(mag_key);
    
    count_total_stars++;
    
    star.x=entry->Value("x"); // x of ref
    star.y=entry->Value("y"); // y of ref
    
    // apply transfo to this star (x,y)=pixels in image
    direct->apply(star,star);
    // check if in image
    if (!W.InFrame(star)) continue; // bye bye
    
    count_ok++;
    
    // set flux
    star.flux=pow(10.,-0.4*mag);
    
    // ok now we copy this star in a refstar and put it in the lclist.Objects
    CountedRef<RefStar> rstar = new RefStar(lclist.RefImage);
    sprintf(name,"calibstar%d",count_ok);
    rstar->name = name;
    rstar->type = 3; // a star with fixed position
    rstar->band = band[0];
    rstar->x = star.x;
    rstar->y = star.y;
    rstar->ra = entry->Value("x"); // ra (deg)
    rstar->dec = entry->Value("y"); // dec (deg)
    rstar->jdmin = -1.e30; // always bright
    rstar->jdmax = 1.e30;
    lclist.Objects.push_back(rstar);
    
    // and also creat a lc (something stupid in the design)
    LightCurve lc(rstar);
    for (ReducedImageCIterator im=lclist.Images.begin(); im != lclist.Images.end(); ++im) {
      lc.push_back(*im, new PhotStar(star)); // add one image and one PhotStar
    }
    lclist.push_back(lc);
    
    // we also want to keep calibration info
    CalibratedStar cstar(star);
    cstar.ra=entry->Value("ra");
    cstar.dec=entry->Value("dec");
    cstar.u=entry->Value("u");
    cstar.g=entry->Value("g");
    cstar.r=entry->Value("r");
    cstar.i=entry->Value("i");
    cstar.z=entry->Value("z");
    cstar.ue=entry->Value("ue");
    cstar.ge=entry->Value("ge");
    cstar.re=entry->Value("re");
    cstar.ie=entry->Value("ie");
    cstar.ze=entry->Value("ze");
    cstar.flux=star.flux;
    cstar.id=count_total;
    cstar.x=star.x;
    cstar.y=star.y;
    
    //keep a link between the two of them
    assocs[(RefStar*)rstar] = cstar;
    //if(count_ok>0)
    //break;
  }
  
  cout << count_total << " objects in the catalog" << endl;
  cout << count_total_stars << " correct stars (with mag in band " << band << ") in the catalog" << endl;
  cout << count_ok << " stars in this image" << endl;
  
  // ok now let's do the fit
  SimFitPhot doFit(lclist,false);
  doFit.bWriteVignets=false; // don't write anything before all is done
  doFit.bWriteLC=false;
  
  // does everything
  // for_each(lclist.begin(), lclist.end(), doFit);
  
  // now we want to write many many things, let's make a list
  ofstream stream(matchedcatalogname.c_str());
  stream << "@NSTARS " << count_ok << endl;
  stream << "@NIMAGES " << lclist.Images.size() << endl;
  stream << "#x :" << endl;
  stream << "#y :" << endl;
  stream << "#flux :" << endl;
  stream << "#error :" << endl;
  stream << "#sky :" << endl;
  stream << "#skyerror :" << endl;
  stream << "#mag :" << endl;
  stream << "#mage :" << endl;
  stream << "#ra : initial " << endl;
  stream << "#dec : initial " << endl;
  stream << "#ix : initial x" << endl;
  stream << "#iy : initial y" << endl;
  stream << "#u : from catalog" << endl;
  stream << "#g : from catalog" << endl;
  stream << "#r : from catalog" << endl;
  stream << "#i : from catalog" << endl;
  stream << "#z : from catalog" << endl;
  stream << "#ue : from catalog" << endl;
  stream << "#ge : from catalog" << endl;
  stream << "#re : from catalog" << endl;
  stream << "#ie : from catalog" << endl;
  stream << "#ze : from catalog" << endl;
  stream << "#img : image number" << endl;
  stream << "#star : start number in the catalog" << endl;
  stream << "#chi2pdf : chi2 of PSF photometry" << endl;
  stream << "#end" <<endl;
  stream << setprecision(12);
  
  
  for(LightCurveList::iterator ilc = lclist.begin(); ilc!= lclist.end() ; ++ilc) { // loop on lc

    doFit(*ilc);


    CalibratedStar cstar=assocs.find(ilc->Ref)->second;
    //cout << "=== " << cstar.r << " " << cstar.flux << " ===" << endl;
    int count_img=0;
    double chi2pdf=ilc->chi2ndf();
    for (LightCurve::const_iterator it = ilc->begin(); it != ilc->end(); ++it) { // loop on points
      count_img++;
      const Fiducial<PhotStar> *fs = *it;
      if(fabs(fs->flux)<0.001) // do not print unfitted fluxes
	continue;
      stream << fs->x << " ";
      stream << fs->y << " ";
      stream << fs->flux << " ";
      if(fs->varflux>0)
	stream << sqrt(fs->varflux) << " ";
      else
	stream << 0 << " ";
      stream << fs->sky << " ";
      if(fs->varsky>0 && fs->varsky<10000000.)
	stream << sqrt(fs->varsky) << " ";
      else
	stream << 0 << " ";
      
      // mag
      if (band=="u") stream << cstar.u << " " << cstar.ue << " ";
      if (band=="g") stream << cstar.g << " " << cstar.ge << " ";
      if (band=="r") stream << cstar.r << " " << cstar.re << " ";
      if (band=="i") stream << cstar.i << " " << cstar.ie << " ";
      if (band=="z") stream << cstar.z << " " << cstar.ze << " ";
      
      stream << cstar.ra << " ";
      stream << cstar.dec << " ";
      stream << cstar.x << " ";
      stream << cstar.y << " ";
      stream << cstar.u << " ";
      stream << cstar.g << " ";
      stream << cstar.r << " ";
      stream << cstar.i << " ";
      stream << cstar.z << " ";
      stream << cstar.ue << " ";
      stream << cstar.ge << " ";
      stream << cstar.re << " ";
      stream << cstar.ie << " ";
      stream << cstar.ze << " ";
      stream << count_img << " "; 
      stream << cstar.id << " "; 
      stream << chi2pdf << " ";      
      stream << endl;
    }
  }
  stream.close();
  return EXIT_SUCCESS;
}

