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
#include <imageutils.h>
#include <apersestar.h>
#include <fastfinder.h>
#include <string>
#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <dbimage1> <dbimage2> <dbimage3> ... -r <referenceimage> -c <catalog>" << endl ;
  cerr << "options:"<< endl;
  cerr << "     -o <catalog> : output catalog name (default is calibration.list)" << endl;
  cerr << "     -n # : max number of images (default is unlimited)" << endl;
  cerr << "     -f # : first star to fit (def is 1, starts at 1)" << endl;
  cerr << "     -l # : last star to fit (def is 1000, included)" << endl;
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
  double neighborDist,neighborFlux;
  double neighborFluxContamination;
  double neighborNsigma;
  int id;
};

static double sqr(const double& x) {return x*x;};

int main(int argc, char **argv)
{
  string referencedbimage = "";
  string catalogname = "";
  string matchedcatalogname = "calibration.list";
  vector<string> dbimages;
  size_t maxnimages = 0;
  int first_star = 1;
  int last_star  = 1000;
  if (argc < 7)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  dbimages.push_back(arg);continue;
	}
      switch (arg[1])
	{
	case 'f' : first_star=atoi(argv[++i]); break;
	case 'l' : last_star=atoi(argv[++i]); break;
	case 'r' : referencedbimage = argv[++i]; break;
	case 'c' : catalogname = argv[++i]; break;
	case 'o' : matchedcatalogname = argv[++i]; break;
	case 'n' : maxnimages = atoi(argv[++i]); break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  

  if(!FileExists(catalogname)) {
    cerr << "cant find catalog " << catalogname << endl;
    usage(argv[0]);
  }
  cout << "catalog          = " << catalogname << endl;
  cout << "referencedbimage = " << referencedbimage << endl;
  cout << "n. dbimages      = " << dbimages.size();
  if( maxnimages > 0 && dbimages.size() > maxnimages) cout << " limited to " << maxnimages;
  cout << endl;

  // put all of this info in a LightCurveList which is the food of the photometric fitter
  LightCurveList lclist;
  lclist.RefImage = new ReducedImage(referencedbimage);
  for (size_t im=0;im<dbimages.size();++im) {
    if ( maxnimages > 0 && im >=  maxnimages ) break;
    lclist.Images.push_back(new ReducedImage(dbimages[im]));
  }
  
  // we know want to put new objects in the list
  lclist.Objects.clear();
  lclist.clear();

  
  FitsHeader header(lclist.RefImage->FitsName());
  // prepare transfo and frames for stars' selection
  Frame W = lclist.RefImage->UsablePart();
  // W = W.Rescale(1.); // remove boundaries
  Gtransfo* Pix2RaDec=0;
  WCSFromHeader(header, Pix2RaDec);
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01,W);
  Frame radecW = (ApplyTransfo(W,*Pix2RaDec)).Rescale(1.1);

  DictFile catalog(catalogname);
  
  //int requiredlevel=2;
  // get keys for mag
  string band = header.KeyVal("TOADBAND");
  string mag_key=getkey("m"+band,catalog);
  
  double mag_limit = 99;
  
  if(band=="g") mag_limit = 21.;
  if(band=="r") mag_limit = 21.;
  if(band=="i") mag_limit = 21.;
  if(band=="z") mag_limit = 21.;
  

  BaseStar star;
  int count_total=0;
  int count_total_stars=0;
  int count_ok=0;
  
  double mag;
  char name[100];
  
  map<RefStar*,CalibratedStar> assocs;
  
  // open catalog for match
  string refimage_catalog = lclist.RefImage->AperCatalogName();
  if( !FileExists(refimage_catalog)) {
    cerr << "need aperse catalog of reference image to calibrate" << endl;
    exit(EXIT_FAILURE);
  }
     
  AperSEStarList starlist(refimage_catalog);
  FastFinder finder((const BaseStarList&)starlist);
  
  for(DictFileCIterator entry=catalog.begin();entry!=catalog.end();++entry) {

    count_total++;
    
    //if(int(entry->Value("level"))<requiredlevel)  continue; // not a star with correct level 
    mag=entry->Value(mag_key);

    

    count_total_stars++;
    
    if(mag>mag_limit) continue; // ignore dim stars unused for calibration to save CPU


    star.x=entry->Value("ra"); //star.x=entry->Value("x"); // ra (deg)
    star.y=entry->Value("dec");//star.y=entry->Value("y"); // dec (deg)
    
    // now check if star is in image
    if(!radecW.InFrame(star)) continue; // bye bye
    // apply transfo to this star (x,y)=pixels
    RaDec2Pix->apply(star,star);
    // check again
    if (!W.InFrame(star)) continue; // bye bye
    
    
    // now do the match to get better coordinates
    const BaseStar * closest_basestar = NULL;  
    const AperSEStar * second_closest_star 
      = dynamic_cast<const AperSEStar *>(finder.SecondClosest(star,50.,closest_basestar)); // 10 pixels
    if ( !closest_basestar) {
      cerr << "ERROR cannot find star at " << star.x << " " << star.y << endl;
      continue;
    }
    count_ok++;
    
    if(count_ok<first_star || count_ok>last_star) {
      cout << "warning, skipping star number " << count_ok << endl;
      continue;
    }

    // mod coordinates
    star.x = closest_basestar->x;
    star.y = closest_basestar->y;
    
    // compute rough flux contamination
    const AperSEStar * closest_star = dynamic_cast<const AperSEStar *>(closest_basestar);
    double dist = 0;
    double max_flux_contamination = 0;
    double nsigma = 1.e6;
    
    // based on D2 ccd_35 : average is 2.378, 90%<3.1
    double calibration_seeing = 3.1;
    double calibration_aperture_radius = 5.*calibration_seeing;
    double reference_seeing = lclist.RefImage->Seeing();
    double seeing_scale = sqr(calibration_seeing/reference_seeing);
    
    if ( second_closest_star ) {
      
      const double &aper_radius = calibration_aperture_radius;
      double mxx = second_closest_star->gmxx;
      double myy = second_closest_star->gmyy;
      double mxy = second_closest_star->gmxy;

      // scale with best image seeing
      mxx *= seeing_scale;
      myy *= seeing_scale;
      mxy *= seeing_scale;
      
      
      Point dP =  *second_closest_star - *closest_basestar;
      dist = sqrt(dP.x*dP.x+dP.y*dP.y);
      
      nsigma = 0;
      bool is_inside = (dist<aper_radius);
      double distance_scale = (dist - aper_radius)/dist;
      if(dist<aper_radius) {
	//cout << "whao , second_nearest position is inside aperture radius" << endl;
	distance_scale *= -1;
      }
	
      dP.x *= distance_scale;
      dP.y *= distance_scale;
      
      // calcul du moment de l'objet   
      double det = mxx*myy-mxy*mxy;
      nsigma = sqrt((myy*dP.x*dP.x-2.*mxy*dP.x*dP.y+mxx*dP.y*dP.y)/det);
      // for gaussian
      double frac_gaus= erfc(nsigma/sqrt(2.))/2.;
      
      // for a moffat 2.5
      // ------------------------------------------
      double nsigma_moffat = nsigma*0.6547; // inside ~4 sigma radius 
      double nsigma2 = nsigma_moffat*nsigma_moffat;
      double frac_moffat = 0.1591549431*(-2*nsigma_moffat+(3.141592654-2*atan(nsigma_moffat))*(1.+nsigma2))/(1.+nsigma2);
      
      if(is_inside) {
	frac_gaus   = 1.-frac_gaus;
	frac_moffat = 1.-frac_moffat;
      }

      if(frac_moffat>frac_gaus)
	max_flux_contamination = frac_moffat*second_closest_star->flux;
      else
	max_flux_contamination = frac_gaus*second_closest_star->flux;
      
      cout << "x,y,dist,ndist,mxx,mxy,myy,nsigma,fg,fm,rfrac="
	   << closest_star->x << " "
	   << closest_star->y << " " 
	   << dist << " " 
	   << dist - aper_radius << " " 	
	   << mxx  << " " 
	   << mxy  << " " 
	   << myy  << " " 
	   << nsigma  << " " 
	   << frac_gaus << " " 
	   << frac_moffat << " " 
	   << max_flux_contamination/closest_star->flux << " "
	   << endl;
    }    
    
    // set flux
    star.flux=pow(10.,-0.4*mag);
    
    // ok now we copy this star in a refstar and put it in the lclist.Objects
    CountedRef<RefStar> rstar = new RefStar(lclist.RefImage);
    sprintf(name,"calibstar%d",count_ok);
    rstar->name = name;
    rstar->type = 1; // a star
    rstar->band = band[0];
    rstar->x = star.x;
    rstar->y = star.y;
    rstar->ra = entry->Value("ra"); // entry->Value("x"); // ra (deg)
    rstar->dec = entry->Value("dec"); // entry->Value("y"); // dec (deg)
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
    
    cstar.ra=entry->Value("ra"); // entry->Value("x");
    cstar.dec=entry->Value("dec"); // entry->Value("y");
    
    if(entry->HasKey("mu")) cstar.u=entry->Value("mu"); else cstar.u=99;
    if(entry->HasKey("mg")) cstar.g=entry->Value("mg"); else cstar.g=99;
    if(entry->HasKey("mr")) cstar.r=entry->Value("mr"); else cstar.r=99;
    if(entry->HasKey("mi")) cstar.i=entry->Value("mi"); else cstar.i=99;
    if(entry->HasKey("mz")) cstar.z=entry->Value("mz"); else cstar.z=99;
    if(entry->HasKey("emu")) cstar.ue=entry->Value("emu"); else cstar.ue=99;
    if(entry->HasKey("emg")) cstar.ge=entry->Value("emg"); else cstar.ge=99;
    if(entry->HasKey("emr")) cstar.re=entry->Value("emr"); else cstar.re=99;
    if(entry->HasKey("emi")) cstar.ie=entry->Value("emi"); else cstar.ie=99;
    if(entry->HasKey("emz")) cstar.ze=entry->Value("emz"); else cstar.ze=99;
    
    cstar.flux=star.flux;
    cstar.id=count_total;
    cstar.x=star.x;
    cstar.y=star.y;
    if(second_closest_star) {
      cstar.neighborDist = dist;
      cstar.neighborFlux = second_closest_star->flux;
    }else{
      cstar.neighborDist = closest_star->neighborDist;
      cstar.neighborFlux = closest_star->neighborFlux;
    }
    cstar.neighborFluxContamination = max_flux_contamination;
    cstar.neighborNsigma = nsigma;

    //keep a link between the two of them
    assocs[(RefStar*)rstar] = cstar;
    //if(count_ok>0)
    //break;

  } // FIN BOUCLE FOR

  
  cout << count_total << " objects in the catalog" << endl;
  cout << count_total_stars << " correct stars (with mag in band " << band << ") in the catalog" << endl;
  cout << count_ok << " stars in this image" << endl;
  
  //exit(0); // DEBUG
  
  // ok now let's do the fit
  SimFitPhot doFit(lclist,false);
  doFit.bWriteVignets=false; // don't write anything before all is done
  doFit.bWriteLC=false;
  
  // does everything
  // for_each(lclist.begin(), lclist.end(), doFit);
  
  // now we want to write many many things, let's make a list
  ofstream stream(matchedcatalogname.c_str());
  stream << "@CALIBCATALOG " << catalogname << endl;
  stream << "@NSTARS " << count_ok << endl;
  stream << "@NIMAGES " << lclist.Images.size() << endl;
  stream << "#x :" << endl;
  stream << "#y :" << endl;
  stream << "#flux :" << endl;
  stream << "#error :" << endl;
  stream << "#sky :" << endl;
  stream << "#skyerror :" << endl;
  stream << "#xerror :" << endl;
  stream << "#yerror :" << endl;
  stream << "#name :" << endl;
  stream << "#mjd :" << endl;
  stream << "#seeing :" << endl;
  stream << "#exptime :" << endl;
  stream << "#phratio :" << endl;
  stream << "#gseeing :" << endl;
  stream << "#sesky :" << endl;
  stream << "#sigsky :" << endl;
  stream << "#sigscale :" << endl; 
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
  stream << "#chi2pdf : chi2 pdf of total PSF photometry" << endl;
  stream << "#satur : 1 if some pixels are saturated" << endl;
  stream << "#nsatur : number of pixels  saturated" << endl;
  stream << "#neid : distance to nearest neighbor" << endl;
  stream << "#neif : flux of nearest neighbor" << endl;
  stream << "#neic : flux contamination due to nearest neighbor" << endl;
  stream << "#neins : nsigma of nearest neighbor" << endl;
  stream << "#end" <<endl;
  stream << setprecision(12);
  
  
  /* reserve_images
     first call simfit::load with only_reserve_images = true
     which calls simfitvignet::prepareautoresize
     which in turn calls  reserve_dimage_in_server(...)

     after, with the std call to simfit::load with only_reserve_images = false
     vignet::resize calls vignet::load and there get_dimage_from_server(...)

   */
  for(LightCurveList::iterator ilc = lclist.begin(); ilc!= lclist.end() ; ++ilc) { // loop on lc
    // each entry is a star
    doFit.zeFit.Load(*ilc,false,true);
  }


  
  for(LightCurveList::iterator ilc = lclist.begin(); ilc!= lclist.end() ; ++ilc) { // loop on lc
    
    doFit(*ilc);
    
    CalibratedStar cstar=assocs.find(ilc->Ref)->second;
    //cout << "=== " << cstar.r << " " << cstar.flux << " ===" << endl;
    int count_img=0;
    double chi2pdf=ilc->chi2ndf();
    for (LightCurve::const_iterator it = ilc->begin(); it != ilc->end(); ++it) { // loop on points
      count_img++;
      const Fiducial<PhotStar> *fs = *it;

      // ###########"" WWWWWWWWARNING
      // on l'enleve pour test 15/02/2010
      //if(fabs(fs->flux)<0.001) // do not print unfitted fluxes
      //continue;

      double sigposX = 0;
      double sigposY=0;

      string aligned_name = fs->Name() ; 
      size_t align_pos = aligned_name.find("enlarged");
      string dbim_name;
      if(align_pos!=string::npos) dbim_name = aligned_name.erase(0,align_pos+8); 
      size_t pos = dbim_name.find("p");

      if(pos!=string::npos) dbim_name.replace(pos,1,"");

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
      if(fs->varx>0)
	stream << sqrt(fs->varx) << " ";
      else
	stream << 0 << " ";
      if(fs->vary>0)
	stream << sqrt(fs->vary) << " ";
      else
	stream << 0 << " ";
      stream << dbim_name << " ";
      stream << fs->ModifiedJulianDate() << " ";
      stream << fs->Seeing() << " ";
      stream << fs->ExposureTime() << " ";
      stream << fs->photomratio << " ";
      stream << fs->GFSeeing() << " ";
      stream << fs->SESky() << " ";
      stream << fs->SIGSky() << " ";
      stream << fs->sigscale_varflux << " ";
      
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
      
      if(fs->has_saturated_pixels)
	stream << 1 << " ";  
      else
	stream << 0 << " ";  
      stream << fs->n_saturated_pixels << " ";
      stream << cstar.neighborDist << " "; 
      stream << cstar.neighborFlux << " "; 
      stream << cstar.neighborFluxContamination << " ";  
      stream << cstar.neighborNsigma << " ";        
      stream << endl;
    }
  }
  stream.close();
  return EXIT_SUCCESS;
}

