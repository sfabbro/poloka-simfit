#include <iostream>
#include <fstream>
#include <lightcurve.h>
#include <simfitphot.h>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <filename> " << endl ;
}

int main(int argc, char **argv)
{
  if (argc < 2)  {usage(argv[0]);  exit(1);}

  ifstream lightfile(argv[1]);
  if (!lightfile) return EXIT_FAILURE;

  // make list of lightcurves
  LightCurveList fids(lightfile);
  
  for (LightCurveList::iterator si=fids.begin(); si != fids.end(); ++si) {
    cout << "Number of lightcurves to do : " << si->size() << endl;
    for (LightCurve::iterator it=si->begin(); it!= si->end() ; ++it) {
      const PhotStar *star = (*it);
      cout << "make_lightcurve main : " << *star << " " << (*it)->Image()->Name() << endl;
    }
  }

  // freach lightcurve, do the fit
  for (LightCurveList::iterator si=fids.begin(); si != fids.end(); ++si) {
    LightCurve &lc = (*si);
    
    
    SimFit zeFit;
    zeFit.CreateAndLoad(lc);
 //    // reserve space for all vignets
//     zeFit.reserve(fids.Images.size());
  
//     // now create 4 fwhm maximum size vignettes
//     const double rad = 4.*2.3548;
    
//     // reference vignette
//     const PhotStar *star = lc.Ref;
//     zeFit.VignetRef = SimFitRefVignet(star, fids.RefImage, 
// 				      int(ceil(fids.RefImage->Seeing()*rad)));
    
//     // all other vignets
//     double worstSeeing = -1e29;
//     for (ReducedImageCIterator it=fids.Images.begin(); it != fids.Images.end(); ++it) {
//       const double curSeeing = (*it)->Seeing();
//       SimFitVignet *vig = new SimFitVignet(star, *it, fids.RefImage, int(ceil(curSeeing*rad)));
//       zeFit.push_back(vig);
//       if (worstSeeing < curSeeing) worstSeeing = curSeeing;
//     }
    
//     zeFit.FindMinimumScale(worstSeeing);
    
    
  }
  
  
  //SimFitPhot doFit(fids);
  //for_each(fids.begin(), fids.end(), doFit);
  
  return EXIT_SUCCESS;
}

