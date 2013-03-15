#include <math.h>
#include <iostream>
#include <fstream>

#include <poloka/matvect.h>
#include <poloka/polokaexception.h>
#include <poloka/fileutils.h>
#include <poloka/dictfile.h>

#include <poloka/lightcurvepoint.h>

static double sq(const double x) { return x*x;}


static void usage(const char *progname) {
  cerr << "Usage: " << progname << " [OPTION] DIRECTORY\n"
       << "Fit one light curve point per night of for a light curve DIRECTORY\n\n"
       << "    -1 : fit one single constant flux for all points\n\n";
  exit(EXIT_FAILURE);
}


int main(int argc, char **argv) {

  bool fitsingleflux = false;
  string lcdir;

  for (int i=1; i<argc; ++i) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      if (lcdir.empty()) {
	cerr << argv[0] << ": unexpected argument " << arg << endl;
	usage(argv[0]);
      }
      lcdir = arg;
      continue;
    }
    switch (arg[1]) {
    case '1': 
      fitsingleflux = true;
      break;
    default : 
      cerr << argv[0] << ": unknown option " << arg << endl;
      usage(argv[0]);
      break;
    }
  }

 try {
   
  Vect FluxVec;
  {
    Mat m;
    if (m.readFits(lcdir + "/vec_sn.fits") != 0) {
      cerr << argv[0] << ": error reading vec_sn.fits from " << lcdir << endl;
      return EXIT_FAILURE;
    }
    FluxVec = m;
  }

  unsigned int nflux = FluxVec.Size();
  
  // matrice de covariance
  Mat CovarianceMat; 
  Mat fluxCovarianceMat;
  Mat FluxCovarianceMat;
  if (FileExists(lcdir+"/flux_pmat_sn.fits")) {
    cout << argv[0] << ": getthing flux covariance matrix from simphot\n";
    fluxCovarianceMat.readFits(lcdir+"/flux_pmat_sn.fits");
    fluxCovarianceMat.Symmetrize("R");
    FluxCovarianceMat = fluxCovarianceMat;
  } else {
    if (!FileExists(lcdir+"/pmat_sn.fits")) {
      cerr << argv[0] << ": missing pmat_sn.fits in " << lcdir << endl;
      return EXIT_FAILURE; 
    }
    cout << argv[0] << ": getting flux covariance matrix from lc";
    CovarianceMat.readFits(lcdir+"/pmat_sn.fits");
    FluxCovarianceMat = CovarianceMat.SubBlock(0,nflux-1,0,nflux-1);
    FluxCovarianceMat.Symmetrize("L");
  }
  

  Mat A;
  if (!fitsingleflux) {
    A.readFits(lcdir+"/nightmat_sn.fits");
  } else {
    A.allocate(1,nflux);
    for(unsigned int i=0; i<nflux; ++i)
      A(0,i)=1;
  }

#ifdef DEBUG
  cout << "A before cleaning:"  << endl;
  cout << A << endl;

  cout << "FluxVec before cleaning:"  << endl;
  cout << FluxVec << endl;
  
  cout << "FluxCovarianceMat before cleaning:"  << endl;
  cout << FluxCovarianceMat << endl;
#endif 

  // ==== remove points without data ====
  for (unsigned int i=0; i< FluxVec.Size(); i++) {
    if (fabs(FluxVec(i))<1.e-30) {
      cout << argv[0] << ": removing " << i << " : ";
      cout << FluxCovarianceMat.SizeX() << " => ";
      FluxCovarianceMat = FluxCovarianceMat.WithoutRows(i,i);
      FluxCovarianceMat = FluxCovarianceMat.WithoutColumns(i,i);
      cout << FluxCovarianceMat.SizeX() << endl;
      Mat mFluxVec = FluxVec;
      FluxVec = mFluxVec.WithoutRows(i,i);
      A = A.WithoutRows(i,i);
      i--;
    }
  }
  
  Mat Abis = A; // we save a copy for output 

#ifdef DEBUG
  cout << "FluxVec after cleaning:"  << endl;
  cout << FluxVec << endl;
  cout << "A after cleaning:"  << endl;
  cout << A << endl;
  cout << "FluxCovarianceMat after cleaning:"  << endl;
  cout << FluxCovarianceMat << endl;
#endif
  
  vector<int> suppressedfluxes;
  Vect flux_per_night;
  Mat FluxPerNightCovMat;
  Mat AtWA;
  Mat FluxWeightMat;

  int noutliers = 0;

  while (true) {
  
    FluxWeightMat = FluxCovarianceMat;
    FluxWeightMat.CholeskyInvert("L");
    //cout << "FluxWeightMat" << endl;
    //cout << FluxWeightMat << endl;
    AtWA = A.transposed()*FluxWeightMat*A;
    FluxPerNightCovMat= AtWA; 
    FluxPerNightCovMat.CholeskyInvert("L");
    flux_per_night = (FluxPerNightCovMat*(A.transposed()*FluxWeightMat))*FluxVec;
#ifdef DEBUG
    cout << "Mean flux per night" << endl;
    cout << flux_per_night << endl;
    cout << "Covariance matrix" << endl;
    cout <<FluxPerNightCovMat<< endl;
#endif
    // now compute chi2
    Vect B = FluxVec - A*flux_per_night;

    double chi2 = B*(FluxWeightMat*B);     
    int ndf = A.SizeY()-A.SizeX();

#ifdef DEBUG
    cout << "chi2 = " << chi2 << endl;   
    cout << "ndf = " << ndf << endl;
    cout << "chi2/ndf = " << chi2/ndf << endl;
#endif
   
    if (ndf==0 || chi2/ndf < 1.5 || suppressedfluxes.size() >= 6) {
      if (ndf == 0) {
	cout << argv[0] << ": degrees of freedom is zero\n";
      }
      break;
    }
    
    int outlier = -1;
    double flux_chi2;
    double chi2_max = 0;

    for (unsigned int iflux = 0; iflux<B.Size(); ++iflux) {      
      flux_chi2 = sq(B(iflux))*FluxWeightMat(iflux,iflux);
      if(flux_chi2>chi2_max) {
	chi2_max = flux_chi2;
	outlier = iflux;
      }
    }
    
    cout << outlier << " " << sqrt(chi2_max) << endl;
    
    if (sqrt(chi2_max)<3.)
      break;
    noutliers ++;
    suppressedfluxes.push_back(outlier);
    // on vire cet outlier de toutes les matrices
    A = A.WithoutRows(outlier,outlier);
    {
      Mat mFluxVec = FluxVec;
      FluxVec = mFluxVec.WithoutRows(outlier,outlier);
      //cout << FluxVec << endl;
    }
    FluxCovarianceMat = FluxCovarianceMat.WithoutRows(outlier,outlier);
    FluxCovarianceMat = FluxCovarianceMat.WithoutColumns(outlier,outlier);
  }
  
  // OUTPUT
  // ============================================================================================
  
  // compute chi2
  Vect B = FluxVec - A*flux_per_night;
  double chi2 = B*(FluxWeightMat*B);     
  int ndf = A.SizeY() - A.SizeX();
  double chi2ndf = 0;  
  if (ndf>0) chi2ndf = chi2/ndf;

  cout << argv[0] << ": SUMMARY "
       << FluxVec.Size() << " " 
       << flux_per_night.Size() << " "
       << noutliers << " " 
       << chi2ndf << endl;
  
  // if chi2dof>1 scale all errors
  
  if (chi2ndf>1) {
    for (unsigned int j= 0; j<FluxPerNightCovMat.SizeY();j++)
      for (unsigned int i= 0; i<FluxPerNightCovMat.SizeX();i++) {
	FluxPerNightCovMat(i,j) *= chi2ndf;
	AtWA(i,j) /= chi2ndf;
      }
  }
  
  // save these results in ASCII files
  // ... TODO ...
  {
    {
      ofstream st((lcdir+"/flux_per_expo_covmat.dat").c_str());
      st.setf(ios::fixed);
      st << FluxCovarianceMat;
      st.close();
    }
    FluxCovarianceMat.writeFits(lcdir+"/flux_per_expo_covmat.fits");
    
    {
      ofstream st((lcdir+"/flux_per_night_covmat.dat").c_str());
      st.setf(ios::fixed);
      st <<  FluxPerNightCovMat;
      st.close();
    }
    FluxPerNightCovMat.writeFits(lcdir+"/flux_per_night_covmat.fits");
    
    {
      ofstream st((lcdir+"/flux_per_expo_weightmat.dat").c_str());
      st.setf(ios::fixed);
      st << FluxWeightMat;
      st.close();
    }
    FluxWeightMat.writeFits(lcdir+"/flux_per_expo_weightmat.fits");
    
    {
      ofstream st((lcdir+"/flux_per_night_weightmat.dat").c_str());
      st.setf(ios::fixed);
      st <<  AtWA;
      st.close();
    }
    AtWA.writeFits(lcdir+"/flux_per_night_weightmat.fits");
    
  }
  
  /*
  // now we read the lightcurve point list to get julian dates and zero point
  vector< CountedRef<LightCurvePoint> > lcpoints;
  {
    vector< CountedRef<LightCurvePoint> > lcpoints_all;
    obj_input<xmlstream> oi("lc.xml");
    oi >> lcpoints_all;
    oi.close();
    
    // no keep only points with flux!=0 (i.e. fitted points)
    vector< CountedRef<LightCurvePoint> >::iterator it = lcpoints_all.begin();
    for(;it!=lcpoints_all.end();++it) {
      //cout << (*it)->julianday << " " << (*it)->flux;
      if((*it)->flux!=0) {
	lcpoints.push_back(*it);
      }
      //cout << endl;
    }
  }
  */
  DictFile lcdata(lcdir+"/lc2fit.dat");
  string instrumentName = lcdata.GlobalValue("INSTRUMENT");
  string bandName = lcdata.GlobalValue("BAND");
  string magSystem = lcdata.GlobalValue("MAGSYS");
  
  
  // get zero point 
  double zp = lcdata.front().Value("ZP");
  //cout << "zp=" << zp << endl;
  
  // read light curve points
  vector< CountedRef<LightCurvePoint> > lcpoints;
  for (DictFileCIterator line = lcdata.begin(); line != lcdata.end(); ++line) {
    CountedRef<LightCurvePoint> lcp = new LightCurvePoint();
    lcp->modifiedjulianday = line->Value("Date");
    lcp->flux = line->Value("Flux");
    lcpoints.push_back(lcp);
  }

  ofstream outputlc((lcdir+"/lc2fit_per_night.dat").c_str());
  outputlc << "#Date : (Modified julian date! days since January 1st, 2003)\n"
	   << "#Flux : \n"  
	   << "#Fluxerr : \n"
	   << "#ZP : elixir zp\n"
	   << "#chi2ndf : \n";
  outputlc << "@INSTRUMENT " << instrumentName << endl;
  outputlc << "@BAND " << bandName << endl;
  outputlc << "@MAGSYS " << magSystem << endl;
  outputlc.setf(ios::fixed);
  
  
  vector< CountedRef<LightCurvePoint> > newlcpoints;
  double mjd;
  int nexpo;
  for (unsigned int night = 0; night < flux_per_night.Size(); ++ night) {
    CountedRef<LightCurvePoint> newpoint = new LightCurvePoint();
    newpoint->flux = flux_per_night(night);
    newpoint->eflux = sqrt(FluxPerNightCovMat(night,night));
    newpoint->computemag(zp);
    // now get julian day (mean of all exposures)
    mjd=0;
    nexpo=0;

    if (Abis.SizeY() != lcpoints.size()) {
      char message[1000];
      sprintf(message,"not same number of exposures nightmat=%d lc=%d",int(Abis.SizeY()),int(lcpoints.size()));
      throw(PolokaException(message));
    }


    for (unsigned int expo=0; expo < Abis.SizeY(); ++expo) {
      if (Abis(night,expo) > 0.5) {
	// VERIFIER QUE ABIS.SIZEY est = LC.POINTS.SIZE
	//cout << "A.sizeY:" << A.SizeY() << "  Abis.sizeY:" << Abis.SizeY() << "   lcpoints.size:" << lcpoints.size() << "   newlcpoints.size:" << newlcpoints.size() << endl;
	mjd += lcpoints[expo]->modifiedjulianday;
	//cout << night << " " 
	//     << expo << " " 
	//     << lcpoints[expo]->julianday-2452854.0 << endl;
	nexpo++;
      }
    }
    newpoint->modifiedjulianday = mjd/nexpo;
    
    outputlc << (*newpoint) << " " << chi2ndf << endl;
    newlcpoints.push_back(newpoint);
  }
  outputlc.close();
 } catch (PolokaException p) {
   p.PrintMessage(cerr);
   return EXIT_FAILURE; 
 }
 return EXIT_SUCCESS;
}

