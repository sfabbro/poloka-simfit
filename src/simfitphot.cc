#include <algorithm> // min_element

#include <fastfinder.h>
#include <reducedutils.h>
#include <fstream>
#include "simfit.h"
#include "simfitphot.h"

SimFitPhot::SimFitPhot(LightCurveList& Fiducials,bool usegal)
{
#ifdef FNAME
  cout << " > SimFitPhot::SimFitPhot(LightCurveList& Fiducials)" << endl;
#endif
  //init_phot(Fiducials);
#ifdef DEBUG
  cout << "  zeFit.reserve ... " << endl;
#endif
  bWriteVignets=false;
  bWriteLC=true;
  bOutputDirectoryFromName=false;
  
  zeFit.reserve(Fiducials.Images.size());
  
  zeFit.VignetRef = new SimFitRefVignet(Fiducials.RefImage,usegal); //  no data will be read cause no star is defined
#ifdef DEBUG
  SimFitRefVignet* toto = zeFit.VignetRef;
  cout << " in SimFitPhot::SimFitPhot,  zeFit.VignetRef     = " << toto << endl;
#endif
  for (ReducedImageCIterator it=Fiducials.Images.begin(); it != Fiducials.Images.end(); ++it)
    {
      cout << flush << " > SimFitPhot::SimFitPhot() : making vignets for " << (*it)->Name() << "\r";
      SimFitVignet *vig = new SimFitVignet(*it,zeFit.VignetRef);
      zeFit.push_back(vig);
    }
  // cannot do this since we do not have any data loaded yet 
  //zeFit.FindMinimumScale(worstSeeing);
  
}

#define DEBUG0
void SimFitPhot::operator() (LightCurve& Lc)
{

#ifdef DEBUG0
  switch (Lc.Ref->type)
    {
    case -1:  cout << " ============= SimFitPhot::operator() Fitting a Sn + Galaxy with FIXED position \"" << Lc.Ref->name << "\" =============" << endl; break;
    case 0: cout << " ============= SimFitPhot::operator() Fitting a Sn + Galaxy \"" << Lc.Ref->name << "\" =============" << endl; break;
    case 1: cout << " ============= SimFitPhot::operator() Fitting a star (without galaxy) \"" << Lc.Ref->name << "\"  =============" << endl; break;
    case 3: cout << " ============= SimFitPhot::operator() Fitting a star (without galaxy) with fixed pos \"" << Lc.Ref->name << "\"  =============" << endl; break;
    case 2: cout << " ============= SimFitPhot::operator() Fitting a galaxy (without star) \"" <<  Lc.Ref->name << "\"  =============" << endl; break;
    default: cout << " SimFitPhot::operator() : Error : unknown star type :" << Lc.Ref->type << endl; return;
    }
#endif

#ifdef DEBUG0
  cout << " ============= SimFitPhot::operator() zeFit.Load =============" << endl;
#endif  
  zeFit.Load(Lc);
  
  //string dir = "./lc";
  string dir = ".";
  if(bOutputDirectoryFromName) {
    dir = Lc.Ref->name;
  }
#ifdef DEBUG0
  cout << "DEBUG bOutputDirectoryFromName " << bOutputDirectoryFromName << endl;
  cout << "DEBUG bWriteVignets " << bWriteVignets << endl;
  cout << "DEBUG bWriteLC " << bWriteLC << endl;
  cout << "DEBUG isdir  " << dir << " " << IsDirectory(dir) << endl;
#endif  
  if( bOutputDirectoryFromName && (bWriteVignets || bWriteLC) && !IsDirectory(dir))
    MKDir(dir.c_str());
  
  
  //============================================================
  // galaxy 
  //============================================================
  if(Lc.Ref->type == 2) {
    zeFit.SetWhatToFit(FitGal);
    zeFit.UseGalaxyModel(true);
    if(! zeFit.DoTheFit(50,0.005)) {
      return;
    }
  }

  //============================================================
  // star with galaxy 
  //============================================================
  if(Lc.Ref->type == 0) {
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() First FitFlux | FitGal | FitSky =============" << endl;
#endif	
    zeFit.SetWhatToFit(FitFlux | FitGal | FitSky); 
    if(! zeFit.DoTheFit(0,0.1)) return;  
    zeFit.write("sn_init",dir, WriteGalaxy);
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() First FitFlux | FitPos  =============" << endl;
#endif
    zeFit.SetWhatToFit(FitFlux | FitPos);
    if(! zeFit.DoTheFit(3,0.1)) return;
    
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() Now FitFlux | FitPos | FitGal | FitSky =============" << endl;
#endif	
    zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky); // then everything    
    if(! zeFit.DoTheFit(30,0.1)) return;     
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
       (*itVig)->CheckWeight();
     }
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitPos | FitGal | FitSky =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitPos | FitSky);
     if(! zeFit.DoTheFit(30,0.05)) return;
  }


  //============================================================
  // star with galaxy BUT with fixed position
  //============================================================
  if(Lc.Ref->type == -1) {
#ifdef DEBUG0
    cout << " ============= SimFitPhot::operator() FitInitialGalaxy =============" << endl;
#endif
     zeFit.SetWhatToFit(FitFlux | FitGal | FitSky); // then everything    
     if(! zeFit.DoTheFit(30,0.05)) return;     
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() Robustify  =============" << endl;
#endif
     for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
       (*itVig)->KillOutliers();
       (*itVig)->CheckWeight();
     }
#ifdef DEBUG0
     cout << " ============= SimFitPhot::operator() refit FitFlux | FitGal  | FitSky =============" << endl;
#endif	
     zeFit.SetWhatToFit(FitFlux | FitGal | FitSky);
     if(! zeFit.DoTheFit(30,0.005)) return;
  }
  
  //============================================================
  // star without galaxy
  //============================================================
  if(Lc.Ref->type==1){
    zeFit.SetWhatToFit(FitFlux);
    zeFit.UseGalaxyModel(false);
    if(! zeFit.DoTheFit()) return;
    zeFit.SetWhatToFit(FitFlux  | FitPos | FitSky );
    zeFit.UseGalaxyModel(false);
    if(! zeFit.DoTheFit(10,1)) return;
    // robustify to get rid of other stars in the vignet
    
    for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
      (*itVig)->KillOutliers();
      (*itVig)->CheckWeight();
    }
    zeFit.SetWhatToFit(FitFlux  | FitPos | FitSky );
    zeFit.UseGalaxyModel(false);
    if(! zeFit.DoTheFit(10,0.5)) return;
  }
  //============================================================
  // star without galaxy with fixed pos
  //============================================================
  if(Lc.Ref->type==3){
    zeFit.SetWhatToFit(FitFlux);
    zeFit.UseGalaxyModel(false);
    if(! zeFit.DoTheFit()) return;
    zeFit.SetWhatToFit(FitFlux  | FitSky );
    zeFit.UseGalaxyModel(false);
    if(! zeFit.DoTheFit(10,0.05)) return;
    // robustify to get rid of other stars in the vignet
    for (SimFitVignetIterator itVig = zeFit.begin(); itVig != zeFit.end(); ++itVig) {
      (*itVig)->KillOutliers();
      (*itVig)->CheckWeight();
    }
    zeFit.SetWhatToFit(FitFlux  | FitSky );
    zeFit.UseGalaxyModel(false);
       if(! zeFit.DoTheFit(10,0.05)) return;
    zeFit.GetCovariance();
  }

  // assign results
  Lc.totflux    = zeFit.TotFlux();
  Lc.totsky     = zeFit.TotSky();
  Lc.galflux    = zeFit.GalFlux();
  Lc.vargalflux = zeFit.VarGalFlux();
  Lc.vartotflux = zeFit.VarTotFlux();
  Lc.vartotsky  = zeFit.VarTotSky();
  Lc.chi2       = zeFit.Chi2();
  Lc.ndf        = zeFit.Dof();
  Lc.resmean    = zeFit.MeanMedianRms(Lc.resmed, Lc.resrms, Lc.resadev);

  
 
  if(bWriteVignets) {
    zeFit.write("sn",dir, WriteGalaxy|WriteResid|WriteWeight|WriteData|WritePsf);
  }

  
  if(bWriteLC) {
    zeFit.write("sn",dir,WriteLightCurve|WriteVignetsInfo|WriteMatrices);
    ofstream lstream((string(dir+"/lc2fit.dat")).c_str());
    Lc.write_lc2fit(lstream);
  }

}



