#include <iostream>
#include <fstream>
#include <algorithm>

#include <poloka/lightcurve.h>
#include <poloka/simfitphot.h>


using namespace std;

static void usage(const char *pgname) {

  cerr << pgname << " <filename> " << endl ;
  std::cerr << "-d : create a subdirectory per object" << std::endl; 
  std::cerr << "-v : write all vignets" << std::endl; 
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

  if (argc < 2) 
    usage(argv[0]);
  
  string lightfilename = "";
  bool subdirperobject = false;
  bool WriteVignets = false;
  for (int i=1; i<argc; ++i) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      if(lightfilename!="") {
	cerr << "unexpected argument " << arg << endl;
	usage(argv[0]);
      }
      lightfilename = arg;
      continue;
    }
    switch (arg[1]) {
    case 'd': 
      subdirperobject = true;
      break;
    case 'v': 
      WriteVignets = true;
      break;
    default : 
      cerr << "unknown option " << arg << endl;
      usage(argv[0]);
      break;
    }
  }

  ifstream lightfile(lightfilename.c_str());
  if (!lightfile) return EXIT_FAILURE;

  LightCurveList fids(lightfile);
  SimFitPhot doFit(fids);
  doFit.bOutputDirectoryFromName = subdirperobject;
  doFit.bWriteVignets = WriteVignets;

  for_each(fids.begin(), fids.end(), doFit);
  fids.write("lightcurvelist.dat");

  return EXIT_SUCCESS;
}

