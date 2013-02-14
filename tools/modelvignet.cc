#include <iostream>
#include <fstream>
#include <string>

#include <poloka/fileutils.h>
#include <poloka/fitsimage.h>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <galaxy_fits> <psf_fits> <size (even)> <flux> <result> " << endl ;
}

int main(int argc, char **argv)
{
  if (argc != 6)  {usage(argv[0]);  exit(1);}
  
  string filename;
  filename = argv[1]; if(! FileExists(filename)) { cerr << "FATAL " << filename << " does not exists" << endl; return EXIT_FAILURE;}
  FitsImage galaxy(filename);
  filename = argv[2]; if(! FileExists(filename)) { cerr << "FATAL " << filename << " does not exists" << endl; return EXIT_FAILURE;}
  FitsImage psf(filename);
  filename = argv[2]; if(! FileExists(filename)) { cerr << "FATAL " << filename << " does not exists" << endl; return EXIT_FAILURE;}
  
  
  unsigned int nx = atoi(argv[4]);
  if ((nx%2)==0 ) {
    cout << "size must be even" << endl;
    return EXIT_FAILURE;
  }

  

  int hx = (nx-1)/2;  
  double flux = atof(argv[5]);
  string output = argv[6];
  
  Image result(2*nx+1,nx);
  int hx_gal = (galaxy.Nx()-1)/2;
  int margin_to_remove = 4;

  int hx_psf = (psf.Nx()-1)/2;
  
  double val;
  double maxval=0;

  // galaxy
  for(int j=0;j<nx;j++) {
    int jgal = j-hx;
    if(abs(jgal)>hx_gal-margin_to_remove) continue;
    for(int i=0;i<nx;i++) {
      int igal = i-hx;
      if(abs(igal)>hx_gal-margin_to_remove) continue;
      val = galaxy(igal+hx_gal,jgal+hx_gal);
      if(val>0) {
	result(i,j)=val; // left
	result(i+nx+1,j)=val; // right
	if(maxval<val) maxval=val;
      }
    }
  }
  // add sn
  for(int j=0;j<nx;j++) {
    int jpsf = j-hx;
    if(abs(jpsf)>hx_psf) continue;
    for(int i=0;i<nx;i++) {
      int ipsf = i-hx;
      if(abs(ipsf)>hx_psf) continue;
      val = psf(ipsf+hx_psf,jpsf+hx_psf);
      if(val>0) {
	val*=flux;
	result(i+nx+1,j)+=val; // right
	if(maxval<result(i+nx+1,j)) maxval=result(i+nx+1,j);
      }
    }
  }
  // add a line
  for(int j=0;j<nx;j++) {
    result(nx,j)=maxval;
  }
  
  // write output
  {FitsImage(output,result);}
  
  return EXIT_SUCCESS;
}

