#include <list>
#include <map>

#include <vignetserver.h>
#include <fitsimage.h>


#define SERVERDEBUG true
#define SUPERSHORT unsigned char

using namespace std;

// for storage purpose only
class KernelStorage : public Kernel {

private :
  
  float *float_data;
  SUPERSHORT *unsigned_data;
  

public :
  KernelStorage() :
    Kernel(),float_data(0),unsigned_data(0) {}
  
  KernelStorage(const Kernel& kern) :
    Kernel(kern),float_data(0),unsigned_data(0) {}
  
  void FreeMem() {
    delete [] data; // free mem 
    data   = 0;
    data00 = 0;
  }

  void SaveAsFloat() {
    // allocate float mem.
    int size = nx*ny;
    float_data = new float[size];
    
    // copy
    const DPixel *double_it = data;
    float *float_it         = float_data;
    for(int i=0; i<size;i++,double_it++,float_it++)
      *float_it = float(*double_it);

    // free int. Kernel mem
    FreeMem();
  }

  void SaveAsSUPERSHORT() {
    // allocate float mem.
    int size = nx*ny;
    unsigned_data = new SUPERSHORT[size];
    
    // copy
    const DPixel *double_it = data;
    SUPERSHORT *unsigned_it   = unsigned_data;
    for(int i=0; i<size;i++,double_it++,unsigned_it++)
      *unsigned_it = (SUPERSHORT)(*double_it);

    // free int. Kernel mem
    FreeMem();
  }

  
  void RestoreKernel() {
    
    if(data) return;

    if(float_data) {
      
      // allocate
      int size = Nx()*Ny();
      data = new DPixel[size];
      minindex=0; maxindex = max(nx*ny-1,0);
      
      // copy data
      DPixel *double_it     = data;
      const float *float_it = float_data;
      
      for(int i=0; i<size;i++,double_it++,float_it++)
	*double_it = double(*float_it);
      
      // need to refix some pointers
      data00 = &(*this)(hSizeX,hSizeY);
      minindex = begin()-data00; 
      maxindex = minindex + Nx()*Ny()-1;
    }
    if(unsigned_data) {
      
      // allocate
      int size = Nx()*Ny();
      data = new DPixel[size];
      minindex=0; maxindex = max(nx*ny-1,0);
      
      // copy data
      DPixel *double_it     = data;
      const SUPERSHORT *unsigned_it = unsigned_data;
      
      for(int i=0; i<size;i++,double_it++,unsigned_it++)
	*double_it = double(*unsigned_it);
      
      // need to refix some pointers
      data00 = &(*this)(hSizeX,hSizeY);
      minindex = begin()-data00; 
      maxindex = minindex + Nx()*Ny()-1;
    }      
  }
  
};



struct VignetData {
  
  Window window;
  KernelStorage kernel_storage;
  VignetData(const Window& window_i) 
    : window(window_i) {};
};

struct VignetDataForImage : public std::list<VignetData> {
  bool read;
  VignetDataForImage() : read(false) {};
};

static std::map<string, VignetDataForImage >* toto = 0; // image filename as key


static std::map<string, VignetDataForImage >& VignetServer() {
  if(toto == 0) toto = new std::map<string, VignetDataForImage >;
  return *toto;
}

void reserve_vignet_in_server(const std::string& fitsfilename, const Window& window) {
  if(SERVERDEBUG) cout << "SERVERDEBUG:  reserve " << fitsfilename << " for vignet" << endl;
  VignetServer()[fitsfilename].push_back(VignetData(window));
}

#define VALUE_WHEN_OUTSIDE -1.1e60

#include "fileutils.h"

void read_reserved_vignets(const string& fitsfilename) {
  
  if(SERVERDEBUG) cout << "SERVERDEBUG:  read all vignets from " << fitsfilename << endl;
  
  VignetDataForImage& vignets = VignetServer()[fitsfilename];
  string internalFileName = "/tmp/"+BaseName(DirName(fitsfilename))+"_"+CutExtension(BaseName(fitsfilename))+".fits"; 
  //   toto(mem://) is documented but does not works in real life;

  
  ImageCopy(fitsfilename, internalFileName);

  if(SERVERDEBUG) cout << "SERVERDEBUG:  copy " << fitsfilename << " in " << internalFileName << endl;
 
  for(VignetDataForImage::iterator v = vignets.begin(); v!= vignets.end(); ++v) {
    if (v->kernel_storage.Nx() == 0) {
      Kernel kern; 
      kern.readFromImage(internalFileName,v->window,VALUE_WHEN_OUTSIDE);
      v->kernel_storage = kern;
      
      // save data in different formats
      if(fitsfilename.find("satur") != string::npos) {
	//if(SERVERDEBUG) cout << "SERVERDEBUG:  save " << internalFileName << " as SUPERSHORT" << endl;
	v->kernel_storage.SaveAsSUPERSHORT();
      }else{
	//if(SERVERDEBUG) cout << "SERVERDEBUG:  save " << internalFileName << " as float" << endl;
	v->kernel_storage.SaveAsFloat();
      }
      
      // v->kern.readFromImage(internalFileName,v->window,VALUE_WHEN_OUTSIDE);
    }
  }
  unlink(internalFileName.c_str());
  vignets.read = true;
}

void get_vignet_from_server(const std::string& fitsfilename, const Window& window, Kernel& kern, double value_when_outside_fits) {
  
  
  // dimage.readFromImage(fitsfilename,window,value_when_outside_fits);


  VignetDataForImage& vignets    = VignetServer()[fitsfilename];
  VignetDataForImage::iterator v = vignets.begin();
  for(; v!= vignets.end(); ++v) {
    if(v->window == window)
      break;
  }

  if(v == vignets.end() ) {
    kern.readFromImage(fitsfilename,window,value_when_outside_fits);
    cout << "WARNING get_vignet_from_server no such window in server for file " << fitsfilename << endl;
    // if it happens to be requested once again, we'll have it...
    vignets.push_back(VignetData(window));
    vignets.back().kernel_storage = kern;
    return;
    //
    //abort();
  }
  
  if(! vignets.read)
    read_reserved_vignets(fitsfilename);
  
  //if(SERVERDEBUG) cout << "SERVERDEBUG:  restore kernel for " << fitsfilename << endl;
	
  v->kernel_storage.RestoreKernel();
  kern = v->kernel_storage;
  v->kernel_storage.FreeMem();

  //if(SERVERDEBUG) cout << "SERVERDEBUG:  fix value_when_outside_fits " << endl;

  // need to take care of value_when_outside_fits
  DPixel* the_end = kern.end();
  for(DPixel *pixel = kern.begin(); pixel != the_end; ++pixel) {
    if(*pixel == VALUE_WHEN_OUTSIDE) *pixel = value_when_outside_fits;
  }
}

