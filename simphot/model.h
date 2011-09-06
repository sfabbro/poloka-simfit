#ifndef MODEL__H
#define MODEL__H

#include <string>


/* the class that provides basic characteristics of the model
   that is going to be fitted. 
   - model size
   - model PSF (or mode precisely, the kernel to match a given image)
   - coordinate system

   it could be a virtual class with actual implementations, because 
   there are many ways to implement it. For the moment we stick to 
   an actual implementation.


*/

#include "gtransfo.h"
#include "reducedimage.h"
#include "intpoint.h"
#include "pixelblock.h"
#include "imagepsf.h"
#include "sestar.h"

class Mat;
class Vect;

/*! this class contains the subpart of SimPhotFit that deals with the galaxy model
  and the reference PSF. It is inherited by SimPhotFit. */

class Model
{
 private:
  const ReducedImageRef refImage;
  Point objectPosInImage; // in image coordinates, initial value
  IntPoint refPix; // integer offset
  const double overSampling;
  ImagePSF refPSF;
  SEStarList seRef;

 protected :
  PixelBlock galaxyPixels;
  Point objectPos; // in reduced coordinates (objectPosInImage-refPix)
  bool hasGalaxy;

 public :
  Model(const ReducedImage *RefImage, const Point &ObjectPos,
	const double OverSampling=1);

  //! position in reduced coordinates (eventually fitted)
  const Point &ObjectPos() const { return objectPos;}


  //! position in non-reduced coordinates (eventually fitted)
  const Point ObjectPosInImage() const { return refPix+objectPos;}

  const ReducedImage *RefImage() const { return refImage;}

  //!
  bool HasGalaxy() const { return hasGalaxy;}

  //  int HSizeX() const { return hSizeX;}
  // int HSizeY() const { return hSizeY;}
  bool Solve(Mat &A, Vect &B, const string &U_or_L, 
	     const bool FittingGalaxy, const int HalfKernelSize) const;

  const PixelBlock& GalaxyPixels() const { return galaxyPixels;}

  //! PSF centered on current object position, in ref coordinates
  void RefPSFPixels(PixelBlock &PSFPixels) const;

  bool FindTransfos(const ReducedImage* Current, 
		    GtransfoRef &Transfo2Ref,
		    GtransfoRef &TransfoFromRef, 
		    Point &ObjectPosInCurrent,
		    double &PhotomRatio);


  bool FindKernel(const ImagePSF &CurrentPSF,
		  const IntPoint &CurrentIntOffset,
		  const Gtransfo *Vignette2Model,
		  const Gtransfo *Model2Vignette,
		  PixelBlock &Kernel);

};


void ComputePSFPixels(const ImagePSF &PSF, 
		      const Point &PosInImage,
		      const IntPoint &IntOffset, 
		      PixelBlock &PSFPixels,
		      PixelBlock *PSFPixelsXDer=NULL,
		      PixelBlock *PSFPixelsYDer=NULL);



#endif /* MODEL__H */
