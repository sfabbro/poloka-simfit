// This may look like C code, but it is really -*- C++ -*-
#ifndef LCUTILS__H
#define LCUTILS__H

#include <reducedimage.h>
#include "refstar.h"

//! read a list of objects and a ReducedImageList from a light curve file 
istream& lc_read(istream& Stream, RefStarList& Objects, ReducedImageList& Images);

//! create a light curve file from a list of objects and a ReducedImageList
ostream& lc_write(ostream& Stream, const RefStarList& Objects, const ReducedImageList& Images);

/*!
 \page lightfile Syntax of the "lightfile"
 
 A light curve file consists of objects and images.
 Comment are coded with the '#' at the beginning of the line.
 The expected format is the following:
 
 \code
 
 # X Y [DATE_MIN=MJD_MIN DATE_MAX=MJD_MAX NAME=toto BAND=xyz]
 # X Y : starting coordinates of the supernova
 # Sn flux will be fitted on images within [date_min, date_max], and forced to zero outside. 
 # NAME and BAND are used to name files on output. 
 OBJECTS
 510.20 1043.98 DATE_MIN=54xxx DATE_MAX=54yyy NAME=R4D14-4
 
 # The list of all DbImages in a single filter
 IMAGES
 DbImage1
 DbImage2
 DbImage3
 DbImage4
   .
   .
   .
 # photometric reference. It should be the best seeing image 
 PHOREF
 PhoRefDbImage
 
 \endcode
*/ 


#endif // LCUTILS__H
