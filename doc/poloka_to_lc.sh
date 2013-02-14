#!/bin/bash 

git rm autogen.sh ChangeLog config.h.in configure.ac Makefile.am README stamp-h.in 
git rm -r bindings
git rm -r cern_stuff cern_utils
git rm -r cmt
git rm -r dao
git rm -r flat
git rm -r imagemagick_utils
git rm -r m4
git rm -r mc
git rm -r simphot
git rm -r src
git rm -r src_base
git rm -r src_utils
git rm -r lapack_stuff
git rm -r telinst 
git rm -r psf 

mkdir tools
git add tools

# utils 
git mv utils/calibrate.cc utils/fitnight.cc utils/make_single_lightcurve.cc utils/calibratefixedpos.cc utils/make_lightcurve.cc utils/modelvignet.cc tools/
git rm -r utils

# lc
git rm lc/*dao*
git rm lc/*guru*
git rm lc/psfphot.cc
git rm lc/Makefile.am 
git mv lc src 
