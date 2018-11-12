===========================  OVERVIEW  =================================

Get the latest version at

   src    -

   app    -

   test   -

The current release includes all the above packages plus a few 

supplementary packages:

   lapack - a collection of a few linear algebra algorithms (Real*8). 

            We recommend to use the corresponding system library if you

            have one.

   blas   - a collection of basic linear algebra Real*8 subroutines such

            as x*y, a*x+y, etc where x and y are vectors. We recommend 

            to use the corresponding system library if you have one.

=======================  INSTALLATION DETAILS  =========================

We uses a simplified 'Makefile' to compile and install the package. 

User should set up compiler names (FORTRAN and C) and compilation options 

in the file config/make.inc. Type 'make help' for more details.

$ vi config/make.inc

provide correct names and options

$ make libs

$ make install

=========================  RUNNING DEMOs  ==============================

