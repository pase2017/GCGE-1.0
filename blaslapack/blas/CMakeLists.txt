#
#  KONA
#  Basic linear algebra package BLAS-3.2
#

add_library(blas-3.2 STATIC
            "dsymm.f"
            "daxpy.f"  "dcopy.f"  "ddot.f"   "dnrm2.f"  
            "dgemm.f"  "dgemv.f"  "dger.f"   "dtrsv.f"
            "dtrsm.f"  "dscal.f"  "dsyr.f"            
            "dswap.f"  "dsymv.f"  "dsyr2.f"  "dsyr2k.f"  "dsyrk.f" 
            "drot.f"   "drotg.f"              
            "dtrmm.f"  "dtrmv.f"  "idamax.f" "lsame.f"
            "xerbla.f")

install(TARGETS blas-3.2 DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




