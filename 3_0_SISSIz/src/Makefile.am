bin_PROGRAMS = SISSIz

SISSIz_SOURCES=SISSIz.c rnaz_utils.c expfit.c bionj_custom.c bionj_custom.h sissiz.h rnaz_utils.h expfit.h cmdline_sissiz.c cmdline_sissiz.h treeML.h treeML.c

LDADD = -L../librna -lRNA -lm -L../levmar -llevmar -L../sissiz -lSISSI -L../phyml -lPHYML
RNAz_LINK = $(CXX) -o $@



if DISABLE_OPENMP
AM_CPPFLAGS  = -I../librna -I../sissiz -I../levmar -I../phyml
else
AM_LDFLAGS   = -fopenmp
AM_CPPFLAGS  = -fopenmp -I../librna -I../sissiz -I../levmar -I../phyml 

endif

