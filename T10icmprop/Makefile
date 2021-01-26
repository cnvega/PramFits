#######
# Makefile for icmpropertyfinder
#######

####### Select target computer
SYSTYPE = 'FCAGLP-horus'
#SYSTYPE = 'generic'
#SYSTYPE = 'FCAGLP-seminare'
#SYSTYPE = 'i686'

####### Select input format
# Uncomment IOFORMAT = 'HDF5' if the galaxies input files are in HDF5 
# format (either orbit or SAG special dump files). ICM data will also be 
# stored in HDF5 format. 
# Comment it out to use as input the special dump files in binary format;
# the output will also be in binary. 
# For some systems (horus) the location of the HDF5 library is known and is
# specified below, when setting the optimizations. For other systems, 
# specify here the path to the HDF5 library
IOFORMAT = 'HDF5'
#HDF5PATH = /usr/local/hdf5-1.8.8-linux-static

####### Use galaxy orbit files
# Read input from orbit files created by suborbitfinder, instead of using
# the special dump from SAG
ORBITS = 'ON'

####### Files
EXEC   = icmpropertyfinder

OBJS   = icmpropertyfinder.o allvars.o allocate.o allocategas.o \
         readparameterfile.o ngb.o ngbtree.o init.o io_input_G2.o \
         reordergas.o age.o rhoNFW.o rhoSIS.o \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/sort2_flt_int.o \
         nrsrc/trapzd.o nrsrc/polint.o
ifeq ($(IOFORMAT),'HDF5')
  OBJS += read_galaxies_HDF5.o dumpHDF5_icm.o indexrm.o
endif
ifeq ($(ORBITS),'ON')
  OBJS += io_orbits.o
endif

INCL   = proto.h nrsrc/nrutil.h allvars.h readparameterfile.h

.KEEP_STATE:

####### Compiler settings
CC       = gcc
CFLAGS   = -O2 -mtune=native
LIBS     = -lm
OPT      =

ifeq ($(ORBITS),'ON')
  OPT += -DORBITS
endif
ifeq ($(IOFORMAT),'HDF5')
  OPT += -DHDF5OUTPUT
endif

ifeq ($(SYSTYPE),'i686')
  CFLAGS += $(OPT) -Wno-unused-result
  ifeq ($(IOFORMAT),'HDF5')
    CC       = $(HDF5PATH)/bin/h5cc
    CFLAGS   += -I$(HDF5PATH)/include
    LIBS     += -L$(HDF5PATH)/lib 
  endif	
endif

ifeq ($(SYSTYPE),'FCAGLP-horus')
  #CFLAGS += $(OPT) -Wno-unused-result
  CFLAGS += $(OPT) -march=native -Wno-unused-result -w -std=c99 
  ifeq ($(IOFORMAT),'HDF5')
#    HDF5PATH = /usr/local/hdf5-1.8.9-linux-x86_64-static
    HDF5PATH = /usr/
    CC       = $(HDF5PATH)/bin/h5cc
    CFLAGS   += -I$(HDF5PATH)/include
    LIBS     += -L$(HDF5PATH)/lib 
  endif	
endif

ifeq ($(SYSTYPE),'FCAGLP-seminare')
  CFLAGS   += -Wall -pedantic -std=c99 -D_BSD_SOURCE \
              -I/opt/apps/gsl/gsl-1.15/include -Wno-unused-variable
  LIBS     += -L/opt/apps/gsl/gsl-1.15/lib
  ifeq ($(IOFORMAT),'HDF5')
    CC       = $(HDF5PATH)/bin/h5cc
    CFLAGS   += -I$(HDF5PATH)/include
    LIBS     += -L$(HDF5PATH)/lib 
  endif	
endif

ifeq ($(SYSTYPE),'generic')
  ifeq ($(OUTPUT),'HDF5')
    CC       = $(HDF5PATH)/bin/h5cc
    CFLAGS  += -I$(HDF5PATH)/include
    LIBS    += -L$(HDF5PATH)/lib
  endif
  CFLAGS += -march=native -Wno-unused-result -w -std=c99
endif


$(EXEC): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS)  -o $(EXEC) $(LIBS)  

$(OBJS): $(INCL) Makefile

clean:
	rm -f $(OBJS) $(EXEC)

tidy:
	rm -f $(OBJS)
