/**
 * @file allvars.h
 * @brief Header file with macros and global variables for use with
 * file @a icmpropertyfinder.c
 * @author Tomas E. Tecce
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>

#include "nrsrc/nrutil.h"
#include "nrsrc/nricm.h"
#include "proto.h"

#ifdef HDF5OUTPUT
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#define OUTPUTS   93 

#define STEPS 25 /*50, SOFIA (8/1/15) lo cambie para ser consistentes con el archivo reducido de orbitas
                                      ya generado*/ 

#define MAXGROUPS    80000

#define MAXGALAXIES   80000

#define MAX_NGB 20000

/* Number of particles around each gx among which metals are distributed */
#define NMAX 64               

/* Number of particles around each gx allowed when fixing radius */
#define NMAXR 7000 

/* Number of particles around each gx allowed when fixing radius for ram 
   pressure*/
#define NMAXR_RP 7000

/* Number of particles around each gx allowed for DM part. */
#define NMAXDM 5000

#define TOTMAXNEIGHBOURS_RP_FIXEDPART (NMAX*MAXGALAXIES)
#define TOTMAXNEIGHBOURS_RP_FIXEDRADIUS (NMAXR_RP*MAXGALAXIES)

#define  PI 3.141592654

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7


extern char PreTag[100];
extern char OrdTag[100];

/* System of units  */
extern double UnitTime_in_s,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears,
  GravityConstantInternal,
  Hubble;

extern double Time;

extern int     NumPart;
extern int     NumPart0;
extern double  PartMass;
extern double  PartMassType;

extern double  GPMass;

extern struct particle_data 
{
  float  Mass;
  float  Pos[3];
  float  Vel[3];
  int    Nextnode;
} *P;

extern struct NODE
{
  union
  {
    int suns[8];		/* ancestor nodes */
    struct 
    {
      float len; 	        /* sidelength of treenode */
      float center[3];	        /* geometrical center */
      int   sibling;
      int   nextnode;
    } d;
  } u;
} *Nodes, *Nodes_base;

extern float *R2list;
extern int   *IndexList;

/* Properties of gas particles  ############################# */
extern struct gaspart_data
{
  float  U;
  float  Rho;
  float  Ne;
  float  Temp; 
} *PG;

extern int *Id;

/*
 * This version of struct GALAXY is heavily reduced (includes only
 * what is necessary to store ICM properties)
 */
struct GALAXY_ICM
{
  int PaIndex;
  int ParentSubGroup;
  int ParentGroup;
  int Id;
  int Type;
  int CentralGal;

  float Pos[3];
  float Vel[3];
  /*float PosComov[3], VelComov[3];*/
  float Rvir;
  float Mvir;
  float Rscale;

  float ParentGroupMvir;
  float ParentGroupRvir;
  float ParentGroupVvir;
  float ParentGroupCvir;

  float Rgroup;
  float RhoICM_SIS;
  float RhoICM_NFW;

#ifdef ORBITS
  int   OrbitType[STEPS];
  float OrbitX[STEPS];
  float OrbitY[STEPS];
  float OrbitZ[STEPS];
  float OrbitVx[STEPS];
  float OrbitVy[STEPS];
  float OrbitVz[STEPS];
#endif
  
  /* For mergers */
  //  float MergTime;
  //  float MergDist;
  //  float DistSat;
  //  float VcSat;
  
  /* To store the results */
  int   IPG_startRP;   /* First entry in neighbour table */
  int   IPG_numberRP;  /* Number of neighbours */
  float meanRhoICM;
  float medianRhoICM;
  float relvel2ICM;
  float HubbleType;
  float FGravHI;
#ifdef ORBITS
  float PeakRPMedianRhoICM;
  float PeakRPrelvel2ICM;
  float PosPeakRP[3];
  float VelPeakRP[3];
  float OrbitMedianRhoICM[STEPS];
  float OrbitVrel2ICM[STEPS];
#endif
} *GalA, *GalB;  /* Galaxy population */
/* B is the earlier point of time, i.e. z(B) > z(A) */

/* This table will hold all neighbour indices for RP */
extern int *IPG_RP;   

extern int TotNumGalB;
extern int TotNumGalA;

extern double RhoCrit;

extern float ZZ[OUTPUTS], Age[OUTPUTS];

extern float Zcurr;

extern double NeighbourSearchRadius;

extern char gaspart_fname[FILENAME_MAX];

extern int Sx, Sy, Sz, Sc;

extern float Hubble_of_z;
