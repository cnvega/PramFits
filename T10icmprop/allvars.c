/**
 * @file allvars.c
 * @brief Global variables for use with @a icmpropertyfinder.c
 * @author Tomas E. Tecce
 */
#include "allvars.h"


char PreTag[100];
char OrdTag[100];

/* System of units  */
double UnitTime_in_s,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears,
  GravityConstantInternal,
  Hubble;

double  Time;
int     NumPart;
int     NumPart0;
double  PartMass;
double  PartMassType;

double  GPMass;

struct particle_data *P;

struct NODE *Nodes, *Nodes_base;

float *R2list;
int   *IndexList;

struct gaspart_data *PG;

int *Id;

/* B is the earlier point of time, i.e. z(B) > z(A) */
struct GALAXY_ICM *GalA, *GalB;  /* Galaxy population */

int *IPG_RP;   /* This table will hold all neighbour indices for RP*/

int TotNumGalB;
int TotNumGalA;

double RhoCrit;

float ZZ[OUTPUTS], Age[OUTPUTS];

float Zcurr;

double NeighbourSearchRadius;

char gaspart_fname[FILENAME_MAX]; /*for gas TOMAS */

int Sx, Sy, Sz, Sc;

float Hubble_of_z;

