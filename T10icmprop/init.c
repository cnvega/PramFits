/*
  init.c
  
Version for use with icmpropertyfinder.c

*/
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void init(void)
{
  double a;
  double aext[MaxSnapshot];
  int snapshot;
  int totMaxNgb;  
  FILE  *fd; 
  char  buf[111];
  long long bytes;
  
  
  if (OutputListOn == 1) { /* Reading outputs time from external file 
			      (Dolag Sim) */     
    sprintf(buf, "%s/%s", path2, NameOutputsSelection);
    
    printf("Reading outputs from file '%s'\n", buf);
    
    if (!(fd = fopen(buf, "r"))) {
      printf("Error (init): cannot read file '%s': %s (%u)\n", 
	     buf, strerror(errno), errno); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    
    for (snapshot = 0; snapshot <= MaxSnapshot; snapshot++) {
      fscanf(fd," %lf ",&aext[snapshot]);
      ZZ[snapshot]  =  1.0/aext[snapshot] - 1;
      Age[snapshot] = time_to_present(ZZ[snapshot]);
    }
  }
  else {
    for (snapshot = 0, a = TimeOfFirstSnapshot; snapshot <= MaxSnapshot; 
	 snapshot++, a *= TimeBetSnapshot) {
      ZZ[snapshot]  =  1.0/a - 1;
      Age[snapshot] = time_to_present(ZZ[snapshot]);
    }
  }
  
  /* Memory allocation for galaxies */
  if (!(GalA = malloc(bytes=(sizeof(struct GALAXY_ICM)*(MAXGALAXIES+1))))) {
    printf("Error (init): cannot allocate memory for galaxies: %s (%u)\n",
	   strerror(errno), errno);
    printf("bytes=%d%09d sizeof(struct GALAXY)=%d MAXGALAXIES=%d %d \n",
	   (int) (bytes / 1000000000),
	   (int) (bytes % 1000000000),
	   (int)sizeof(struct GALAXY_ICM), MAXGALAXIES, (int)sizeof(size_t));
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  
  /* Memory allocation for gas */
  totMaxNgb = TOTMAXNEIGHBOURS_RP_FIXEDRADIUS;
  if (!(IPG_RP = malloc(sizeof(int)*totMaxNgb))) {
    printf("Error (init): cannot allocate memory for IPG_RP: %s (%u)\n", 
	   strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }

  return;
}


void set_units(void)   /* Set some units */
{
  
  /*  
   * Standard choice of units:
   * UnitLength_in_cm         = 3.085678e21;
   * UnitMass_in_g            = 1.989e43; 
   * UnitVelocity_in_cm_per_s = 1e5;
   */
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s/SEC_PER_MEGAYEAR;
   
  printf("Unit_time_in_Megayears:%e\n", UnitTime_in_Megayears);

  G = GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
  
  UnitDensity_in_cgs = UnitMass_in_g/pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs = UnitMass_in_g/UnitLength_in_cm/pow(UnitTime_in_s,2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs/UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
  
  /* Convert some physical input parameters to internal units */
  Hubble = HUBBLE * UnitTime_in_s;

  /* Compute a few quantitites */
  RhoCrit = 3*Hubble*Hubble/(8*PI*G);

  return;
}
