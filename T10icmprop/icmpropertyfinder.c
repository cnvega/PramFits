/**
 * @file icmpropertyfinder.c
 * @brief Program to find ICM properties such as density, velocities of 
 * galaxies relative to ICM gas (to determine the ram pressure exerted on 
 * the galaxy).
 * @author Tomas E. Tecce
 *
 * This program reads input from the files created running SAG with the
 * option SpecialDumpOn = 1, or if compilation option ORBITS is enabled
 * then the input is read from the galaxy orbit files created by the
 * suborbitfinder program.

 * Use option ComovingToPhysicalOn = 1 if the positions and velocities in 
 * the input files are given in comoving coordinates. All other data (e.g. 
 * virial radii) are assumed to be in physical coordinates. 

 * Values of ICM density and relative velocity are returned in physical 
 * (code) units. If ORBITS is enabled, the ICM density and galaxy relative
 * velocity is obtained at each point stored along the galaxy's orbit.
 *
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


int main(int argc, char **argv)
{
  FILE *fp;
  char parameterfile[FILENAME_MAX];
  char input_fname0[FILENAME_MAX], input_fname[FILENAME_MAX];
  char catalogue_fname[FILENAME_MAX], dump_fname[FILENAME_MAX];
  char orderedGP_fname[FILENAME_MAX];
#ifdef ORBITS
  char suborbits_fname[FILENAME_MAX];
#endif
  char hostname[256];
  int i, p, snapshot0, snapshot, type;
  double partMassGas;
  time_t currtime, t0, t1;
  clock_t c0, c1;
  struct tm *loctime;


  t0 = time(NULL);
  c0 = clock();
  printf("#######\n"
	 "# icmpropertyfinder\n"
	 "#######\n\n");
  gethostname(hostname, 256);
  currtime = time(NULL);
  loctime = localtime(&currtime);
  printf("Running on %s - Run started %s\n", hostname, asctime(loctime)); 
  fflush(stdout);

  
  if(argc != 2 && argc != 3) {
    fprintf(stderr, "Error: no parameter file given\n");
    exit(EXIT_FAILURE);
  }
  
  if (argc >= 3) strcpy(PreTag, "itf_");
  else PreTag[0] = 0;
  
  sprintf(parameterfile, "%s", *(argv+1));
  readparameterfile(parameterfile); 

  printf("\nModel Parameters\n"
         "----------------\n"
	 "Seed                      %d\n"
	 "Path1                     %s\n"
         "Path                      %s\n"
         "Path2                     %s\n"
         "Path3                     %s\n"
         "FileBasename              %s\n"
         "Files                     %d\n"
	 "StartSnapshot             %d\n"
         "MaxSnapshot               %d\n"
         "Omega                     %g\n"
         "OmegaLambda               %g\n"
         "Hubble_h                  %g\n"
	 "H0                        %g\n"
	 "G                         %g\n"
	 "BaryonFrac                %g\n"
	 "GasSearchRadiusGx0        %g\n"
	 "GasSearchRadius           %g\n"
	 "GasConstPartNum           %d\n"
	 "GasMinimum                %d\n"
	 "Niterations               %d\n"
	 "Fracmax                   %g\n"
         "Identifier                %s\n"
         "TimeBetSnapshot           %g\n"
         "TimeOfFirstSnapshot       %g\n"
         "UnitLength_in_cm          %g\n"                  
         "UnitMass_in_g             %g\n"
         "UnitVelocity_in_cm_per_s  %g\n",        
         Seed, path1, path, path2, path3, filebasename, Files, 
	 StartSnapshot, MaxSnapshot, Omega, OmegaLambda, 
	 Hubble_h, H0, G, BaryonFrac, GasSearchRadiusGx0, 
	 GasSearchRadius, GasConstPartNum, GasMinimum, 
	 Niterations, Fracmax, Identifier, TimeBetSnapshot, 
	 TimeOfFirstSnapshot, UnitLength_in_cm, UnitMass_in_g, 
	 UnitVelocity_in_cm_per_s);
  
  printf("\nCode Settings\n"
         "-------------\n"
	 "OutputListOn              %d\n"
	 "ReorderOn                 %d\n\n",
	 OutputListOn, ReorderOn);
  fflush(stdout);
    
  snapshot0 = 0;
  sprintf(input_fname0, "%s%s%s%03d", path1, NameSnapshot, 
	  filebasename, snapshot0);
  NumPart0 = getNumPart(input_fname0, Files, type=0);
  printf("NumPart0 = %d\n", NumPart0);

  /*
   * Keeps the number of gas particles at snapshot = 0
   * necessary for simulations with SF where the number of gas
   * particles is reduced during evolution
   */
  loadpositions(input_fname0, Files, type=0);

  partMassGas = PartMassType; /* Mass of gas particles. Taken from
				 snapshot = 0 because it defines 
				 the baryon fraction for the 
				 simulation. In later snapshots 
				 the gas particle mass diminishes 
				 if there is SF */
  printf("NumPartGas = %d, PartMassGas = %f\n\n", NumPart, partMassGas);
  fflush(stdout);
  
  init();
  set_units();
  
  
  /**************************************
   * Begin loop over selected snapshots *
   **************************************/
  for (snapshot = StartSnapshot; snapshot <= MaxSnapshot; snapshot++) {
 // for (snapshot = 85; snapshot <= MaxSnapshot; snapshot++) {
    
    /* Determine Hubble constant for this redshift */
    Zcurr = ZZ[snapshot];
    Hubble_of_z = H0 * sqrt(Omega*pow(1+Zcurr, 3) +
			    (1-Omega-OmegaLambda)*(1+Zcurr)*(1+Zcurr) +
			    OmegaLambda);
   
    printf("\n-------\n");
    printf("Snapshot %d, Redshift %g, H(z) = %g\n\n", snapshot, Zcurr, 
           Hubble_of_z); fflush(stdout);

#ifndef ORBITS
    /* Load galaxies data (read from SpecialDump files) */
#ifdef HDF5OUTPUT
    sprintf(catalogue_fname, "%s%s/Subgalaxies/specialDump_%s%s_%03d.hdf5",
	    path3, Identifier, PreTag, Identifier, snapshot);
    get_galaxies_data_HDF5(snapshot, catalogue_fname);
#else
    sprintf(catalogue_fname, "%s%s/Subgalaxies/specialDump_%s%s_%03d.dat",
	    path3, Identifier, PreTag, Identifier, snapshot);
    get_galaxies_data(snapshot, catalogue_fname);
#endif
#else /*ORBITS*/
    /* Orbits data files are only in HDF5 format */
#ifdef HDF5OUTPUT
    sprintf(suborbits_fname, "%s%s_%s%03d.hdf5", path, NameSuborbits,
	    PreTag, snapshot);
    read_suborbits_A(suborbits_fname,snapshot);
#else
    fprintf(stderr, "Error (main): to read orbits data files, the "
	    "code must be compiled with HDF5 support - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
#endif
#endif /*ORBITS*/
    
    if (TotNumGalA > 0) {
      if (ReorderOn == 1) { /* Gas particles are reordered here */
	sprintf(input_fname, "%s%s%s%03d", path1, NameSnapshot, filebasename,
		snapshot);
	loadpositions(input_fname, Files, type=0);
	
	printf("\nNumPartGas = %d, NumPart0 = %d\n", NumPart, NumPart0);
	
	reorderinggas();
      }
      else {                     /* Ordered gas particles read from file */
	allocate_memory();         /* Generate memory for Pos, Vel, Mass */
	strcpy(OrdTag, "ord");
	sprintf(orderedGP_fname, "%s%s%s%s_%03d", path1, 
		NameSPHordered, filebasename, OrdTag, snapshot);
	//sprintf(orderedGP_fname, "%s%s%s%03d", path1, 
	//	NameSPHordered, filebasename, snapshot);
	
	if (!(fp = fopen(orderedGP_fname, "r"))) {
	  printf("Error: cannot open file '%s' for input of "
		 "ordered gas particles: %s (%u)\n", 
		 orderedGP_fname, strerror(errno), errno); fflush(stdout);
	  exit(EXIT_FAILURE);
	}
	printf("\nReading from file with %d ordered gas "
	       "particles: %s\n", NumPart0, orderedGP_fname);
	
	for (p=1; p<=NumPart0; p++)
	  fread(&P[p].Pos[0], sizeof(float), 3, fp);
	for (p=1; p<=NumPart0; p++)
	  fread(&P[p].Vel[0], sizeof(float), 3, fp);
	for (p=1; p<=NumPart0; p++)
	  P[p].Mass = partMassGas;
	
	allocate_gasmemory();
	
	for (p=1; p<=NumPart0; p++)
	  fread(&PG[p].U, sizeof(float), 1, fp);
	for (p=1; p<=NumPart0; p++)
	  fread(&PG[p].Rho, sizeof(float), 1, fp);
	
	/* Free unused memory */
	printf("Freeing space for particle ID\n"); 
	Id++;
	free(Id);
      }

      /* Convert snapshot particle data from comoving to physical
	 coordinates: positions, velocities and particle densities */
      /* 
       * To convert Gadget comoving quantities to physical:
       * Positions: divide by 1+z
       * Velocities: divide by sqrt(1+z)
       * Densities: multiply by (1+z)^3
       */
      /* TOMAS 2013-02-12: not anymore */
      /*for (p=1; p<=NumPart0; p++) {
	for (i=0; i<3; i++) {
	  P[p].Pos[i] /= 1.+Zcurr;
	  P[p].Vel[i] /= sqrt(1. + Zcurr);*/   /* To peculiar velocity */
	  //P[p].Vel[i] += Hubble_of_z*P[p].Pos[i]; /* Add Hubble flow */
      /*}
	PG[p].Rho *= pow(1.+Zcurr, 3);
      }*/
      
      /*
       * Defining a grid and finding the cells in which each 
       * galaxy and gas particle reside in order to get the 
       * index of the NMAX nearest gas particles around each 
       * galaxy or the nearest particles that reside inside a 
       * sphere of radius NeighbourSearchRadius around each 
       * galaxy (max allowed value: NMAXR), depending on whether 
       * the option FixedRadiusOn is off or on, respectively.
       */
      find_gaspart_withtree(snapshot);
    }
    
    /* Save data to file */
#ifdef HDF5OUTPUT
    sprintf(dump_fname, "%s%s_%s%03d.hdf5", path, NameICMproperties, 
	    PreTag, snapshot);
    dumpHDF5_icm(dump_fname,snapshot,Zcurr); 
#else
    sprintf(dump_fname, "%s%s_%s%03d", path, NameICMproperties, PreTag, 
	    snapshot);
    dump_icm(dump_fname,snapshot);  
#endif
    
    if (TotNumGalA > 0) {
      free_memory();
      free_gasmemory();
    }
  }  
  
  printf("Program finished - reached snapshot %d\n", snapshot);
  fflush(stdout);

  t1 = time(NULL);
  c1 = clock();
  printf("Elapsed wall clock time: %ld\n", (long)(t1-t0));
  printf("Elapsed CPU time:        %f\n", (float)(c1-c0)/CLOCKS_PER_SEC);
  fflush(stdout);
  
  return EXIT_SUCCESS;
}
  

/**
 * @brief Read input data from the special dump files generated by running 
 * @a subgalaxyfinder with the option SpecialDumpOn = 1 selected 
 *
 * Use this function for the original binary format; for the new files in 
 * HDF5 format, see function @a get_galaxies_data_HDF5
 */
void get_galaxies_data(int snapshot, char *catalogue_fname)
{
  FILE *fp;
  int p, i;
  float a, e, t;
  
  printf("Reading galaxies data from file '%s'...\n", catalogue_fname);
  if (!(fp = fopen(catalogue_fname, "r"))) {
    printf("Error (get_galaxies_data): cannot open file '%s': %s (%u)\n", 
	   catalogue_fname, strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }

  fread(&TotNumGalA, sizeof(int), 1, fp);
  printf("Total galaxies: %d\n", TotNumGalA);

  fread(&a, sizeof(int), 1, fp);
  fread(&e, sizeof(int), 1, fp);
  fread(&t, sizeof(int), 1, fp);
  printf("Alpha = %f, Epsilon = %f, ThreshMerger = %f\n", a, e, t);
  
  fread(&Seed, sizeof(int), 1, fp);
  fread(&Sx, sizeof(int), 1, fp);
  fread(&Sy, sizeof(int), 1, fp);
  fread(&Sz, sizeof(int), 1, fp);
  fread(&Sc, sizeof(int), 1, fp);
  printf("Seed = %d\n", Seed);  
  
  for (p = 1; p <= TotNumGalA; p++)
    fread(&GalA[p].ParentGroup, sizeof(int), 1, fp);
  
  for (p = 1; p <= TotNumGalA; p++) {
    fread(&GalA[p].Type, sizeof(int), 1, fp);
    if (GalA[p].Type > 3 || GalA[p].Type < 0) {
      printf("Read error (get_galaxies_data): "
	     "GalA[%d].Type = %d cannot be\n", p, GalA[p].Type);
      fflush(stdout);
      exit(EXIT_FAILURE);
    }
  }
  
  for (p = 1; p <= TotNumGalA; p++)
    fread(&GalA[p].Mvir, sizeof(float), 1, fp);
  
  for (p = 1; p <= TotNumGalA; p++)
    fread(&GalA[p].Rvir, sizeof(float), 1, fp);
  
  for (p = 1; p <= TotNumGalA; p++)
    fread(&GalA[p].Rscale, sizeof(float), 1, fp);
  
  for (p = 1; p <= TotNumGalA; p++) 
    fread(&GalA[p].Pos[0], sizeof(float), 3, fp);
  
  for (p = 1; p <= TotNumGalA; p++)
    fread(&GalA[p].Vel[0], sizeof(float), 3, fp);
  
  if (TotNumGalA > 0) {
    printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (comoving)\n",
  	   GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
    printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (comoving)\n",
	   GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);
    
    /* TOMAS 2013-02-12
       Only convert to physical coordinates at the end */

    /* Conversion of position and velocity to physical coordinates */
    /*for (p=1; p<=TotNumGalA; p++) {
      for (i=0; i<3; i++) {
	GalA[p].PosComov[i] = GalA[p].Pos[i];
	GalA[p].VelComov[i] = GalA[p].Vel[i];
	GalA[p].Pos[i] /= 1.+Zcurr;
	GalA[p].Vel[i] /= sqrt(1.+Zcurr);*/      /* To peculiar velocity */
	/* Add Hubble flow */
	/*GalA[p].Vel[i] += Hubble_of_z*GalA[p].Pos[i];*/
    /*}
      printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (physical)\n",
             GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
      printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (physical)\n",
             GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);
	     }*/
  }
  else 
    printf("No galaxies found in this snapshot\n");
  
  fclose(fp);
  return;
}


/**
 * @brief Write ICM data to output files
 * @param dump_fname Name of output file
 * @param snapshot Number of current snapshot
 *
 * Remember to reflect changes to this function in the function 
 * @a get_icm_properties in file @a environment_v4.c in SAG
 *
 * @warning All positions, velocities and densities are always
 * converted to physical units; only Pos and Vel are converted back to
 * comoving (to match the snapshots and SAG)
 */
void dump_icm(char *dump_fname, int snapshot)
{
  FILE *fp;
  int count_total, count_total_satellites;
  int count_cluster, count_cluster_satellites;
  int i, n, p, check;
  float buffer;


  count_total   = count_total_satellites   = 0;
  count_cluster = count_cluster_satellites = 0;
  
  for (i=1; i<=TotNumGalA; i++) {
    count_total++;
    if (GalA[i].Type == 2)
      count_total_satellites++;
    
    if (GalA[i].ParentGroup == 1) {
      count_cluster++;
      if( GalA[i].Type == 2)
	count_cluster_satellites++;
    }
  }

  printf("Galaxies: total = %d  satellites = %d\n",
	 count_total, count_total_satellites);
  printf("  in cluster: total = %d  satellites = %d\n\n",
	 count_cluster, count_cluster_satellites);
  
  if (!(fp = fopen(dump_fname, "w"))) {
    printf("Error(dump_icm): cannot open file '%s' for output: %s (%u)\n",
	   dump_fname, strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  
  n = count_total;
  check = fwrite(&n, sizeof(int), 1, fp);    checkwrite(check, 1);
  
  check = fwrite(&Seed, sizeof(int), 1, fp); checkwrite(check, 1);
  check = fwrite(&Sx, sizeof(int), 1, fp);   checkwrite(check, 1);
  check = fwrite(&Sy, sizeof(int), 1, fp);   checkwrite(check, 1);
  check = fwrite(&Sz, sizeof(int), 1, fp);   checkwrite(check, 1);
  check = fwrite(&Sc, sizeof(int), 1, fp);   checkwrite(check, 1);
  
  check = fwrite(&GasSearchRadius, sizeof(float), 1, fp);
  checkwrite(check, 1);
  
  check = fwrite(&GasSearchRadiusGx0, sizeof(float), 1, fp);
  checkwrite(check, 1);
  
  check = fwrite(&GasConstPartNum, sizeof(int), 1, fp);
  checkwrite(check, 1);
  
  check = fwrite(&Niterations, sizeof(int), 1, fp);
  checkwrite(check, 1);
  
  check = fwrite(&Fracmax, sizeof(float), 1, fp);
  checkwrite(check, 1);
  
  if (TotNumGalA != 0) {
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].ParentGroup, sizeof(int), 1, fp);
      checkwrite(check, 1);
    }
    
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].Type, sizeof(int), 1, fp);
      checkwrite(check, 1);
    }
    
    /* Physical position and velocity */
    /* TOMAS 2013-02-12
       Now we do not convert the positions and velocities to physical */
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].Pos[0], sizeof(float), 3, fp);
      checkwrite(check, 3);
    }
      
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].Vel[0], sizeof(float), 3, fp);
      checkwrite(check, 3);
    }
    
    /* Original comoving position */
    /*for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].PosComov[0], sizeof(float), 3, fp);
      checkwrite(check, 3);
      if (p == 1) 
	printf("GalA[1].PosComov[0] = %g\n", GalA[p].PosComov[0]);
	}*/

    /* Original comoving velocity */
    /*for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].VelComov[0], sizeof(float), 3, fp);
      checkwrite(check, 3);
      if (p == 1) 
	printf("GalA[1].VelComov[0] = %g\n", GalA[p].VelComov[0]);
	}*/

    /* Convert density to physical units */
    for (p=1; p<=TotNumGalA; p++) {
      buffer = GalA[p].meanRhoICM*pow(1+ZZ[snapshot],3);
      check = fwrite(&buffer, sizeof(float), 1, fp);
      checkwrite(check, 1);
    }
    
    for (p=1; p<=TotNumGalA; p++) {
      buffer = GalA[p].medianRhoICM*pow(1+ZZ[snapshot],3);
      check = fwrite(&buffer, sizeof(float), 1, fp);
      checkwrite(check, 1);
    }
    
    /* Convert relative velocity to physical units 
       NOTE: this is a squared velocity, so the conversion factor is
       squared as well */
    for (p=1; p<=TotNumGalA; p++) {
      buffer = GalA[p].relvel2ICM/(1+ZZ[snapshot]);
      check = fwrite(&buffer, sizeof(float), 1, fp);
      checkwrite(check, 1);
    }
    
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].IPG_numberRP, sizeof(int), 1, fp);
      checkwrite(check, 1);
    }

    /* Analytic ICM densities already are in physical units */
    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].RhoICM_SIS, sizeof(float), 1, fp);
      checkwrite(check, 1);
    }

    for (p=1; p<=TotNumGalA; p++) {
      check = fwrite(&GalA[p].RhoICM_NFW, sizeof(float), 1, fp);
      checkwrite(check, 1);
    }
  }
  fclose(fp);

  return;
}


/**
 * @brief Check data written to file
 */
void checkwrite(int check, int length)
{
  if (check != length) {
    printf("ERROR (checkerror) during writing. check=%d, length=%d. Stop\n",
	   check,length);
    exit(EXIT_FAILURE);
  }  
}
