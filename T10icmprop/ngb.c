/**
 * @file ngb.c
 * @brief Find neighbouring gas particles. Version for use with
 * @a icmpropertyfinder.c
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Find neigbouring gas particles using a tree search
 * @param snapshot Snapshot to process
 */
void find_gaspart_withtree(int snapshot)
{
  int i, k, p, numngb, totnumngbRP;
  int totMaxNgb;
  int istep, nstep, cont;
  float pram;
  
  ngb_treeallocate(0.7*NumPart0);
  force_treebuild();
  
/*  printf("\nFinding neighbours...\n");*/
  
  totnumngbRP = 0; 

  /* 
   * Find a constant number GasConstPartNum of neighbours 
   * (only if GasSearchRadius is set to 0 in the parameters file)
   */
  if (GasSearchRadius == 0.0) {
    
    if (GasConstPartNum == 0) {
      printf("Error (find_gaspart_withtree): GasSearchRadius and "
	     "GasConstPartNum cannot be both 0\n");
      printf("Please correct values in parameters file and "
	     "rerun - Exit\n"); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    
    totMaxNgb = TOTMAXNEIGHBOURS_RP_FIXEDPART;

    /******* Begin loop over all galaxies *******/
    for (p=1; p<=TotNumGalA; p++) {
/*
      printf("p = %d / %d, Type = %d ParentGroup = %d\n", p, TotNumGalA, 
	     GalA[p].Type, GalA[p].ParentGroup);
      printf("Pos = [%g %g %g] (comoving)\n", GalA[p].Pos[0], 
	     GalA[p].Pos[1], GalA[p].Pos[2]);
      if (GalA[p].Type > 0) 
	printf("Central galaxy = %d\n", GalA[p].CentralGal);
      printf("Radius within FOF group = %g (%g comoving)\n", 
	     GalA[p].Rgroup, GalA[p].Rgroup*(1+ZZ[snapshot]));
      printf("Virial radius of FOF group = %g (%g comoving)\n",
	     GalA[p].ParentGroupRvir, 
	     GalA[p].ParentGroupRvir*(1+ZZ[snapshot]));
      printf("Searching for %d first neighbours...\n", GasConstPartNum);
      fflush(stdout);
  */    
#ifdef ORBITS
      /* 
       * If ORBITS is enabled, we find the properties in all the
       * nstep outputs of the orbit (or until the galaxy is flagged as
       * merged) 
       * NOTE: type 0 and 1 galaxies have a single position stored, so they
       * have to be treated separately 
       */
      pram = 0.0;
      for (nstep=0; nstep<STEPS; nstep++) {
	/* The positions and velocities along the orbit should be
	   in physical coordinates */
	GalA[p].Pos[0] = GalA[p].OrbitX[nstep];
	GalA[p].Pos[1] = GalA[p].OrbitY[nstep];
	GalA[p].Pos[2] = GalA[p].OrbitZ[nstep];
	GalA[p].Vel[0] = GalA[p].OrbitVx[nstep];
	GalA[p].Vel[1] = GalA[p].OrbitVy[nstep];
	GalA[p].Vel[2] = GalA[p].OrbitVz[nstep];
	
	/* Proceed for all galaxies only in the last step;
	   in all previous steps, calculate only if the galaxy
	   is of type 2 (or 4) */
	cont = 0;
	if (nstep<STEPS-1) {
	  if (GalA[p].OrbitType[nstep]==2 || GalA[p].OrbitType[nstep]==4)
	    cont = 1;
	}
	else 
	  cont = 1;
	
	if (cont > 0) {
	  ngb_treefind(GalA[p].Pos, GasConstPartNum, 0);
	  
	  GalA[p].IPG_startRP = totnumngbRP;
	  
	  for (k=0; k<NMAX; k++) {
	    if (totnumngbRP >= totMaxNgb) {
	      printf("GAS SEARCH: maximum allowed total neighbour %d "
		     "reached for galaxy p=%d, nstep=%d\n", 
		     totMaxNgb, p, nstep); fflush(stdout);
	      exit(EXIT_FAILURE);
	    }
	    
	    IPG_RP[GalA[p].IPG_startRP + k]= IndexList[k];
	    /* TOMAS 2012-08-09: when calculating ICM properties along
	       a galaxy's orbit, this number becomes huge in the larger
	       simulations. I think it's because it only has to be
	       increased once in the STEPS */
	    if (nstep == STEPS-1) totnumngbRP++;
	  }
	  
	  GalA[p].IPG_numberRP = GasConstPartNum;
	  
	  if (nstep == 0 || nstep == STEPS-1) {
	    //printf("Computing ICM properties (nstep=%d)...\n", nstep); 
	    fflush(stdout);
	  }
	  
	  get_ICM_properties(p,nstep,snapshot);
	  
	  /* Find the peak RP along the orbit, and store the 
	     corresponding position and velocity */
	  if (GalA[p].OrbitMedianRhoICM[nstep]*
	      GalA[p].OrbitVrel2ICM[nstep] > pram) {
	    pram = GalA[p].OrbitMedianRhoICM[nstep]*
	      GalA[p].OrbitVrel2ICM[nstep];
	    GalA[p].PeakRPMedianRhoICM = GalA[p].OrbitMedianRhoICM[nstep];
	    GalA[p].PeakRPrelvel2ICM = GalA[p].OrbitVrel2ICM[nstep];
	    GalA[p].PosPeakRP[0] = GalA[p].Pos[0];
	    GalA[p].PosPeakRP[1] = GalA[p].Pos[1];
	    GalA[p].PosPeakRP[2] = GalA[p].Pos[2];
	    GalA[p].VelPeakRP[0] = GalA[p].Vel[0];
	    GalA[p].VelPeakRP[1] = GalA[p].Vel[1];
	    GalA[p].VelPeakRP[2] = GalA[p].Vel[2];
	  }
	}
	else if (nstep>0) {
	  /* For merged galaxies, copy the last value of the 
	     properties */
	  if (GalA[p].OrbitType[nstep]==3 || 
	      GalA[p].OrbitType[nstep]==5) {
	    GalA[p].OrbitMedianRhoICM[nstep] = 
	      GalA[p].OrbitMedianRhoICM[nstep-1];
	    GalA[p].OrbitVrel2ICM[nstep] = GalA[p].OrbitVrel2ICM[nstep-1];
	  }
	  
	  /* For type 0 and 1 galaxies, copy backwards the final value */
	  if (nstep==STEPS-1) {
	    if (GalA[p].OrbitType[nstep] <= 1) {
	      GalA[p].PeakRPMedianRhoICM = GalA[p].medianRhoICM;
	      GalA[p].PeakRPrelvel2ICM = GalA[p].relvel2ICM;
	      GalA[p].PosPeakRP[0] = GalA[p].Pos[0];
	      GalA[p].PosPeakRP[1] = GalA[p].Pos[1];
	      GalA[p].PosPeakRP[2] = GalA[p].Pos[2];
	      GalA[p].VelPeakRP[0] = GalA[p].Vel[0];
	      GalA[p].VelPeakRP[1] = GalA[p].Vel[1];
	      GalA[p].VelPeakRP[2] = GalA[p].Vel[2];
	      
	      for (istep=STEPS-2; istep<0; istep--) {
		GalA[p].OrbitMedianRhoICM[istep] = 
		  GalA[p].OrbitMedianRhoICM[STEPS-1];
		GalA[p].OrbitVrel2ICM[istep] = 
		  GalA[p].OrbitVrel2ICM[STEPS-1];
	      }
	    }
	  }
	} /* Close if cont = 0 and nstep > 0 */
      } /* Close STEPS */
#else 
      /* 
       * Procedure for special dump files (used when the positions of
       * galaxies are set within the semianalytic code, e.g. when using
       * DM particles to trace positions)
       */
      ngb_treefind(GalA[p].Pos, GasConstPartNum, 0);
      
      GalA[p].IPG_startRP = totnumngbRP;
      
      for (k = 0; k < NMAX; k++) {
	if (totnumngbRP >= totMaxNgb) {
	  printf("GAS SEARCH: maximum allowed total neighbour %d reached "
		 "for galaxy p=%d\n", totMaxNgb, p); fflush(stdout);
	  exit(EXIT_FAILURE);
	}
	
	IPG_RP[GalA[p].IPG_startRP + k] = IndexList[k];
	totnumngbRP++;
      }
      
      GalA[p].IPG_numberRP = GasConstPartNum;
      
      //printf("Computing ICM properties...\n"); fflush(stdout);
      get_ICM_properties(p,STEPS-1,snapshot);
#endif /*ORBITS*/
      //printf("\n");
    } /* Close loop over galaxies */
  } /* End of search for constant particle number */
    
  /* 
   * Find gas particles within a given search radius (enabled when 
   * GasSearchRadius > 0; in that case, GasConstPartNum is used as 
   * the minimum number of gas particles to find)
   */
  if (GasSearchRadius > 0) { 
    totMaxNgb = TOTMAXNEIGHBOURS_RP_FIXEDRADIUS;
    
    for (p=1; p<=TotNumGalA; p++) {

      /* For central galaxies we use a fraction of Rvir to avoid
	 searching the entire halo */
      /* TOMAS 2013-02-12: Rvir is in physical units, we should convert
	 to comoving distance here */
      if (GalA[p].Type == 0) {
        /*NeighbourSearchRadius = GasSearchRadiusGx0*GalA[p].Rvir;*/
        NeighbourSearchRadius = GasSearchRadiusGx0 *
	  GalA[p].Rvir*(1+ZZ[snapshot]);

	/* Check if radius is < 0, but do not exit on error if the search
	   radius is exactly 0 (the search radius can be set to 0 to 
	   force searching for a fixed number of particles) */
	if (NeighbourSearchRadius < 0.0) {
	  fprintf(stderr, "Error (find_gaspart_withtree): "
		  "NeighbourSearchRadius = %g < 0! - Stop\n", 
		  NeighbourSearchRadius);
	  fprintf(stderr, "p = %d (type 0), GasSearchRadiusGx0 = %g, "
		  "GalA[p].Rvir = %g, ZZ[%d] = %g\n", 
		  p, GasSearchRadiusGx0, GalA[p].Rvir, 
		  snapshot, ZZ[snapshot]); fflush(stderr);
	  exit(EXIT_FAILURE);
	}
      }
      else {  
	/* For satellites, the idea is to go beyond its own Rvir to
	   get the properties of what is outside its attached 
	   substructure */
        /*NeighbourSearchRadius = GasSearchRadius*GalA[p].Rvir;*/
	NeighbourSearchRadius = GasSearchRadius*GalA[p].Rvir*
	  (1+ZZ[snapshot]); 

	if (NeighbourSearchRadius < 0.0) {
	  fprintf(stderr, "Error (find_gaspart_withtree): "
		  "NeighbourSearchRadius = %g < 0! - Stop\n", 
		  NeighbourSearchRadius);
	  fprintf(stderr, "p = %d, GasSearchRadius = %g, GalA[p].Rvir = %g,"
		  " ZZ[%d] = %g\n", p, GasSearchRadius, GalA[p].Rvir, 
		  snapshot, ZZ[snapshot]); fflush(stderr);
	  exit(EXIT_FAILURE);
	}
      } /* Close setting the search radius */
/*
      printf("p = %d / %d, Type = %d, ParentGroup = %d\n", p, TotNumGalA, 
	     GalA[p].Type, GalA[p].ParentGroup);
      printf("Pos = [%g %g %g] (comoving)\n", GalA[p].Pos[0], 
	     GalA[p].Pos[1], GalA[p].Pos[2]);
      if (GalA[p].Type > 0)
	printf("Central galaxy = %d\n", GalA[p].CentralGal);
      printf("Radius within FOF group = %g (%g comoving)\n", 
	     GalA[p].Rgroup, GalA[p].Rgroup*(1+ZZ[snapshot]));
      printf("Virial radius of FOF group = %g (%g comoving)\n",
	     GalA[p].ParentGroupRvir, 
	     GalA[p].ParentGroupRvir*(1+ZZ[snapshot]));
      printf("Searching for particles within R = %f...\n",
	     NeighbourSearchRadius); fflush(stdout);
  */    
#ifdef ORBITS
      pram = 0.0;
      for (nstep=0; nstep<STEPS; nstep++) {
	/* These should be physical positions and velocities */
	/* TOMAS 2013-02-13 not anymore, now we convert to comoving units
	   as used within the snapshot */
	GalA[p].Pos[0] = GalA[p].OrbitX[nstep];
	GalA[p].Pos[1] = GalA[p].OrbitY[nstep];
	GalA[p].Pos[2] = GalA[p].OrbitZ[nstep];
	GalA[p].Vel[0] = GalA[p].OrbitVx[nstep];
	GalA[p].Vel[1] = GalA[p].OrbitVy[nstep];
	GalA[p].Vel[2] = GalA[p].OrbitVz[nstep];
	/*printf("Pos = [%g %g %g] (comoving from orbits)\n", GalA[p].Pos[0], 
	  GalA[p].Pos[1], GalA[p].Pos[2]);*/

	/* Keep updating only for surviving satellites */	
	cont = 0;
	if (nstep<STEPS-1) {
	  if (GalA[p].OrbitType[nstep]==2 || GalA[p].OrbitType[nstep]==4)
	    cont = 1;
	}
	else 
	  cont = 1;

	if (cont>0) {
	  /* Find neighbours within a given radius */
	  numngb = ngb_treefind_variable(GalA[p].Pos,NeighbourSearchRadius);

	  GalA[p].IPG_startRP = totnumngbRP;
	  if (GalA[p].IPG_startRP+numngb >= totMaxNgb) {
	    printf("GAS SEARCH: "
		   "Maximum allowed total neighbours %d reached"
		   " for galaxy p=%d, nstep=%d\n", totMaxNgb, p, nstep); 
	    fflush(stdout);
	    exit(EXIT_FAILURE);
	  }
	  
	  for (k=0; k<numngb; k++) {
	    /*if (p >= 450) {
	      printf("p = %d numngb = %d k = %d\n", p, numngb, k);
	      printf("GalA[p].IPG_startRP+k = %d IndexList[k] = %d\n",
	      GalA[p].IPG_startRP+k, IndexList[k]);
	      }*/
	    IPG_RP[GalA[p].IPG_startRP+k] = IndexList[k];
	    
	    /* TOMAS 2012-08-09: when calculating ICM properties along
	       a galaxy's orbit, this number becomes huge in the larger
	       simulations. I think it's because it only has to be
	       increased once in the STEPS */
	    if (nstep == STEPS-1) totnumngbRP++;	      
	  }
	  
	  GalA[p].IPG_numberRP = numngb;
	  
	  /* Option: if the number of neighbours found is less 
	     than GasMinPartNum, tries to find that many particles */
	  if (GalA[p].IPG_numberRP < GasConstPartNum) {
/*
	    if (nstep == 0) {
	      printf("Numngb = %d < GasConstPartNum = %d; "
		     "searching for additional neighbours...\n", 
		     GalA[p].IPG_numberRP, GasConstPartNum); 
	      fflush(stdout);
	    }
*/
	    totnumngbRP -= GalA[p].IPG_numberRP;
	    
	    ngb_treefind(GalA[p].Pos, GasConstPartNum, 0);
	    
	    GalA[p].IPG_startRP = totnumngbRP;
	    
	    if (GalA[p].IPG_startRP+GasConstPartNum >= totMaxNgb) {
	      printf("GAS SEARCH: maximum allowed total neighbours "
		     "%d reached for galaxy p=%d\n", totMaxNgb, p);
	      fflush(stdout);
	      exit(EXIT_FAILURE);
	    }
	    
	    for (k=0; k<GasConstPartNum; k++) {
	      IPG_RP[GalA[p].IPG_startRP + k]= IndexList[k];
	      totnumngbRP++;
	    }
	    
	    GalA[p].IPG_numberRP = GasConstPartNum;
	  }
	 
/* 
	  if (nstep == 0 || nstep == STEPS-1) {	    
	    printf("Computing ICM properties (nstep=%d)...\n",nstep); 
	    fflush(stdout);
	  }
*/	  
	  get_ICM_properties(p,nstep,snapshot);
	  
	  /* Find the peak RP along the orbit, and store the 
	     corresponding position and velocity */
	  if (GalA[p].OrbitMedianRhoICM[nstep]*
	      GalA[p].OrbitVrel2ICM[nstep] > pram) {
	    pram = GalA[p].OrbitMedianRhoICM[nstep]*
	      GalA[p].OrbitVrel2ICM[nstep];
	    GalA[p].PeakRPMedianRhoICM = GalA[p].OrbitMedianRhoICM[nstep];
	    GalA[p].PeakRPrelvel2ICM = GalA[p].OrbitVrel2ICM[nstep];
	    GalA[p].PosPeakRP[0] = GalA[p].Pos[0];
	    GalA[p].PosPeakRP[1] = GalA[p].Pos[1];
	    GalA[p].PosPeakRP[2] = GalA[p].Pos[2];
	    GalA[p].VelPeakRP[0] = GalA[p].Vel[0];
	    GalA[p].VelPeakRP[1] = GalA[p].Vel[1];
	    GalA[p].VelPeakRP[2] = GalA[p].Vel[2];
	  }
	}
	else if (nstep>0) {
	  /* For merged galaxies, copy the last value of the
	     properties */
	  if (GalA[p].OrbitType[nstep]==3 || 
	      GalA[p].OrbitType[nstep]==5) {
	    GalA[p].OrbitMedianRhoICM[nstep] = 
	      GalA[p].OrbitMedianRhoICM[nstep-1];
	    GalA[p].OrbitVrel2ICM[nstep] = GalA[p].OrbitVrel2ICM[nstep-1];
	  }

	  /* For type 0 and 1 galaxies, copy backwards the final value */
	  if (nstep==STEPS-1) {
	    if (GalA[p].OrbitType[nstep] <= 1) {
	      GalA[p].PeakRPMedianRhoICM = GalA[p].medianRhoICM;
	      GalA[p].PeakRPrelvel2ICM = GalA[p].relvel2ICM;
	      GalA[p].PosPeakRP[0] = GalA[p].Pos[0];
	      GalA[p].PosPeakRP[1] = GalA[p].Pos[1];
	      GalA[p].PosPeakRP[2] = GalA[p].Pos[2];
	      GalA[p].VelPeakRP[0] = GalA[p].Vel[0];
	      GalA[p].VelPeakRP[1] = GalA[p].Vel[1];
	      GalA[p].VelPeakRP[2] = GalA[p].Vel[2];
	      
	      for (istep=STEPS-2; istep<0; istep--) {
		GalA[p].OrbitMedianRhoICM[istep] = 
		  GalA[p].OrbitMedianRhoICM[STEPS-1];
		GalA[p].OrbitVrel2ICM[istep] = 
		  GalA[p].OrbitVrel2ICM[STEPS-1];
	      }
	    }
	  }
	} /* Close if cont = 0 and nstep > 0 */
      } /* Close STEPS */
#else 
      /*
       * Procedure for when special dump files are used
       */

      /* Find neighbours within a given radius */
      numngb = ngb_treefind_variable(GalA[p].Pos,NeighbourSearchRadius);
      GalA[p].IPG_startRP = totnumngbRP;
      if (GalA[p].IPG_startRP+numngb >= totMaxNgb) {
	printf("GAS SEARCH: Maximum allowed total neighbour %d reached"
	       " for galaxy p=%d\n", totMaxNgb, p); fflush(stdout);
	exit(EXIT_FAILURE);
      }
      
      for (k=0; k<numngb; k++) {
	/*if (p >= 450) {
	  printf("p = %d numngb = %d k = %d\n", p, numngb, k);
	  printf("GalA[p].IPG_startRP+k = %d IndexList[k] = %d\n",
	  GalA[p].IPG_startRP+k, IndexList[k]);
	  }*/
	IPG_RP[GalA[p].IPG_startRP+k] = IndexList[k];
	
	totnumngbRP++;
      }
      
      GalA[p].IPG_numberRP = numngb;
      
      /* Option: if the number of neighbours found is less 
	 than GasMinPartNum, tries to find that many particles */
      if (GalA[p].IPG_numberRP < GasConstPartNum) {
/*
	printf("Numngb = %d < GasConstPartNum = %d; "
	       "searching for additional neighbours...\n", 
	       GalA[p].IPG_numberRP, GasConstPartNum); fflush(stdout);
*/
	totnumngbRP -= GalA[p].IPG_numberRP;
	
	ngb_treefind(GalA[p].Pos, GasConstPartNum, 0);
	
	GalA[p].IPG_startRP = totnumngbRP;
	
	if (GalA[p].IPG_startRP+GasConstPartNum >= totMaxNgb) {
	  printf("GAS SEARCH: maximum allowed total neighbours "
		 "%d reached for galaxy p=%d\n", totMaxNgb, p);
	  exit(EXIT_FAILURE);
	}
	
	for (k=0; k<GasConstPartNum; k++) {
	  IPG_RP[GalA[p].IPG_startRP + k]= IndexList[k];
	  totnumngbRP++;
	}
	
	GalA[p].IPG_numberRP = GasConstPartNum;
      }
 
/*     
      printf("Computing ICM properties...\n"); 
      fflush(stdout);
  */    
      get_ICM_properties(p,STEPS-1,snapshot);
#endif /*ORBITS*/
      /*printf("\n");*/
    } /* End of loop over galaxies */
     //printf("TERMINO LOOP OVER GALAXIES\n");
    
  } /* End of search at fixed radius */
  
  ngb_treefree();
  
  return;
}


int comp_func(const void *a, const void *b)
{
  float *aa, *bb;
  
  aa = (float *)a;
  bb = (float *)b;
  
  if (*aa < *bb)
    return -1;
  if (*aa > *bb)
    return +1;
  
  return 0;
}


void check_neighbours(int p)
{
  int i, ig;
  float *rlist;
  float r, rmax, dx, dy, dz;
  
  //printf("Checking galaxy p=%d\n", p);
  
  rlist = vector(1, NumPart0);
  
  for (i = 1; i <= NumPart0; i++){
    dx = GalA[p].Pos[0] - P[i].Pos[0];
    dy = GalA[p].Pos[1] - P[i].Pos[1];
    dz = GalA[p].Pos[2] - P[i].Pos[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    
    rlist[i] = r;
  }
  
  qsort(&rlist[1], NumPart0, sizeof(float), comp_func);

  rmax = rlist[NMAX];
  /*
  for (i = 1; i <= NMAX; i++) {
    printf("r=%g\n", rlist[i]);
  }
  */
  for (i = 0; i < GalA[p].IPG_numberRP; i++) {
    ig = IPG_RP[GalA[p].IPG_startRP+i];
    
    if (ig < 1 || ig > NumPart0)
      exit(2);
    
    dx = GalA[p].Pos[0] - P[ig].Pos[0];
    dy = GalA[p].Pos[1] - P[ig].Pos[1];
    dz = GalA[p].Pos[2] - P[ig].Pos[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (r > rmax)
      printf("Strange! %g  %g  (%g|%g|%g) ig=%d  (%g|%g|%g) \n", 
	     rmax, r, GalA[p].Pos[0], GalA[p].Pos[1], GalA[p].Pos[2], ig, 
	     P[ig].Pos[0], P[ig].Pos[1], P[ig].Pos[2]);
  }
  
  free_vector(rlist, 1, NumPart0);

  return;
}


void min_max_dist(int p)
{
  int i, ig;
  float *rlist;
  float r, dx, dy, dz;
  

  rlist = vector(1, GalA[p].IPG_numberRP);
  
  for (i = 0; i < GalA[p].IPG_numberRP; i++) {
    ig = IPG_RP[GalA[p].IPG_startRP+i];
    
    dx = GalA[p].Pos[0] - P[ig].Pos[0];
    dy = GalA[p].Pos[1] - P[ig].Pos[1];
    dz = GalA[p].Pos[2] - P[ig].Pos[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    
    rlist[i+1] = r;
  }
  
  qsort(&rlist[1], GalA[p].IPG_numberRP, sizeof(float), comp_func);
  
  
  free_vector(rlist, 1, GalA[p].IPG_numberRP); 

  return;
}


/**
 * @brief Use the properties of the neighbouring particles found to 
 * determine local median density and velocity of the galaxy relative to the 
 * selected particles. 
 * @param p Index of galaxy
 * @param nstep Subdivision of snapshot (necessary for ORBITS only)
 * @param snapshot Current snapshot number
 *
 * Results are stored in elements @a medianRhoICM, @a meanRhoICM, 
 * @a IPG_numberRP and @a relvel2ICM of struct @a GalA (in PHYSICAL units).
 */
void get_ICM_properties(int p, int nstep, int snapshot)
{
  int i, j, ig, gascount, gascount_save, moveinlist, median;
  int *sortId;
  float *rlist;
  float meanRho, frac;
  /*float zcurr, H_of_a; TOMAS 2011-12-12 apparently not used here */
  float vxgas, vygas, vzgas;
  float vgroup;
  
  
  gascount = GalA[p].IPG_numberRP;
  
  /* 
   * Order gas particles by density to find the median, and so discard the 
   * high density tail of the distribution in the iterative process used to 
   * eliminate high density gas substructures. This allows the estimation 
   * of the ICM properties local to each galaxy
   */
  gascount_save = gascount;
  
  /* Proceed only if gas particles were found */
  if (gascount > 0) {    
    
    sortId = ivector(1, gascount);
    rlist  = vector(1, gascount);
    
    for (i = 1; i <= gascount; i++) {
      ig = IPG_RP[GalA[p].IPG_startRP + i - 1]; 
      sortId[i] = ig;
      
      /* PG.Rho is converted to physical density in the main program,
	 after reading particle properties */
      rlist[i] = PG[ig].Rho;
      /*printf("i: %d, sortId[i]: %d, rlist[i]: %g \n",i,sortId[i],rlist[i]);*/
    }
    
    sort2_flt_int(gascount, rlist, sortId); /* Reordering */
    
    /*for(i=1; i<= gascount; i++)
      printf("ORDENADO: i: %d, sortId[i]: %d, rlist[i]: %g \n",i,sortId[i],rlist[i]);*/
    
    /* Iteratively discard high density particles 
       (if Niterations is set to >0) */
    if (Niterations > 0) {
      for (j=0; j<Niterations; j++) {
	
	if (gascount >= 2)
	  median = (int)(gascount/2.);
	else
	  median = 1;  /* To avoid problems when only one 
			  particle is found and no minimum 
			  particle number is selected */
	
	GalA[p].medianRhoICM = PG[sortId[median]].Rho;
	
	/*printf("gx: %d, niter: %d, gascount:%d, median: %d, PG[sortId[median]].Rho: %g\n",p,j,gascount,median,PG[sortId[median]].Rho);*/
	/*printf("eliminando particulas de alta densidad ...\n\n");*/
	
	moveinlist = gascount;
	for (i=moveinlist; i>=1; i--) {
	  ig = sortId[i];
	  frac = PG[ig].Rho / GalA[p].medianRhoICM;
	  /*printf("i: %d, ig: %d, PG[ig].Rho: %g, frac: %g, gascount: %d\n",
	    i,ig, PG[ig].Rho,frac,gascount);*/

	  if (frac > Fracmax) {
	    gascount--;
	    /*printf("elimino....i: %d, ig: %d, PG[ig].Rho: %g, frac: %g, gascount: %d\n",i,ig, PG[ig].Rho,frac,gascount);*/
	  }
	  else {
	    /*printf("las que siguen no cumplen condicion de eliminar ...salgo del loop \n\nn");*/
	    break;
	  }
	}
	
	if (gascount <= GasMinimum) {
	  printf("GAS: WARNING: gx %d has ICM density estimated "
		 "only using %d particles\n", p, gascount);
	  fflush(stdout);
	  break;
	}
      } /* End iterations of discarding process */
    } /* End if Niterations > 0 */
  } /* Close if gascount > 0 */
  else { 
    /* No gas found */
    GalA[p].meanRhoICM   = 0.0;
    GalA[p].medianRhoICM = 0.0;
    printf("GAS: WARNING: no gas particles found around gx %d "
	   "with the given search criteria\n", p); fflush(stdout);
  }
  
  GalA[p].IPG_numberRP = gascount; /* Save how many particles were
				      used to compute ICM properties */

  /* 
   * Find mean and median density, and velocity of galaxy relative 
   * to ICM gas
   */
  vxgas = 0.;
  vygas = 0.;
  vzgas = 0.;
  meanRho = 0.;
  
  if (gascount > 0) {
    for (i = 1; i <= gascount; i++) {
      ig = sortId[i]; /* Evaluate in the list of indexes 
			 ordered by density */
    
      /* These velocities and density are converted to physical units
	 in the main program, after reading the particle data. Here we
	 keep the comoving units used in the snapshot */
      vxgas   += P[ig].Vel[0];
      vygas   += P[ig].Vel[1];
      vzgas   += P[ig].Vel[2];
      meanRho += PG[ig].Rho;
      
      /* printf("i: %d, ig: %d, P[ig].Vel: %g|%g|%g, PG[ig].Rho:%g \n\n",
	 i,ig,P[ig].Vel[0],P[ig].Vel[1],P[ig].Vel[2],PG[ig].Rho); */
    }
    
    vxgas /= gascount;
    vygas /= gascount;
    vzgas /= gascount;
    
    meanRho /= gascount;
    
    if (gascount >= 2)
      median =(int)(gascount/2.);
    else
      median = 1;   /* To avoid errors when only one particle 
		       is found and no minimum particle number 
		       is set */
    
    GalA[p].meanRhoICM   = meanRho;
    GalA[p].medianRhoICM = PG[sortId[median]].Rho;
    
    GalA[p].relvel2ICM = (vxgas-GalA[p].Vel[0])*(vxgas-GalA[p].Vel[0]) +
      (vygas-GalA[p].Vel[1])*(vygas-GalA[p].Vel[1]) +
      (vzgas-GalA[p].Vel[2])*(vzgas-GalA[p].Vel[2]);
    
    /*printf("RAM PRESSURE GX: \n"
	   "GalA[p].Type: %d\n"
	   "GalA[p].Rvir: %g\n"
	   "GalA[p].Pos: %g %g %g\n",
	   "GalA[p].Vel: %g %g %g\n",
	   GalA[p].Type, GalA[p].Rvir, 
	   GalA[p].Pos[0], GalA[p].Pos[1], GalA[p].Pos[2],
	   GalA[p].Vel[0], GalA[p].Vel[1], GalA[p].Vel[2]);
    printf("RAM PRESSURE: p: %d, GalA[p].Type: %d, vxgas: %g, vygas: %g, vzgas: %g, GalA[p].relvel2ICM:%g \n\n",
	   p,GalA[p].Type,GalA[p].vxgasICM,GalA[p].vygasICM,GalA[p].vzgasICM,GalA[p].relvel2ICM); 
    */
    
    free_vector(rlist, 1, gascount_save);
    free_ivector(sortId, 1, gascount_save);
  }
  else {
    /* If no gas particles were found, just store the galaxy's velocity 
       (squared); densities were set to 0 above */
    GalA[p].relvel2ICM = (GalA[p].Vel[0]*GalA[p].Vel[0]) +
      (GalA[p].Vel[1]*GalA[p].Vel[1]) + (GalA[p].Vel[2]*GalA[p].Vel[2]);
  }

#ifdef ORBITS
  /* For ORBITS, store the values found along the orbit */
  /*if (GalA[p].Type <= 1) {
    nstep = 0;
    for (nstep=0; nstep<STEPS; nstep++) {
      GalA[p].OrbitMedianRhoICM[nstep] = GalA[p].medianRhoICM;
      GalA[p].OrbitVrel2ICM[nstep] = GalA[p].relvel2ICM;
    }
  }
  else {*/
    GalA[p].OrbitMedianRhoICM[nstep] = GalA[p].medianRhoICM;
    GalA[p].OrbitVrel2ICM[nstep] = GalA[p].relvel2ICM;
    /*}*/
#endif

  /*
   * Now calculate the density of gas using an analytic profile, 
   * assuming that the ICM is distributed parallel to the DM.
   * These densities are saved now in physical units. 
   * Note that we cannot directly calculate for type 0 galaxies (r = 0);
   * in those cases, use the density at 0.1rvir 
   * TOMAS 2013-02-13: in special dump files, there are type 2 galaxies
   * which have Rgroup = 0... merged galaxies? In ORBITS this happens for
   * type 3 galaxies
   */
  if (GalA[p].Rgroup > 0) { 
    GalA[p].RhoICM_SIS = BaryonFrac*density_SIS(GalA[p].Rgroup,
						GalA[p].ParentGroupMvir,
						GalA[p].ParentGroupRvir);
    if (GalA[p].RhoICM_SIS < 0) {
      fprintf(stderr, "Error (get_icm_properties): "
	      "function density_SIS returned an error - Exit\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

#ifndef ORBITS
    /* TOMAS 2013-02-13: halo concentration is not yet available for
       option ORBITS */
    GalA[p].RhoICM_NFW = BaryonFrac*density_NFW(GalA[p].Rgroup,
						GalA[p].ParentGroupMvir,
						GalA[p].ParentGroupRvir,
						GalA[p].ParentGroupCvir,
						ZZ[snapshot]);
    if (GalA[p].RhoICM_NFW < 0) {
      fprintf(stderr, "Error (get_icm_properties): "
	      "function density_NFW returned an error - Exit\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
#endif
  }
  else {
    GalA[p].RhoICM_SIS = BaryonFrac*density_SIS(0.1*GalA[p].ParentGroupRvir,
						GalA[p].ParentGroupMvir,
						GalA[p].ParentGroupRvir);
    if (GalA[p].RhoICM_SIS < 0) {
      fprintf(stderr, "Error (get_icm_properties): "
	      "function density_SIS returned an error - Exit\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

#ifndef ORBITS
    GalA[p].RhoICM_NFW = BaryonFrac*density_NFW(0.1*GalA[p].ParentGroupRvir,
						GalA[p].ParentGroupMvir,
						GalA[p].ParentGroupRvir,
						GalA[p].ParentGroupCvir,
						ZZ[snapshot]);
    if (GalA[p].RhoICM_SIS < 0) {
      fprintf(stderr, "Error (get_icm_properties): "
	      "function density_SIS returned an error - Exit\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
#endif
  }

  if (nstep == STEPS-1) {
/*
    printf("Median ICM density: %g (%g h^2 g cm^-3 physical)\n",
	   GalA[p].medianRhoICM, 
	   GalA[p].medianRhoICM*pow(1+ZZ[snapshot],3)*UnitDensity_in_cgs);
    printf("ICM density (SIS):  %g (%g h^2 g cm^-3 physical)\n",
	   GalA[p].RhoICM_SIS/pow(1+ZZ[snapshot],3), 
	   GalA[p].RhoICM_SIS*UnitDensity_in_cgs);
    printf("ICM density (NFW):  %g (%g h^2 g cm^-3 physical)\n",
	   GalA[p].RhoICM_NFW/pow(1+ZZ[snapshot],3), 
	   GalA[p].RhoICM_NFW*UnitDensity_in_cgs);
    printf("Relative velocity: %g (%g km s^-1 physical)\n",
	   sqrt(GalA[p].relvel2ICM), 
	   sqrt(GalA[p].relvel2ICM/(1+ZZ[snapshot])));
*/

    if (GalA[p].Type > 0) {
      vgroup = sqrt((GalA[p].Vel[0]-GalA[GalA[p].CentralGal].Vel[0])*
		    (GalA[p].Vel[0]-GalA[GalA[p].CentralGal].Vel[0]) + 
		    (GalA[p].Vel[1]-GalA[GalA[p].CentralGal].Vel[1])*
		    (GalA[p].Vel[1]-GalA[GalA[p].CentralGal].Vel[1]) + 
		    (GalA[p].Vel[2]-GalA[GalA[p].CentralGal].Vel[2])*
		    (GalA[p].Vel[2]-GalA[GalA[p].CentralGal].Vel[2]));
/*
      printf("Velocity within FOF group: %g (%g km s^-1 physical)\n", 
	     vgroup, vgroup/sqrt(1+ZZ[snapshot]));
*/
    }
/*
    printf("Virial velocity: %g km s^-1 (physical)\n", 
	   GalA[p].ParentGroupVvir*UnitVelocity_in_cm_per_s/1e5);
    fflush(stdout);
*/
  }
	 
  return;
}
