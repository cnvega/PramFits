/**
 * @file io_orbits.c
 * @brief Read galaxy input data from orbit files
 * @author Tomas E. Tecce
 */
#ifdef HDF5OUTPUT
#ifdef ORBITS
#include <hdf5.h>
#include <hdf5_hl.h>

#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Reads galaxy orbits data from file, to set the positions of 
 * satellite galaxies and the mass loss of their host subhaloes via tidal 
 * stripping
 * @param suborbits_fname Name of file with orbits data
 * @param snap Current snapshot number
 */
void read_suborbits_A(char *suborbits_fname, int snap)
{
  int i, j, ibuf, p;
  int *ibuffer, *ibufferN;
  float factor;
  float *buffer, *buffer3d, *bufferN;
  hid_t file_id;
  herr_t status;
  hsize_t dims[2];
  size_t nrow, n_values;

  factor = 1.0 + ZZ[snap];

  printf("Reading galaxy orbits from file '%s'...\n", suborbits_fname);
  fflush(stdout);

  /* Open orbits file (HDF5 format only) */
  file_id = H5Fopen(suborbits_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  /* Sanity check: see if the snapshot matches */
  status = H5LTget_attribute_int(file_id, "/", "Snapshot", &ibuf);
  if (ibuf != snap) {
    fprintf(stderr, "Error (read_suborbits_A): snapshot number does not "
	    "match (orbits %d, SAG %d) - Exit\n", ibuf, snap);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  /* Read the ParentGroup dataset to get the number of galaxies */
  status = H5LTget_dataset_info(file_id, "/ParentGroup", dims, NULL, NULL);
  TotNumGalA = (size_t)(dims[0]*dims[1]);
  
  /* Proceed if galaxies were found */
  if (TotNumGalA > 0) {
    printf("Reading orbits for %d galaxies.\n\n", TotNumGalA); 
    fflush(stdout);

    /* Allocate memory for reading buffers */
    ibuffer  = malloc(TotNumGalA*sizeof(int));
    ibufferN = malloc(TotNumGalA*STEPS*sizeof(int));
    buffer   = malloc(TotNumGalA*sizeof(float));
    buffer3d = malloc(3*TotNumGalA*sizeof(float));
    bufferN  = malloc(TotNumGalA*STEPS*sizeof(float));

    /* Read dataset */
    status = H5LTread_dataset_int(file_id, "/ParentGroup", ibuffer);
    /* Save it by rows */
    /* Remember that galaxies are indexed from 1 to N */
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].ParentGroup = ibuffer[i*nrow + j];
    }

    status = H5LTread_dataset_int(file_id, "/Type", ibuffer);
    status = H5LTget_dataset_info(file_id, "/Type", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Type = ibuffer[i*nrow + j];
    }

    status = H5LTread_dataset_int(file_id, "/CentralGal", ibuffer);
    status = H5LTget_dataset_info(file_id, "/CentralGal",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].CentralGal = ibuffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Pos", buffer3d);
    status = H5LTget_dataset_info(file_id, "/Pos", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Pos[j] = buffer3d[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/Vel", buffer3d);
    status = H5LTget_dataset_info(file_id, "/Vel", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Vel[j] = buffer3d[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Rvir", buffer);
    status = H5LTget_dataset_info(file_id, "/Rvir", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Rvir = buffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Mvir", buffer);
    status = H5LTget_dataset_info(file_id, "/Mvir", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Mvir = buffer[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitX", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitX",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].OrbitX[j] = bufferN[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitY", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitY",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].OrbitY[j] = bufferN[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitZ", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitZ",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].OrbitZ[j] = bufferN[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitVx", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitVx",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++) 
	GalA[i+1].OrbitVx[j] = bufferN[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitVy", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitVy",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].OrbitVy[j] = bufferN[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Orbits/OrbitVz", bufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitVz",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++) 
	GalA[i+1].OrbitVz[j] = bufferN[i*nrow + j];
    }
    
    status = H5LTread_dataset_int(file_id, "/Orbits/OrbitType", ibufferN);
    status = H5LTget_dataset_info(file_id,"/Orbits/OrbitType",dims,NULL,NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<(int)dims[0]; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].OrbitType[j] = ibufferN[i*nrow + j];
    }

    /* Conversion of position and velocity to physical coordinates
       (Pos and Vel; the others already are in physical units) */
    /* TOMAS 2013-02-12:
       Use the comoving coordinates as read from the snapshot files */
    /*for (p=1; p<=TotNumGalA; p++) {
      for (i=0; i<3; i++) {
	GalA[p].PosComov[i] = GalA[p].Pos[i];
	GalA[p].VelComov[i] = GalA[p].Vel[i];
	GalA[p].Pos[i] /= 1.+Zcurr;
	GalA[p].Vel[i] /= sqrt(1.+Zcurr);*/    /* To peculiar velocity */
	/* Add Hubble flow */
	/*GalA[p].Vel[i] += Hubble_of_z*GalA[p].Pos[i];*/
    /* }  
    }*/
    printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (comoving)\n",
	   GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
    printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (comoving)\n",
	   GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);
    /*printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (physical)\n",
	   GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
    printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (physical)\n",
    GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);*/
    
    /* For ORBITS, OrbitX, OrbitVx etc store positions and velocities
       along the orbit relative to the central galaxy. To search for
       neighbouring gas particles, we need to transform them to the
       simulation reference frame (and transform to comoving units) */
    for (p=1; p<=TotNumGalA; p++) {
      /* Before transforming, calculate the current groupcentric radius */
      GalA[p].Rgroup = GalA[p].OrbitX[STEPS-1]*GalA[p].OrbitX[STEPS-1] +
	GalA[p].OrbitY[STEPS-1]*GalA[p].OrbitY[STEPS-1] +
	GalA[p].OrbitZ[STEPS-1]*GalA[p].OrbitZ[STEPS-1];
      GalA[p].Rgroup = sqrt(GalA[p].Rgroup); /* Already in physical units */

      /* TOMAS 2013-02-13:
	 This has to be done for all galaxies. Otherwise, when using
	 Orbit(X,Y,Z) to assign positions, a central galaxy will have its
	 relative position to itself i.e. [0,0,0] */
      for (i=0; i<STEPS; i++) {
	GalA[p].OrbitX[i] = GalA[GalA[p].CentralGal].Pos[0] + 
	  GalA[p].OrbitX[i]*factor;
	GalA[p].OrbitY[i] = GalA[GalA[p].CentralGal].Pos[1] + 
	  GalA[p].OrbitY[i]*factor;
	GalA[p].OrbitZ[i] = GalA[GalA[p].CentralGal].Pos[2] +
	  GalA[p].OrbitZ[i]*factor;
	
	GalA[p].OrbitVx[i] = GalA[GalA[p].CentralGal].Vel[0] + 
	  GalA[p].OrbitVx[i]*sqrt(factor);
	GalA[p].OrbitVy[i] = GalA[GalA[p].CentralGal].Vel[1] +
	  GalA[p].OrbitVy[i]*sqrt(factor);
	GalA[p].OrbitVz[i] = GalA[GalA[p].CentralGal].Vel[2] +
	  GalA[p].OrbitVz[i]*sqrt(factor);
      }

      /* Set other quantities */
      GalA[p].ParentGroupRvir = GalA[GalA[p].CentralGal].Rvir;
      GalA[p].ParentGroupMvir = GalA[GalA[p].CentralGal].Mvir;
      GalA[p].ParentGroupVvir = sqrt(G*GalA[GalA[p].CentralGal].Mvir/
				     GalA[GalA[p].CentralGal].Rvir);
    }

    free(buffer);
    free(ibuffer);
    free(ibufferN);
    free(bufferN);
    free(buffer3d);
  } /* Close if TotNumGalA > 0 */

  status = H5Fclose(file_id);
  
  printf("Suborbits catalogue '%s' read.\n", suborbits_fname);
  if (TotNumGalA == 0) printf("No galaxies found.\n");

  return;
}
#endif /* ORBITS */
#endif /* HDF5OUTPUT */
  
