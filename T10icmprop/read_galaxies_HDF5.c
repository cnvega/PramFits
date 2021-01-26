/**
 * @file read_galaxies_HDF5.c
 * @brief Read input galaxies data in HDF5 format
 * @author Tomas E. Tecce
 */
#ifdef HDF5OUTPUT
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

#include <hdf5.h>
#include <hdf5_hl.h>


/**
 * @brief Read galaxies data from special dump files created by SAG
 * @param snap Number of the snapshot to process
 * @param fname COMPLETAR
 */
void get_galaxies_data_HDF5(int snap, char *fname)
{
  int p;
  hid_t file_id;
  hsize_t dims[2];
  int *ibuffer;
  float *buffer, *buffer3d;
  herr_t status;
  size_t i, j, nrow, n_values;


  /* Open special dump file from SAG */
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Get the total number of galaxies, stored as an attribute */
  status = H5LTget_attribute_int(file_id, "/", "Numgal", &TotNumGalA);

  /* Get the random number seeds */
  status = H5LTget_attribute_int(file_id, "/", "Seed", &Seed);

  /* TOMAS: not used anymore 2013-02-13 */
  /*status = H5LTget_attribute_int(file_id, "/", "Sx", &Sx);
  status = H5LTget_attribute_int(file_id, "/", "Sy", &Sy);
  status = H5LTget_attribute_int(file_id, "/", "Sz", &Sz);
  status = H5LTget_attribute_int(file_id, "/", "Sc", &Sc);*/

  /* Proceed if galaxies were found */
  if (TotNumGalA > 0) {

    /* Allocate memory for read buffers */
    buffer    = malloc(TotNumGalA*sizeof(float));
    buffer3d  = malloc(3*TotNumGalA*sizeof(float));
    ibuffer   = malloc(TotNumGalA*sizeof(int));
    
    /* Read dataset */
    status = H5LTread_dataset_int(file_id, "/ParentGroup", ibuffer);
    
    /* Get the dimensions of the dataset */
    status = H5LTget_dataset_info(file_id, "/ParentGroup", dims, NULL,NULL);
    
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
    
    status = H5LTread_dataset_float(file_id, "/Mvir", buffer);
    status = H5LTget_dataset_info(file_id, "/Mvir", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Mvir = buffer[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/Rvir", buffer);
    status = H5LTget_dataset_info(file_id, "/Rvir", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Rvir = buffer[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/Rscale", buffer);
    status = H5LTget_dataset_info(file_id, "/Rscale", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Rscale = buffer[i*nrow + j];
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
    
    status = H5LTread_dataset_int(file_id, "/CentralGal", ibuffer);
    status = H5LTget_dataset_info(file_id, "/CentralGal", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].CentralGal = ibuffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/ParentGroupMvir", buffer);
    status = H5LTget_dataset_info(file_id, "/ParentGroupMvir", dims, 
				  NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].ParentGroupMvir = buffer[i*nrow + j];
    }
    
    status = H5LTread_dataset_float(file_id, "/ParentGroupRvir", buffer);
    status = H5LTget_dataset_info(file_id, "/ParentGroupRvir", dims, 
				  NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].ParentGroupRvir = buffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/ParentGroupCvir", buffer);
    status = H5LTget_dataset_info(file_id, "/ParentGroupCvir", dims, 
				  NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].ParentGroupCvir = buffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/ParentGroupVvir", buffer);
    status = H5LTget_dataset_info(file_id, "/ParentGroupVvir", dims, 
				  NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].ParentGroupVvir = buffer[i*nrow + j];
    }

    status = H5LTread_dataset_float(file_id, "/Rgroup", buffer);
    status = H5LTget_dataset_info(file_id, "/Rgroup", dims, NULL, NULL);
    n_values = (size_t)(dims[0] * dims[1]);
    nrow = (size_t)dims[1];
    for (i=0; i<n_values/nrow; i++ ) {
      for (j=0; j<nrow; j++)
	GalA[i+1].Rgroup = buffer[i*nrow + j];
    }

    free(buffer);
    free(buffer3d);
    free(ibuffer);
    
    /* Close file */
    status = H5Fclose(file_id);
    
    printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (comoving)\n",
	   GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
    printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (comoving)\n",
	   GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);

    /* TOMAS 2013-02-12
       Only convert the final results to physical coordinates */
      
    /* Conversion of position and velocity to physical coordinates */
    /*for (p=1; p<=TotNumGalA; p++) {
      for (i=0; i<3; i++) {
	GalA[p].PosComov[i] = GalA[p].Pos[i];
	GalA[p].VelComov[i] = GalA[p].Vel[i];
	GalA[p].Pos[i] /= 1.+Zcurr;
	GalA[p].Vel[i] /= sqrt(1.+Zcurr);*/      /* To peculiar velocity */
	/* Add Hubble flow */
	/*GalA[p].Vel[i] += Hubble_of_z*GalA[p].Pos[i];*/
    /* }  
    }
    printf("Pos: x(0) = %f, y(0) = %f, z(0) = %f (physical)\n",
	   GalA[1].Pos[0], GalA[1].Pos[1], GalA[1].Pos[2]);
    printf("Vel: Vx(0) = %f, Vy(0) = %f, Vz(0) = %f (physical)\n",
	   GalA[1].Vel[0], GalA[1].Vel[1], GalA[1].Vel[2]);*/
  
  } /* Close if galaxies were found */
  else
    printf("No galaxies found in this snapshot\n");
    
  return;
}
#endif
