/**
 * @file dumpHDF5_icm.c
 * @brief Function to write output data in HDF5 format
 * @author Tomas E. Tecce
 */
#ifdef HDF5OUTPUT
#include <hdf5.h>
#include <hdf5_hl.h>

#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Write ICM data to output files, in HDF5 format
 * @author Tomas E. Tecce
 *
 * Remember to reflect changes to this function in the function 
 * get_icm_properties in SAG
 * Pos and Vel are stored in comoving units (recovered from elements
 * PosComov and VelComov, respectively), to match the units in the 
 * snapshots and in SAG. The same for the position and velocity at the
 * RP peak. Densities and velocities relative to ICM are stored in PHYSICAL 
 * units.
 * If ORBITS is enabled, the values stored for medianRhoICM and relvel2ICM
 * are those that correspond to the maximum RP along the orbit for type 2
 * satellites (or merged galaxies).
 * Datasets written in output:
 *   /ParentGroup          INT
 *   /Type                 INT
 *   /Pos                  FLOAT x 3
 *   /Vel                  FLOAT x 3
 *   /PosPeakRP            FLOAT x 3   (for ORBITS only)
 *   /VelPeakRP            FLOAT x 3   (for ORBITS only)
 *   /meanRhoICM           FLOAT     
 *   /medianRhoICM         FLOAT
 *   /relvel2ICM           FLOAT
 *   /PeakRPMedianRhoICM   FLOAT       (for ORBITS only)
 *   /PeakRPrelvel2ICM     FLOAT       (for ORBITS only)
 *   /IPG_numberRP         INT
 *   /RhoICM_SIS           FLOAT
 *   /RhoICM_NFW           FLOAT
 */
void dumpHDF5_icm(char *filename, int snapshot, float zcurr)
{
  int i, n, p, check, index, indrp;
  int count_total, count_total_satellites;
  int count_cluster, count_cluster_satellites;
  char desc[256], label[256];
  float pram, factor;
  double buf;
  int *ibuffer;
  float *buffer, *buffer3d;
  hsize_t dims[2]={TotNumGalA,1};
  hsize_t dim3d[2]={TotNumGalA,3};
  hid_t file_id;
  herr_t status;

  factor = 1.0 + zcurr;

  buffer   = malloc(TotNumGalA*sizeof(float));
  buffer3d = malloc(3*TotNumGalA*sizeof(float));
  ibuffer  = malloc(TotNumGalA*sizeof(int)); 

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
    
  /* Create a new HDF5 file using default properties */
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* 
   * Write file attributes 
   * Store the parameters used and the code and compilation options
   */
  status = dumpHDF5_file_attr(file_id, snapshot, Zcurr);

  /*
   * Write the datasets
   */
  sprintf(label, "/ParentGroup");
  sprintf(desc, "Index of parent FOF group");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalA[p].ParentGroup;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);
  
  sprintf(label, "/Type");
  sprintf(desc, "Galaxy type");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalA[p].Type;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);
  
  /* Original position in comoving coordinates, retrieved from PosComov */
  /* TOMAS 2013-02-12: now we do not convert the position or the 
     velocity */
  sprintf(label, "/Pos");
  sprintf(desc, "Position of galaxy");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      /*buffer3d[index] = GalA[p].PosComov[i];*/
      buffer3d[index] = GalA[p].Pos[i];
      if (p == 1) {
	/*printf("GalA[1].PosComov[%d] = %g\n", i, GalA[p].PosComov[i]);*/
	printf("GalA[1].Pos[%d] = %g\n", i, GalA[p].Pos[i]);
      }
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In comoving units; to convert to "
				   "physical units divide by (1+z). "
				   "Unit length given in h^-1 cm.");

  /* Original velocity in comoving coordinates as written in Gadget, 
     retrieved from VelComov */
  /* TOMAS 2013-02-12: now we do not convert the position or the 
     velocity */
  sprintf(label, "/Vel");
  sprintf(desc, "Velocity of galaxy");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      /*buffer3d[index] = GalA[p].VelComov[i];*/
      buffer3d[index] = GalA[p].Vel[i];
      if (p == 1) {
	/*printf("GalA[1].VelComov[%d] = %g\n", i, GalA[p].VelComov[i]);*/
	printf("GalA[1].Vel[%d] = %g\n", i, GalA[p].Vel[i]);
      }
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In comoving units as written in Gadget."
				   " To convert to physical units divide "
				   "by sqrt(1+z). Unit velocity given in "
				   "cm s^-1.");

#ifdef ORBITS
  /* Position at peak RP, converted to comoving coordinates */
  sprintf(label, "/PosPeakRP");
  sprintf(desc, "Position of galaxy at the point of peak RP");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      buffer3d[index] = GalA[p].PosPeakRP[i]*factor;
      if (p == 1)
	printf("GalA[1].PosPeakRP[%d] = %g\n", i, GalA[p].PosPeakRP[i]);
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In comoving units; to convert to "
				   "physical units divide by (1+z). "
				   "Unit length given in h^-1 cm.");

  /* Velocity at peak RP, converted to comoving coordinates as written in 
     Gadget */
  sprintf(label, "/VelPeakRP");
  sprintf(desc, "Velocity of galaxy at the point of peak RP");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      buffer3d[index] = GalA[p].VelPeakRP[i]*sqrt(factor);
      if (p == 1)
	printf("GalA[1].VelPeakRP[%d] = %g\n", i, GalA[p].VelPeakRP[i]);
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In comoving units as written in Gadget."
				   " To convert to physical units divide "
				   "by sqrt(1+z). Unit velocity given in "
				   "cm s^-1.");
#endif
  
  sprintf(label, "/meanRhoICM");
  sprintf(desc, "Mean density of the ICM");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].meanRhoICM*pow(factor,3);
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitDensity_in_cgs, 
				   "Unit density given in h^2 g cm^-3 "
				   "(physical units).");

  sprintf(label, "/medianRhoICM");
  sprintf(desc, "Median density of the ICM");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].medianRhoICM*pow(factor,3);
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitDensity_in_cgs, 
				   "Unit density given in h^2 g cm^-3 "
				   "(physical units).");
        
  sprintf(label, "/relvel2ICM");
  sprintf(desc, "Velocity of the galaxy relative to the ICM, squared");
  for (p=1; p<=TotNumGalA; p++) {
    /* This is the square of a velocity, so the conversion factor
       is also squared */
    buffer[p-1] = GalA[p].relvel2ICM/factor;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  buf = UnitVelocity_in_cm_per_s*UnitVelocity_in_cm_per_s;
  status = set_float_variable_attr(file_id, label, desc, &buf, 
				   "Unit squared velocity in cm^2 s^-2 "
				   "(physical units).");

#ifdef ORBITS
  sprintf(label, "/PeakRPMedianRhoICM");
  sprintf(desc, "Median density of the ICM at the point of peak RP");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].PeakRPMedianRhoICM;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitDensity_in_cgs, 
				   "Unit density given in h^2 g cm^-3 "
				   "(physical units).");
        
  sprintf(label, "/PeakRPrelvel2ICM");
  sprintf(desc, "Velocity of the galaxy relative to the ICM, squared, "
	  "at the point of peak RP");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].PeakRPrelvel2ICM;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  buf = UnitVelocity_in_cm_per_s*UnitVelocity_in_cm_per_s;
  status = set_float_variable_attr(file_id, label, desc, &buf, 
				   "Unit squared velocity in cm^2 s^-2 "
				   "(physical units).");
#endif
  
  sprintf(label, "/IPG_numberRP");
  sprintf(desc, "Number of gas particles used to determine ICM properties");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalA[p].IPG_numberRP;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  /* Analytical densities are already in physical units */
  sprintf(label, "/RhoICM_SIS");
  sprintf(desc, "Density of the ICM using a SIS profile");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].RhoICM_SIS;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitDensity_in_cgs, 
				   "Unit density given in h^2 g cm^-3 "
				   "(physical units).");

  sprintf(label, "/RhoICM_NFW");
  sprintf(desc, "Density of the ICM using a NFW profile");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalA[p].RhoICM_NFW;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitDensity_in_cgs, 
				   "Unit density given in h^2 g cm^-3 "
				   "(physical units).");

  free(buffer);
  free(buffer3d);
  free(ibuffer);

  status = H5Fclose(file_id);

  return;
}


/**
 * @brief Set the attributes of the output variables
 *
 * For each attribute set, store a description, the value of the 
 * corresponding unit in the code, and a comment. 
 */
int set_float_variable_attr(hid_t file_id, char *label, char *desc, 
			    double *unit, char *comment)
{
  herr_t status;

  status = H5LTset_attribute_string(file_id, label, "Description", desc);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Description'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  status = H5LTset_attribute_double(file_id, label, "Unit", unit, 1);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Unit'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  status = H5LTset_attribute_string(file_id, label, "Comment", comment);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Comment'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


/**
 * @brief Check for errors in creating data
 */
void checkwrite_HDF5(herr_t status)
{
  if (status < 0) {
    fprintf(stderr, "Error (checkwrite_HDF5): cannot create attribute - "
	    "Exit|\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return;
}


/**
 * @brief Writes the parameters used in the run as attributes to the root
 * group of the output file
 */
herr_t dumpHDF5_file_attr(hid_t file_id, int snapshot, float zcurr)
{
  herr_t status;
  char hostname[256];
  time_t currtime;
  struct tm *loctime;
  float buf;


  status = H5LTset_attribute_string(file_id, "/", "Format", "ICM-HDF5");
  checkwrite_HDF5(status);

  gethostname(hostname, 256);
  status = H5LTset_attribute_string(file_id, "/", "System", hostname);
  checkwrite_HDF5(status);
  
  currtime = time(NULL);
  loctime = localtime(&currtime);
  status = H5LTset_attribute_string(file_id, "/", "Created", 
				    asctime(loctime));
  
  status = H5LTset_attribute_int(file_id, "/", "Snapshot", &snapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Numgal", &TotNumGalA, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "Redshift", &zcurr, 1);
  checkwrite_HDF5(status);

  /*
   * Write the parameters used in the run
   */
  status = H5LTset_attribute_int(file_id, "/", "Seed", &Seed, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Sx", &Sx, 1);
  checkwrite_HDF5(status);
  
  status = H5LTset_attribute_int(file_id, "/", "Sy", &Sy, 1);
  checkwrite_HDF5(status);
 
  status = H5LTset_attribute_int(file_id, "/", "Sz", &Sz, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Sc", &Sc, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path1", path1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path", path);
  checkwrite_HDF5(status);
  
  status = H5LTset_attribute_string(file_id, "/", "Path2", path2);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path3", path3);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "FileBasename", 
				    filebasename);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Files", &Files, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "MaxSnapshot", 
				 &MaxSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Omega", &Omega, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "OmegaLambda", 
				    &OmegaLambda, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Hubble_h", 
				    &Hubble_h, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "H0", &H0, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "G", &G, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "GasSearchRadiusGx0", 
				   &GasSearchRadiusGx0, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "GasSearchRadius", 
				   &GasSearchRadius, 1);
  checkwrite_HDF5(status);
  
  status = H5LTset_attribute_int(file_id, "/", "GasConstPartNum", 
				 &GasConstPartNum, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "GasMinimum", 
				 &GasMinimum, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Niterations", 
				 &Niterations, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "Fracmax", &Fracmax, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Identifier", 
				   Identifier);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "TimeBetSnapshot", 
				    &TimeBetSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "TimeOfFirstSnapshot", 
				    &TimeOfFirstSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "UnitLength_in_cm", 
				    &UnitLength_in_cm, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "UnitMass_in_g", 
				    &UnitMass_in_g, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", 
				    "UnitVelocity_in_cm_per_s", 
				    &UnitVelocity_in_cm_per_s, 1);
  checkwrite_HDF5(status);


  /*
   * Code options
   */
  status = H5LTset_attribute_uint(file_id, "/", "OutputListOn", 
				 &OutputListOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "ReorderOn", 
				 &ReorderOn, 1);
  checkwrite_HDF5(status);

  /*status = H5LTset_attribute_uint(file_id, "/", "ComovingToPhysicalOn", 
				  &ComovingToPhysicalOn, 1);
				  checkwrite_HDF5(status);*/
#ifdef ORBITS
  status = H5LTset_attribute_string(file_id, "/", "ORBITS", "YES");
#else
  status = H5LTset_attribute_string(file_id, "/", "ORBITS", "NO");
#endif
  checkwrite_HDF5(status);
  
  return status;
}
#endif
