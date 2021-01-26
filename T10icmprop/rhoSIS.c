/**
 * @file rhoSIS.c
 * @brief Function to calculate density (using a singular isothermal sphere 
 * profile)
 * @author Tomas E. Tecce
 * @date 2013-02-12
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

/**
 * @brief Calculate density as a function of radius for a singular 
 * isothermal sphere profile
 * @param r Radius in code units
 * @param mhalo Host halo virial mass in code units
 * @param rhalo Host halo virial radius in code units
 * @return Density at radius @a r in code units (>= 0) on success, -1 if
 * an error occurs
 */
float density_SIS(float r, float mhalo, float rhalo)
{
  if (r<=0) {
    fprintf(stderr, "Error (density_SIS): cannot calculate for "
	    "r = %g <= 0!\n", r); fflush(stderr);
    return -1;
  }

  if (mhalo<=0) {
    fprintf(stderr, "Error (density_SIS): mhalo = %g <= 0!\n", mhalo); 
    fflush(stderr);
    return -1;
  }

  if (rhalo<=0) {
    fprintf(stderr, "Error (density_SIS): rhalo = %g <= 0!\n", rhalo); 
    fflush(stderr);
    return -1;
  }

  return mhalo/(4*M_PI*r*r*rhalo);
}
