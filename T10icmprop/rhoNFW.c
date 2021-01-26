/**
 * @file rhoNFW.c
 * @brief Function to calculate density (using a NFW profile)
 * @author Tomas E. Tecce
 * @date 2013-02-12
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

/**
 * @brief Calculate density as a function of radius using a NFW profile
 * @param r Radius in code units
 * @param mhalo Host halo virial mass in code units
 * @param rhalo Host halo virial radius in code units
 * @param chalo Host halo concentration
 * @param zcurr Redshift
 * @return Density at radius @a r in code units (>= 0) on success, -1 if
 * an error occurs
 */
float density_NFW(float r, float mhalo, float rhalo, float chalo, 
		  float zcurr)
{
  float cx, delta0;
  double hubble_of_z, rhocrit;

  if (r<=0) {
    fprintf(stderr, "Error (density_NFW): cannot calculate for "
	    "r = %g <= 0!\n", r); fflush(stderr);
    return -1;
  }

  if (mhalo<=0) {
    fprintf(stderr, "Error (density_NFW): mhalo = %g <= 0!\n", mhalo); 
    fflush(stderr);
    return -1;
  }

  if (rhalo<=0) {
    fprintf(stderr, "Error (density_NFW): rhalo = %g <= 0!\n", rhalo); 
    fflush(stderr);
    return -1;
  }

  if (chalo<=0 || chalo>50) {
    fprintf(stderr, "Error (density_NFW): chalo = %g outside of allowed "
	    "range (0,50]!\n", chalo); fflush(stderr);
    return -1;
  }

  /* Critical density at the current redshift */
  /* This is H^2(z) */
  hubble_of_z = Hubble*Hubble*(Omega * pow(1+zcurr,3) + 
			       (1-Omega-OmegaLambda)*pow(1+zcurr,2) + 
			       OmegaLambda);
  rhocrit = 3*hubble_of_z/(8*M_PI*G);

  /* Parameter delta0 from the NFW */
  delta0 = (200.0/3.0)*chalo*chalo*chalo/(log(1+chalo)-(chalo/(1+chalo)));
  
  cx = chalo*r/rhalo;

  return (float)rhocrit*delta0/(cx*(1+cx)*(1+cx));
}
