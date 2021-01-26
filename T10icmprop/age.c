/**
 * @file age.c
 * @brief Functions to calculate time elapsed to present
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

/* Return time to present as a function of redshift */
double time_to_present(double z)  
{
  double time;
  double integrand(double a);
  double qromb(double (*func)(double), double a, double b);

  time = 1.0/Hubble*qromb(integrand, 1/(z+1), 1);

  return time;
}


double integrand(double a)
{
  return 1.0/sqrt(Omega/a + (1-Omega-OmegaLambda) + OmegaLambda*a*a);
}









