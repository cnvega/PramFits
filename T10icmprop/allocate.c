/**
 * @file allocate.c
 * @brief Functions to allocate memory to particle data structures
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void allocate_memory(void)
{
  int totbytes=0, bytes;
  
  printf("Allocating memory... \n"); fflush(stdout);
  
  /* 
   * Due to the reordering which includes gas particles which have
   * disappeared, I must reserve memory using the number of particles in
   * snapshot 0 (NumPart0) 
   */
  if (NumPart0 > 0) {
    
    if (!(P=malloc(bytes=NumPart0*sizeof(struct particle_data)))) {
      fprintf(stderr,"Error (allocate_memory): failed to allocate memory."
	      " (A)\n"); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    totbytes += bytes;
    
    P--;   /* Start with offset 1 */
    
    if (!(Id=malloc(bytes=NumPart0*sizeof(int)))) {
      fprintf(stderr,"Error (allocate_memory): failed to allocate memory."
	      " (B)\n"); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    totbytes += bytes;

    Id--;   /* Start with offset 1 */
  }

  printf("Allocated memory for P and Id: %f Mbyte.\n", 
	 ((float)totbytes)/(1024.0*1024.0)); fflush(stdout);

  return;
}


void free_memory(void)
{
  printf("Freeing space for particle positions P\n");
  P++;
  free(P);
  
  return;
}
