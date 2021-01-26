/**
 * @file allocategas.c
 * @brief Functions to allocate memory to gas particle data structures
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void allocate_gasmemory(void)
{
  long long bytes;
  long long totbytesPG=0;

  /* 
   * Due to the reordering which includes gas particles which have
   * disappeared, I must reserve memory using the number of particles in
   * snapshot 0 (NumPart0) 
   */
  if (NumPart0 > 0) {
    if (!(PG=malloc(bytes=NumPart0*sizeof(struct gaspart_data)))) {
      printf("Error (allocate_gasmemory): failed to allocate gas memory: "
	     "PG. (C)\n"); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    totbytesPG += bytes;
    
    PG--;   /* Start with offset 1 */
  }
  
  printf("Allocated memory for gas particles PG: %f Mbyte.\n", 
	 ((float)totbytesPG)/(1024.0*1024.0)); fflush(stdout);
  
  return;
}


void free_gasmemory(void)
{
  printf("Freeing space for gas particle PG\n");
  PG++;
  free(PG);
  
  return;
}
