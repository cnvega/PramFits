/**
 * @file io_input_G2.c
 * @brief I/O functions to read data from Gadget-2 format files, for use 
 * with @a icmpropfinder.c
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda0;
  double HubbleParam0;
  char fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; /* Fills to 256 bytes */
} header1;


/*  LT code */
struct block_header_data
{
  int first_tag;
  char block_name[4];
  int block_length, second_tag;
} block_header;

#define READ_BLOCK_HEADER fread(&block_header, sizeof(block_header), 1, fd);
/* :) :) :) */

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);


int getNumPart(char *fname, int files, int type)
{
  FILE *fd;
  char buf[200];
  int i,dummy;
  

  i = 0;
  
  if (files>1)
    sprintf(buf,"%s.%d",fname,i);
  else
    sprintf(buf,"%s",fname);
  
  if (!(fd = fopen(buf, "r"))) {
    printf("Error (getNumPart): can't open file `%s`: %s (%u)\n",
	   buf, strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  
  printf("Reading `%s' ...\n", buf); fflush(stdout);
  
  READ_BLOCK_HEADER;
  SKIP;
  fread(&header1, sizeof(header1), 1, fd);
  SKIP;

  return header1.npart[type];
}


/**
 * @brief Load positions for simulation particles
 * @param fname Simulation snapshot file
 * @param files Number of files per snapshot
 * @param type Particle type as in Gadget
 */
void loadpositions(char *fname, int files, int type)
{
  FILE *fd;
  char buf[200];
  int i, k, dummy;
  int n, pc, pc_new, id, check;
  float pos[3];
  float vel[3];
  float masa;
  int ntot_withmasses;


  for (i=0, pc=1; i<files; i++, pc=pc_new) {
    if (files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    if (!(fd=fopen(buf,"r"))) {
      printf("Error (loadpositions): can't open file '%s': %s (%u)\n", 
	     buf, strerror(errno), errno); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    
    printf("Reading '%s' ...\n", buf); fflush(stdout);
    
    READ_BLOCK_HEADER;
    printf("READ_BLOCK_HEADER HEADER \n");
    printf("  first_tag: %d \n", block_header.first_tag);
    printf("  block_name: %s \n", block_header.block_name);
    printf("  block_length: %d \n", block_header.block_length);
    printf("  second_tag: %d \n", block_header.second_tag);
    
    SKIP;
    fread(&header1, sizeof(header1), 1, fd);
    SKIP;
    
    if (files==1)
      NumPart=header1.npart[type];
    else
      NumPart=header1.npartTotal[type];
    
    printf("NumPart: %d \n",NumPart); fflush(stdout);

    for (k=0; k<6; k++)
      printf("type = %d npart = %d mass = %g\n", k, 
	     header1.npart[k], header1.mass[k]);
    printf("header1.time = %g\n", header1.time);
    
    for (k=0, ntot_withmasses=0; k<5; k++) {
      if (header1.mass[k]==0)
	ntot_withmasses+= header1.npart[k];
    }
    
    printf("ntot_withmasses:%d \n",ntot_withmasses);
    fflush(stdout);
    
    printf("header1.redshift = %g\n", header1.redshift);
    printf("header1.flag_sfr = %d\n", header1.flag_sfr);
    printf("header1.flag_feedback = %d\n", header1.flag_sfr);
    printf("header1.npartTotal[0] = %d\n", header1.npartTotal[0]);
    printf("header1.npartTotal[1] = %d\n", header1.npartTotal[1]);
    printf("header1.npartTotal[2] = %d\n", header1.npartTotal[2]);
    printf("header1.npartTotal[3] = %d\n", header1.npartTotal[3]);
    printf("header1.npartTotal[4] = %d\n", header1.npartTotal[4]);
    printf("header1.npartTotal[5] = %d\n", header1.npartTotal[5]);
    printf("header1.flag_cooling = %d\n", header1.flag_cooling);
    printf("header1.num_files = %d\n", header1.num_files);
    printf("header1.BoxSize = %g\n", header1.BoxSize);
    printf("header1.Omega0 = %g\n", header1.Omega0);
    printf("header1.OmegaLambda0 = %g\n", header1.OmegaLambda0);
    printf("header1.HubbleParam0 = %g\n", header1.HubbleParam0);
    
    if (i==0) allocate_memory();
    
    READ_BLOCK_HEADER;
    printf("READ_BLOCK_HEADER POS \n");
    printf("  first_tag: %d \n", block_header.first_tag);
    printf("  block_name: %s \n", block_header.block_name);
    printf("  block_length: %d \n", block_header.block_length);
    printf("  second_tag: %d \n", block_header.second_tag);
    
    SKIP;
    for (k=0, pc_new=pc; k<6; k++) {
      for (n=0; n<header1.npart[k]; n++) {
	check = fread(&pos[0],sizeof(float),3,fd);
	if (check !=3) {
	  printf("Error (loadpositions) when reading positions. Stop\n");
	  fflush(stdout);
	  exit(EXIT_FAILURE);
	}
	
	if (k==type) {
	  P[pc_new].Pos[0] = pos[0];
	  P[pc_new].Pos[1] = pos[1];
	  P[pc_new].Pos[2] = pos[2];
	  pc_new++;
	}
      }
    }
    printf("\n\nP[1].Pos[0]: %g \n", P[1].Pos[0]);
    printf("P[2].Pos[0]: %g \n", P[2].Pos[0]);
    printf("P[NumPart-1].Pos[0]: %g \n", P[NumPart-1].Pos[0]);
    printf("P[NumPart].Pos[0]: %g \n", P[NumPart].Pos[0]);
    SKIP;
    
    READ_BLOCK_HEADER;
    SKIP;
    
    for (k=0, pc_new=pc; k<6; k++) {
      for (n=0; n<header1.npart[k]; n++) {
	check = fread(&vel[0],sizeof(float),3,fd);
	if (check !=3) {
	  printf("Error (loadpositions) when reading velocities. Stop\n");
	  fflush(stdout);
	  exit(EXIT_FAILURE);
	}
	
	if (k==type) {
	  P[pc_new].Vel[0] = vel[0];
	  P[pc_new].Vel[1] = vel[1];
	  P[pc_new].Vel[2] = vel[2];
	  pc_new++;
	}
      }
    }
    SKIP;
      
    printf("\n\nP[1].Vel[0]:%g \n",P[1].Vel[0]);
    printf("P[2].Vel[0]:%g \n",P[2].Vel[0]);
    printf("NumPart:%d \n\n",NumPart);
    printf("P[NumPart-1].Vel[0]:%g \n",P[NumPart-1].Vel[0]);
    printf("P[NumPart].Vel[0]:%g \n",P[NumPart].Vel[0]);
    
    READ_BLOCK_HEADER;
    SKIP;
    for (k=0, pc_new=pc; k<6; k++) {
      for(n=0; n<header1.npart[k]; n++) {
	check = fread(&id,sizeof(float),1,fd); /*decia sizeof(float)!!*/
	if (check != 1) {
	  printf("Error (loadpositions) when reading identities. Stop\n");
	  fflush(stdout);
	  exit(EXIT_FAILURE);
	}
	
	if (k == type) {
	  Id[pc_new] = id;
	  pc_new++;
	}
      }
    }
    printf("\n\nId[1]: %d \n", Id[1]);   
    printf("Id[2]: %d \n", Id[2]);  
    printf("NumPart: %d \n\n", NumPart); 
    printf("Id[NumPart-1]: %d \n", Id[NumPart-1]);
    printf("Id[NumPart]: %d \n", Id[NumPart]);
    SKIP;
      
    if (ntot_withmasses > 0) { 
      READ_BLOCK_HEADER;
      SKIP;
      for (k=0, pc_new=pc; k<6; k++) {
	if (header1.mass[k]==0) {
	  printf("k: %d, header1.mass[k]: %g \n", k, header1.mass[k]);

	  for(n=0; n<header1.npart[k]; n++) {
	    check = fread(&masa,sizeof(float),1,fd); /*decia sizeof(float)!!*/
	    if (check != 1) {
	      printf("Error (loadpositions) when reading masses. Stop\n");
	      fflush(stdout);
	      exit(EXIT_FAILURE);
	    }
	    
	    if (k==type) {
	      P[pc_new].Mass = masa;
	      pc_new++;
	    }
	  }
	}
	else {
	  P[pc_new].Mass = header1.mass[k];
	  pc_new++;
	} 
      }
      
      printf("\n\nP[1].Mass:%g \n", P[1].Mass);   
      printf("P[2],Mass:%g \n", P[2].Mass);   
      printf("P[NumPart-1].Mass:%g \n", P[NumPart-1].Mass);
      printf("P[NumPart].Mass:%g \n", P[NumPart].Mass);
      SKIP;
    } 
    
    if (type==0) /* Variables loaded only for gas particles */
      if (header1.npart[type]>0) {
	if (i==0)
	  allocate_gasmemory();

	SKIP;
	for (n=0, pc_new=pc; n<header1.npart[type];n++) {
	  fread(&PG[pc_new].U, sizeof(float), 1, fd);
	  pc_new++;
	}
	SKIP;
	
	SKIP;
	for (n=0, pc_new=pc; n<header1.npart[type];n++) {
	  fread(&PG[pc_new].Rho, sizeof(float), 1, fd);
	  pc_new++;
	}
	
	SKIP;
	if (header1.flag_cooling) {
	  SKIP;
	  for (n=0, pc_new=pc; n<header1.npart[type];n++) {
	    fread(&PG[pc_new].Ne, sizeof(float), 1, fd);
	    pc_new++;
	  }
	  SKIP;
	}
	else
	  for (n=0, pc_new=pc; n<header1.npart[type];n++) {
	    PG[pc_new].Ne = 1.0;
	    pc_new++;
	  }
      }
    
    fclose(fd);
  }
    
  Time = header1.time;
  
  PartMassType = P[1].Mass; /* Valid for DM and gas particles at 
			       snapshot=0 since they all 
			       have equal masses */
  printf("type: %d, header1.npart[type]: %d, PartMassType:%g \n",
	 type, header1.npart[type], PartMassType);
  
#ifdef PERIODIC
  BoxSize = header1.BoxSize;
#endif
  
  return;
}
