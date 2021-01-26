#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void reorderinggas(void)
{
  int i, j, minid;
  float xyzsave[3], xyzsource[3];
  float usave, usource;
  float rhosave, rhosource;
  float nesave, nesource;
  float tempsave, tempsource;
  /*  float velsave[3], velsource[3]; */
  int idsource, idsave, dest;
  
  
  printf("Reordering gas particles...\n");

  for (i=1, minid=Id[i]; i<=NumPart; i++) {
    if (Id[i]<minid)
      minid = Id[i];
  }
  printf("Minimum ID = %d\n", minid);
  for (i=1; i<=NumPart; i++)
    Id[i]-= (minid-1);
  
  for (i=1; i<=NumPart; i++) {
    if (Id[i] != i) {
      for (j=0;j<3;j++) {
	xyzsource[j] = P[i].Pos[j];
	/* velsource[j]=P[i].Vel[j]; */ 
      }
      usource = PG[i].U;
      rhosource = PG[i].Rho;
      nesource = PG[i].Ne;
      tempsource = PG[i].Temp; 
      idsource = Id[i];
      dest = Id[i];
      
      do {
	for (j=0; j<3; j++) {
	  xyzsave[j] = P[dest].Pos[j];
	  /* velsave[j]=P[dest].Vel[j]; */
	}
	usave = PG[dest].U;
	rhosave = PG[dest].Rho;
	nesave = PG[dest].Ne;
	tempsave = PG[dest].Temp; 
	idsave = Id[dest];
	
	for (j=0; j<3; j++) {
	  P[dest].Pos[j]= xyzsource[j];
	  /* P[dest].Vel[j]= velsource[j]; */
	}
	PG[dest].U = usource;
	PG[dest].Rho = rhosource;
	PG[dest].Ne = nesource;
	PG[dest].Temp = tempsource;
	Id[dest] = idsource;
	
	if (dest == i) 
	  break;
	
	for (j=0; j<3; j++) {
	  xyzsource[j]= xyzsave[j];
	  /* velsource[j]= velsave[j]; */
	}
	usource = usave;
	rhosource = rhosave;
	nesource = nesave;
	tempsource = tempsave;
	idsource = idsave;
	dest = idsource;
      }
      while(1);
    }
  }
  
  printf("done.\n");
  
  Id++;
  free(Id);
  printf("Space for particle ID freed.\n");
  
  return;
}
