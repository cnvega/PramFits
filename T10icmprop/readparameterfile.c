/**
 * @file readparameterfile.c
 * @brief Read input parameters from file. Version for use 
 * with @a icmpropertyfinder.c
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

int  Seed;
char path1[256];
char path[256];
char path2[256];
char path3[256];
char filebasename[256];
int Files;
int StartSnapshot;
int MaxSnapshot;
double Omega;
double OmegaLambda;
double Hubble_h;	/* Hubble constant in units of 100 km/s/Mpc */
double H0;		/* Hubble constant in code units (see gadget output)*/
double G;  		/* and G to compute interparticle separation */
double FracNeighbourSearchRadius; /* Fraction of radius to distribute metals 
				     in gas particles around each gx */ 
float GasSearchRadiusGx0;          /* Para buscar particulas de gas alrededor
				      de gx tipo 0 que pueden ser la central 
				      del cumnulo */
float GasSearchRadius;          /* Para buscar particulas de gas alrededor
				   de las galaxias para calcular ram
				   pressure */
float BaryonFrac;
int GasConstPartNum;
int GasMinimum;

/* These two are parameters for the filtering of gas particles */
int Niterations;       
float Fracmax;

char 	Identifier[256];
double	TimeBetSnapshot;
double	TimeOfFirstSnapshot;
double	UnitLength_in_cm;
double	UnitMass_in_g;
double	UnitVelocity_in_cm_per_s;

char    NameSnapshot[256];
char    NameSuborbits[256];
char    NameSPHordered[256];
char    NameICMproperties[256];
char    NameOutputsSelection[256];

unsigned int OutputListOn;
unsigned int ReorderOn;
/*unsigned int ComovingToPhysicalOn;*/


/**
 * @brief Read input parameters from file
 * @param filename Name of input file
 */
void readparameterfile(char *filename)
{
  FILE *fp;
  char buffer[256], buffer2[256];
  char *error;
  int elements;
  
  
  fp = fopen(filename, "r");
  if (fp==NULL) {
    fprintf(stdout, "Error: cannot open file '%s': %s (%u)\n", 
	    filename, strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);     
  }
  
  fprintf(stdout,"Reading parameters from file '%s'...\n", filename);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Seed%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  Seed = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path1%s",path1);
  checkforerror(error,elements,&buffer[0]);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path%s",path);
  checkforerror(error,elements,&buffer[0]);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path2%s",path2);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path3%s",path3);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FileBasename%s",filebasename);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Files%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  Files = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"StartSnapshot%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  StartSnapshot = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MaxSnapshot%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  MaxSnapshot = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Omega%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  Omega = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);  
  elements = sscanf(buffer,"OmegaLambda%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  OmegaLambda = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Hubble_h%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  Hubble_h = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"H0%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  H0 = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"G%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  G = atof(buffer2);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BaryonFrac%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  BaryonFrac = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GasSearchRadiusGx0%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  GasSearchRadiusGx0 = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GasSearchRadius%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  GasSearchRadius = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GasConstPartNum%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  GasConstPartNum = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GasMinimum%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  GasMinimum = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Niterations%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  Niterations = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Fracmax%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  Fracmax = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Identifier%s",Identifier);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeBetSnapshot%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  TimeBetSnapshot = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeOfFirstSnapshot%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  TimeOfFirstSnapshot = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitLength_in_cm%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitLength_in_cm = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitMass_in_g%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitMass_in_g = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitVelocity_in_cm_per_s%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitVelocity_in_cm_per_s = atof(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSnapshot%s",NameSnapshot);
  checkforerror(error,elements,&buffer[0]);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSuborbits%s",NameSuborbits);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSPHordered%s",NameSPHordered);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameICMproperties%s",NameICMproperties);
  checkforerror(error,elements,&buffer[0]);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameOutputsSelection%s",NameOutputsSelection);
  checkforerror(error,elements,&buffer[0]);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"OutputListOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  OutputListOn = atoi(buffer2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ReorderOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  ReorderOn = atoi(buffer2);
  
  /* TOMAS 2012-08-08: I remove this because it only generates confusion.
     Data from snapshots is always in comoving coordinates, positions
     and velocities from SAG in comoving coordinates, other data in 
     physical */
  /*error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ComovingToPhysicalOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  ComovingToPhysicalOn = atoi(buffer2);*/

  return;
}


void checkforerror(char *error, int elements, char *buf)
{
  char buffer2[256]; 
  
  
  if (error==NULL) {
    fprintf(stderr, "Error: couldn't read from parameter file\n");
    fprintf(stderr, "Check input file and try again. Stop\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  if (elements==0) {
    strncpy(buffer2,buf,strlen(buf));
    /* strncpy does not add a final \0 terminator. It has to be added 
       by hand */
    buffer2[strlen(buf)-1]='\0';
    fprintf(stderr, "Error: couldn't convert line '%s'. Stop.\n", buffer2);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  return;
}
