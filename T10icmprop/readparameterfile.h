extern char     path1[256];
extern char     path[256];
extern char     path2[256];
extern char     path3[256];
extern char     filebasename[256];
extern int      Files;
extern int 	MaxSnapshot;
extern int      StartSnapshot;
extern double 	Omega; 		
extern double   OmegaLambda;
extern double 	Hubble_h;
extern double	H0;
extern double 	G;  		
extern double   FracNeighbourSearchRadius;
extern float    BaryonFrac;
extern float    GasSearchRadiusGx0;
extern float    GasSearchRadius;
extern int      GasConstPartNum;
extern int      GasMinimum;

/* These two are parameters for the filtering of gas particles */
extern int      Niterations;       
extern float    Fracmax;

extern char	Identifier[256];
extern double	TimeBetSnapshot;
extern double	TimeOfFirstSnapshot;
extern double	UnitLength_in_cm;
extern double	UnitMass_in_g;
extern double	UnitVelocity_in_cm_per_s;

extern char     NameSnapshot[256];
extern char     NameSuborbits[256];
extern char     NameSPHordered[256];
extern char     NameICMproperties[256];
extern char     NameOutputsSelection[256];

/*** Code options */
extern int          Seed;
extern unsigned int ReorderOn;
extern unsigned int OutputListOn;
/*extern unsigned int ComovingToPhysicalOn;*/
