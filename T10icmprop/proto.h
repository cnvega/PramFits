/**
 * @file proto.h
 * @brief Function prototypes for use with @a icmpropertyfinder.c
 * @author Tomas E. Tecce
 */
#ifdef HDF5OUTPUT
#include <hdf5.h>
#include <hdf5_hl.h>
#endif


void   allocate_gasmemory();
void   allocate_memory();
void   checkforerror(char *, int, char *);
void   checkwrite(int, int);
float  density_NFW(float, float, float, float, float);
float  density_SIS(float, float, float);
void   dump_icm(char *, int);
void   find_gaspart_withtree(int);
int    force_treebuild();
void   force_update_node_recursive(int, int, double, double, double, double);
void   free_memory();
void   free_gasmemory();
void   get_galaxies_data(int, char *);
void   get_galaxies_data_HDF5(int, char *);
void   get_ICM_properties(int, int, int);
int    getNumPart(char *, int, int);
int    indexrm(int, int, int, int);
void   init();
double integrand(double);
void   loadpositions(char *, int, int);
float  ngb_select_closest(int, int, float *, int *);
void   ngb_treeallocate(int);
float  ngb_treefind(float [3], int, float);
int    ngb_treefind_variable(float [3], double);
void   ngb_treefree();
void   polint(double [], double [], int, double, double *, double *);
void   readparameterfile(char *);
void   reorderinggas();
void   set_units();
void   sort2_flt_int(unsigned long, float [], int []);
double time_to_present(double);

#ifdef HDF5OUTPUT
void   checkwrite_HDF5(herr_t);
herr_t dumpHDF5_file_attr(hid_t, int, float);
void   dumpHDF5_icm(char *, int, float);
int    set_float_variable_attr(hid_t, char *, char *, double *, char *);
#endif

#ifdef ORBITS
void read_suborbits_A(char *, int);
#endif
