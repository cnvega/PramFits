#include "allvars.h"
#include "proto.h"

/*
 * Returns the row-major offset of element (i,j) of a matrix of dimensions
 * sizex X sizey. Element indices are assumed to run from 1 to size. i is
 * the row index, j the column index.
 * To get the offset for indices running from 0 to size-1, add sizey to
 * the result from this function
*/
int indexrm(int i, int j, int sizex, int sizey)
{
  if (i>sizex || j>sizey) {
    fprintf(stderr, "Error (indexrm): index value too large - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
    if (i<0 || j<0) {
    fprintf(stderr, "Error (indexrm): negative index value - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return j + (sizey*(i-1));
} 
