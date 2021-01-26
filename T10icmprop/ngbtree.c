/**
 * @file ngbtree.c
 * @brief Tree search COMPLETAR
 */
#include "allvars.h"
#include "proto.h"
#include "readparameterfile.h"

/* 
   Index convention for accessing tree nodes:
   the indices 0...NumPart0 reference single particles
   the indices All.MaxPart... All.MaxPart+nodes-1 reference tree nodes
 
   `Nodes_base' points to the first tree node
   `nodes' is shifted, such that nodes[All.MaxPart] gives the first 
   tree node
*/

static int MaxNodes;

static int last; /* Auxiliary variable used to set-up non-recursive walk */

double part_dens;

static struct particle_data *PP;


int force_treebuild(void)
{
  int i, j, subnode, numnodes;
  int nfree, th, nn, startnode;
  double xmin[3], xmax[3], len, center[3];
  double Center[3], Len;
  double lenhalf;
  struct NODE *nfreep;
 
  PP = P+1;
  
  startnode = 0; 
  numnodes  = 0;
  Nodes     = Nodes_base - NumPart0;

  /* Find enclosing rectangle */
  for (j = 0; j < 3; j++)	
    xmin[j] = xmax[j] = PP[0].Pos[j];
  
  for (i = 1; i < NumPart0; i++)
    for (j = 0; j < 3; j++) {
      if (PP[i].Pos[j] > xmax[j])
	xmax[j] = PP[i].Pos[j];
      if (PP[i].Pos[j] < xmin[j])
	xmin[j] = PP[i].Pos[j];
    }

  part_dens = NumPart0/((xmax[0]-xmin[0]) * (xmax[1]-xmin[1]) * (xmax[2]-xmin[2]));

  /* Determine maxmimum extension */
  for (j = 1, len = xmax[0] - xmin[0]; j < 3; j++)
    if ((xmax[j] - xmin[j]) > len)
      len = xmax[j] - xmin[j];
  
  len *= 1.00001;

  for (j = 0; j < 3; j++)
    Center[j] = 0.5 * (xmax[j] + xmin[j]);
  Len = len;
  
  /* Create a root node and insert first particle as its leaf */
  nfree = NumPart0 + startnode;	/* index */
  nfreep = &Nodes[nfree];	/* select first node */
  
  for (i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;
  
  numnodes++;
  nfree++;
  nfreep++;
  
  for (i = 0; i < NumPart0; i++)	{ /* Insert all other particles */
    th = NumPart0 + startnode;	/* select index of first node in tree */
    
    len = Len;
    lenhalf = Len/2;
    for(j = 0; j < 3; j++)
      center[j] = Center[j];
    
    while(1) {	  /* "th" will always point to an internal node */
      len     *= 0.5;
      lenhalf *= 0.5;
      
      subnode = 0;
      if (PP[i].Pos[0] > center[0]) {
	center[0] += lenhalf;
	subnode += 1;
      }
      else {
	center[0] -= lenhalf;
      }
      if (PP[i].Pos[1] > center[1]) {
	center[1] += lenhalf;
	subnode += 2;
      }
      else {
	center[1] -= lenhalf;
      }
      if (PP[i].Pos[2] > center[2]) {
	center[2] += lenhalf;
	subnode += 4;
      }
      else {
	center[2] -= lenhalf;
      }
      
      nn = Nodes[th].u.suns[subnode];
      
      if (nn >= 0) {	/* Ok, something is in the daughter slot already, 
			   need to continue */
	if (nn >= NumPart0) { /* The daughter node is an internal node. */
	  th= nn;
	}
	else  {  /* The daughter node is a particle. 
		    Need to generate new node */
	  
	  /* "nn"  is the particle index */
          
	  Nodes[th].u.suns[subnode] = nfree;
	  nfreep->u.suns[0] = -1;
	  nfreep->u.suns[1] = -1;
	  nfreep->u.suns[2] = -1;
	  nfreep->u.suns[3] = -1;
	  nfreep->u.suns[4] = -1;
	  nfreep->u.suns[5] = -1;
	  nfreep->u.suns[6] = -1;
	  nfreep->u.suns[7] = -1;
	  
	  subnode = 0;
	  if (PP[nn].Pos[0] > center[0])
	    subnode += 1;
	  if (PP[nn].Pos[1] > center[1])
	    subnode += 2;
	  if (PP[nn].Pos[2] > center[2])
	    subnode += 4;
	  
	  nfreep->u.suns[subnode] = nn;
	  
	  th = nfree;	/* Resume trying to insert the new particle at 
			   the newly created internal node */
	  
	  numnodes++;
	  nfree++;
	  nfreep++;
	  
	  if ((numnodes + startnode) >= MaxNodes) {
	    printf("Error: Maximum number %d of tree-nodes reached.\n", 
		   MaxNodes);
	    printf("For particle %d\n", i);
	    exit(EXIT_FAILURE);
	  }
	}
      }
      else {
	/* 
	 * Here we have found an empty slot where we can 
	 * attach the new particle as a leaf 
	 */
	Nodes[th].u.suns[subnode] = i;
	break;	/* done for this particle */
      }
    }
  }
  
  /* Now walk the tree ones to set-up centers and faster walk */
  
  last = -1;
  
  force_update_node_recursive(NumPart0 + startnode, -1, Len, Center[0], 
			      Center[1], Center[2]);
  
  if (last >= 0) {
    if (last<NumPart0)
      PP[last].Nextnode = -1;
    else
      Nodes[last].u.d.nextnode =-1;
  }
  
  printf("Used %d nodes out of allocated %d. (filled fraction %g)\n",
         numnodes, MaxNodes, (double)numnodes/MaxNodes);
  
  return numnodes;
}


void force_update_node_recursive(int no, int sib, double len, double cx, 
				 double cy, double cz)
{
  int j, jj, p, pp = 0, nextsib, suns[8];
  double ccx, ccy, ccz;
  
  if (no >= NumPart0)  { /* Internal node */
    for (j = 0; j < 8; j++)
      suns[j] = Nodes[no].u.suns[j];  /* This "backup" is necessary because 
					 the nextnode entry will
					 overwrite one element (union!) */
    Nodes[no].u.d.center[0] = cx;
    Nodes[no].u.d.center[1] = cy;
    Nodes[no].u.d.center[2] = cz;
    Nodes[no].u.d.len = len;
    Nodes[no].u.d.sibling = sib;
    
    if (last >= 0) {
      if (last >= NumPart0)
	Nodes[last].u.d.nextnode = no;
      else
	PP[last].Nextnode = no;
    }
    last = no;
    
    for (j = 0; j < 8; j++) {
      if ((p = suns[j]) >= 0) {
	ccx = cx;
	ccy = cy;
	ccz = cz;
	
	if ((j & 1)) 
	  ccx += len/4;
	else
	  ccx -= len/4;
	if ((j & 2)) 
	  ccy += len/4;
	else
	  ccy -= len/4;
	if ((j & 4)) 
	  ccz += len/4;
	else
	  ccz -= len/4;
	
	/* Check if we have a sibling on the same level */
	for (jj = j + 1; jj < 8; jj++)
	  if ((pp = suns[jj]) >= 0)
	    break;
	
	if (jj < 8)	/* Yes, we do */
	  nextsib = pp;
	else
	  nextsib = sib;
	
	force_update_node_recursive(p, nextsib, len/2, ccx, ccy, ccz);
      }
    }
  }
  else {    /* Single particle */
    if (last >= 0) {
      if (last >= NumPart0)
	Nodes[last].u.d.nextnode = no;
      else
	PP[last].Nextnode = no;
    }
    last = no;
  }

  return;
}


float ngb_treefind(float xyz[3], int desngb, float hguess)
{
  int   numngb;
  float h2max;


  if (hguess==0)
    hguess = pow(3*desngb/(4*M_PI)/(part_dens), 0.33);

  do {
    numngb = ngb_treefind_variable(xyz, hguess);
    
    if (numngb == MAX_NGB) {
      hguess /= 1.1;
      continue;
    }
    
    if (numngb < desngb) {
      hguess *= 1.26;
      continue;
    }
    
    if (numngb>=desngb) {
      h2max = ngb_select_closest(desngb, numngb, R2list-1, IndexList-1);
      break;
    }
    
    hguess *=1.26; 
  }
  while(1);
  
  return h2max;
}


float ngb_select_closest(int k, int n, float * arr, int *ind)
{
#define SWAP(a,b)  temp =(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;
  int i, ir, j, l, mid, ai, tempi;
  float a, temp;
  
  
  l = 1;
  ir = n;
  while (1) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
	SWAP(arr[l], arr[ir]);
	SWAPI(ind[l], ind[ir]);
      }
      return arr[k];
    }
    else {
      mid = (l + ir) >> 1;
      SWAP(arr[mid], arr[l + 1]);
      SWAPI(ind[mid], ind[l + 1]);
      
      if (arr[l] > arr[ir]) {
	SWAP(arr[l], arr[ir]);
	SWAPI(ind[l], ind[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
	SWAP(arr[l + 1], arr[ir]);
	SWAPI(ind[l + 1], ind[ir]);
      }
      if (arr[l] > arr[l + 1]) {
	SWAP(arr[l], arr[l + 1]);
	SWAPI(ind[l], ind[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ai = ind[l + 1];
      
      while (1) {
	do
	  i++;
	while (arr[i] < a);
	
	do
	  j--;
	while (arr[j] > a);
	
	if (j < i)
	  break;
	SWAP(arr[i], arr[j]);
	SWAPI(ind[i], ind[j]);
      }
      
      arr[l + 1] = arr[j];
      arr[j] = a;
      ind[l + 1] = ind[j];
      ind[j] = ai;
      
      if (j >= k)
	ir = j - 1;
      if (j <= k)
	l = i;
    }
  }
#undef SWAP
#undef SWAPI
}


/* 
 * These macros map a coordinate difference
 * to the nearest periodic image
 */
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))


/**
 * @brief This function returns all neighbours (and only those) with 
 * distance <= @a hguess and returns them in @a ngblistback and 
 * @a r2listback
 */
int ngb_treefind_variable(float searchcenter[3], double hguess)
{
  int k, numngb;
  int no, p, no_save;
  double dx, dy, dz, r2, h2;
  struct NODE *this;
  double searchmin[3], searchmax[3];
  
  
  h2 = hguess*hguess;
  
  for (k=0; k<3; k++) {    /* Cube-box window */
    searchmin[k] = searchcenter[k] - hguess;
    searchmax[k] = searchcenter[k] + hguess;
  }

  numngb = 0;
  no = NumPart0;
  
  while (no >= 0) {
    if (no < NumPart0) { /* Single particle */
      p = no;
      no = PP[no].Nextnode;
      
#ifdef PERIODIC
      if (NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]) < (searchmin[0] - searchcenter[0]))
	continue;
      if (NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]) > (searchmax[0] - searchcenter[0]))
	continue;
      if (NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]) < (searchmin[1] - searchcenter[1]))
	continue;
      if (NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]) > (searchmax[1] - searchcenter[1]))
	continue;
      if (NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]) < (searchmin[2] - searchcenter[2]))
	continue;
      if (NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]) > (searchmax[2] - searchcenter[2]))
	continue;
#else
      if (PP[p].Pos[0] < (searchmin[0]))
	continue;
      if (PP[p].Pos[0] > (searchmax[0]))
	continue;
      if (PP[p].Pos[1] < (searchmin[1]))
	continue;
      if (PP[p].Pos[1] > (searchmax[1]))
	continue;
      if (PP[p].Pos[2] < (searchmin[2]))
	continue;
      if (PP[p].Pos[2] > (searchmax[2]))
	continue;
#endif

#ifdef PERIODIC
      dx = NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]);
      dy = NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]);
      dz = NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]);
#else
      dx = PP[p].Pos[0] -  searchcenter[0];
      dy = PP[p].Pos[1] -  searchcenter[1];
      dz = PP[p].Pos[2] -  searchcenter[2];
#endif
      r2=  dx*dx + dy*dy +dz*dz;
      
      if (r2 < h2) {
	if (numngb < MAX_NGB) {
	  IndexList[numngb] = p+1;
	  R2list[numngb++] = r2;
	}
	else
	  return numngb;
      }
    }
    else {
      this = &Nodes[no];
      
      no_save = no;
      no = Nodes[no].u.d.sibling;   /* In case the node can be discarded */
#ifdef PERIODIC
      if ((NGB_PERIODIC(this->u.d.center[0] - searchcenter[0]) + 0.5 * this->u.d.len) < (searchmin[0] - searchcenter[0]))
	continue;
      if ((NGB_PERIODIC(this->u.d.center[0] - searchcenter[0]) - 0.5 * this->u.d.len) > (searchmax[0] - searchcenter[0]))
	continue;
      if ((NGB_PERIODIC(this->u.d.center[1] - searchcenter[1]) + 0.5 * this->u.d.len) < (searchmin[1] - searchcenter[1]))
	continue;
      if ((NGB_PERIODIC(this->u.d.center[1] - searchcenter[1]) - 0.5 * this->u.d.len) > (searchmax[1] - searchcenter[1]))
	continue;
      if ((NGB_PERIODIC(this->u.d.center[2] - searchcenter[2]) + 0.5 * this->u.d.len) < (searchmin[2] - searchcenter[2]))
	continue;
      if ((NGB_PERIODIC(this->u.d.center[2] - searchcenter[2]) - 0.5 * this->u.d.len) > (searchmax[2] - searchcenter[2]))
	continue;
#else
      if ((this->u.d.center[0] + 0.5 * this->u.d.len) < (searchmin[0]))
	continue;
      if ((this->u.d.center[0] - 0.5 * this->u.d.len) > (searchmax[0]))
	continue;
      if ((this->u.d.center[1] + 0.5 * this->u.d.len) < (searchmin[1]))
	continue;
      if ((this->u.d.center[1] - 0.5 * this->u.d.len) > (searchmax[1]))
	continue;
      if ((this->u.d.center[2] + 0.5 * this->u.d.len) < (searchmin[2]))
	continue;
      if ((this->u.d.center[2] - 0.5 * this->u.d.len) > (searchmax[2]))
	continue;
#endif
      no = this->u.d.nextnode;	/* ok, we need to open the node */
    }
  }
  
  return numngb;
}


/**
 * @brief Allocate memory for the neighbour lists
 */
void ngb_treeallocate(int maxnodes)
{
  int totbytes=0, bytes;


  MaxNodes = maxnodes;
  
  if (!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE)))) {
    printf("Error (ngb_treeallocate): failed to allocate memory for %d "
	   "tree-nodes (%d bytes):\n "
	   "%s (%u)\n", MaxNodes, bytes, strerror(errno), errno);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  totbytes += bytes;

  if (!(R2list= malloc(MAX_NGB*sizeof(float)))) {
    printf("Error (ngb_treeallocate): failed to allocate memory for R2list:"
	   " %s (%u)\n", strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  totbytes += MAX_NGB*sizeof(float);
  
  printf("Allocated %f Mbyte for ngb search.\n", 
	 ((double) totbytes) / (1024.0 * 1024.0));
  
  if (!(IndexList= malloc(MAX_NGB*sizeof(float)))) {
    printf("Error (ngb_treeallocate): failed to allocate memory for IndexList:"
	   " %s (%u)\n", strerror(errno), errno); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  totbytes += MAX_NGB*sizeof(float);
  
  printf("Allocated %f Mbyte for ngb search.\n", 
	 ((double) totbytes) / (1024.0 * 1024.0));

  return;
}


/**
 * @brief Free memory allocated for tree building
 */
void ngb_treefree(void)
{
  free(IndexList);
  free(R2list);
  free(Nodes_base);

  return;
}
