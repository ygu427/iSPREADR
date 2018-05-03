/* C implementation of N-distance functions */
#include <R_ext/Print.h>          /* Rprintf */
#include <R.h>          /* for NA_REAL, includes math.h */

/* This function computes N-distance from the pairwise distance
   matrix. n: number of all observations, which is also the dimension
   of distmat; n1: number of obsrvations in group 1; ndist: the result. */
void ndist_from_distmat (double *distmat, int *n, int *n1, double *ndist)
{
  double WGX=0, WGY=0, BG=0;
  int i,j;
  int nn=*n, nx=*n1, ny=nn-nx;
  for (i=1; i<nx; i++)
    for (j=0; j<i; j++)
      WGX += distmat[i*nn+j];

  for (i=nx+1; i<nn; i++)
    for (j=nx; j<i; j++)
      WGY += distmat[i*nn+j];

  for (i=nx; i<nn; i++)
    for (j=0; j<nx; j++)
      BG += distmat[i*nn+j];

  WGX /=(nx*nx*0.5); WGY /=(ny*ny*0.5); BG /=(nx*ny*0.5);
  *ndist = sqrt(BG -WGX -WGY);
}

/* This function is the permutation version of
   ndist_from_distmat. combs: an integer matrix of
   combinations. combs[k] represents the k-th permutation and
   combs[k][i] is the index of the new i-th observation. By C
   convention, it starts with integer 0 instead of 1. Integer K is the
   number of combinations. ndistvec is a length K vector.*/
void ndist_from_distmat_perm (double *distmat, int *combs, int *K, int *n, int *n1, double *ndistvec)
{
  double WGX, WGY, BG;
  int i,j,k,ik,jk;              /* ik, jk: permuted indices */
  int nn=*n, nx=*n1, ny=nn-nx, kk=*K;
  /* compute the total distance. This is a constant for all permutations. */
  double TD=0;
  for (i=1; i<nn; i++)
    for (j=0; j<i; j++)
      TD += distmat[i*nn+j];

  for (k=0; k<kk; k++){
    WGX=0; WGY=0; BG=0;
    for (i=1; i<nx; i++){
      for (j=0; j<i; j++){
        ik=combs[k*nn+i]; jk=combs[k*nn+j];
        WGX += distmat[ik*nn+jk];
      }
    }
    for (i=nx+1; i<nn; i++){
      for (j=nx; j<i; j++){
        ik=combs[k*nn+i]; jk=combs[k*nn+j];
        WGY += distmat[ik*nn+jk];
      }
    }
    BG = TD -WGX -WGY;
    WGX /=(nx*nx*0.5); WGY /=(ny*ny*0.5); BG /=(nx*ny*0.5);
    ndistvec[k] = sqrt(BG -WGX -WGY);
  }
}

