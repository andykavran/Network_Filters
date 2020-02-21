/* Program to perform EM/BP community detection using the
 * degree-corrected SBM on an arbitrary network read from a file
 *
 * Written by Mark Newman  28 NOV 2014
 * Simplified the BP by including only leading-order terms  2 DEC 2014
 * Modified to include continuous metadata read in from a different
 *   file  15 DEC 2014
 * Modified to fit the priors to the metadata  16 DEC 2014
 */

/* Program control */

#define DEBUG
#undef VERBOSE

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "readwg2.h"

/* Constants */

#define N 5            // Number of groups
#define M 5            // Degree of Bernstein polynomials in expansion
#define LABELLEN 100
#define CHECKINT 1000
#define EPSILON 1e-10
#define MMIN 0.0
#define MMAX 1.0
#define GMMA_ACC 1e-3
#define BP_ACC 1e-3
#define BP_MAXSTEPS 20
#define EM_ACC 1e-3
#define EM_MAXSTEP 50
/* Types */

/* Globals */

NETWORK G;       // Struct storing the network
int twom;        // Twice the number of edges

double *meta;    // Metadata array

double gmma[N][M+1];               // Prior parameters, one set for each group
double omega[N][N];                // Parameters

double ***eta;    // Messages
double **q;       // One-point marginals

int binom[M+1];  // Binomial coefficients
int fail = 0;
gsl_rng *rng;


/* Read metadata from separate file */



void get_metadata(FILE *stream)
{
  int u;
  double m;
  char label[LABELLEN];

  meta = malloc(G.nvertices*sizeof(double));
  for (u=0; u<G.nvertices; u++) {
    fscanf(stream,"%s %lf",label,&m);
    if ((m<MMIN)||(m>MMAX)) meta[u] = gsl_rng_uniform(rng);
    else meta[u] = (m-MMIN)/(MMAX-MMIN);   // Convert to [0,1]
  }
}


/* Function to calculate binomial coefficients */

void binoms()
{
  int b,k;

  binom[0] = 1;
  b = M;
  for (k=1; k<M; k++) {
    b /= k;
    binom[k] = b;
    b *= M - k;
  }
  binom[M] = b/M;
}


/* Function to calculate priors for group r, metadata x */

double prior(int r, double x)
{
  int k;
  double result=0.0;
  double factor,ratio;

  if (x<EPSILON) return gmma[r][M];     // Avoid divide-by-zero
  factor = pow(x,M);
  ratio = (1-x)/x;
  for (k=0; k<=M; k++) {
    result += gmma[r][k]*binom[k]*factor;
    factor *= ratio;
  }
  return result;
}


/* Do BP */



int bp()
{
  int i,j;
  int u,v;
  int r,s;
  int steps;
  double deltaeta,maxdelta;
  double logeta,neweta,sum,norm,largest;
  double d[N];
  double logpre[N];
  double logqun[N];
  double ***logetaun;

  // Make space for the new etas

  logetaun = malloc(G.nvertices*sizeof(double**));
  for (u=0; u<G.nvertices; u++) {
    logetaun[u] = malloc(G.vertex[u].degree*sizeof(double*));
    for (i=0; i<G.vertex[u].degree; i++) {
      logetaun[u][i] = malloc(N*sizeof(double));
    }
  }

  steps = 0;
  do {

    /* Calculate the expected group degrees */

    for (r=0; r<N; r++) {
      d[r] = 0.0;
      for (u=0; u<G.nvertices; u++) d[r] += q[u][r]*G.vertex[u].degree;
    }

    /* Calculate the log-prefactors (without the leading factor of d_i or
     * the prior) */

    for (r=0; r<N; r++) {
      logpre[r] = 0.0;
      for (s=0; s<N; s++) logpre[r] -= omega[r][s]*d[s];
    }

    /* Calculate new values for the one-vertex marginals */

    // fprintf(stderr,"Calculating one-vertex marginals...    \n");
    for (u=0; u<G.nvertices; u++) {
      for (r=0; r<N; r++) {
	logqun[r] = log(prior(r,meta[u])) + G.vertex[u].degree*logpre[r];
	for (i=0; i<G.vertex[u].degree; i++) {
	  sum = 0.0;
	  for (s=0; s<N; s++) sum += eta[u][i][s]*omega[r][s];
	  logqun[r] += log(sum);
	}
	if (r==0) largest = logqun[r];
	else if (logqun[r]>largest) largest = logqun[r];
      }

      /* Normalize */

      norm = 0.0;
      for (r=0; r<N; r++) {
	logqun[r] -= largest;
	norm += exp(logqun[r]);
      }
      for (r=0; r<N; r++) q[u][r] = exp(logqun[r])/norm;
    }

    /* Calculate (unnormalized) new values for the (log) messages */

    for (u=0; u<G.nvertices; u++) {
      if (u%CHECKINT==0) {
	// fprintf(stderr,"Calculating messages... node %i    \n",u);
      }
      for (i=0; i<G.vertex[u].degree; i++) {
	v = G.vertex[u].edge[i].target;
	for (r=0; r<N; r++) {
	  logeta = log(prior(r,meta[v])) + G.vertex[v].degree*logpre[r];
	  for (j=0; j<G.vertex[v].degree; j++) {
	    if (G.vertex[v].edge[j].target!=u) {
#if N==2
	      logeta += log(eta[v][j][0]*omega[r][0]
                            + eta[v][j][1]*omega[r][1]);
#else
	      sum = 0.0;
	      for (s=0; s<N; s++) sum += eta[v][j][s]*omega[r][s];
	      logeta += log(sum);
#endif
	    }
	  }
	  logetaun[u][i][r] = logeta;
	}
      }
    }

    /* Normalize and calculate largest change */

    // fprintf(stderr,"Normalizing messages...\n");
    maxdelta = 0.0;
    for (u=0; u<G.nvertices; u++) {
      for (i=0; i<G.vertex[u].degree; i++) {
	norm = 0.0;
	largest = logetaun[u][i][0];
	for (r=1; r<N; r++) {
	  if (logetaun[u][i][r]>largest) largest = logetaun[u][i][r];
	}
	for (r=0; r<N; r++) {
	  logetaun[u][i][r] -= largest;
	  norm += exp(logetaun[u][i][r]);
	}
	for (r=0; r<N; r++) {
	  neweta = exp(logetaun[u][i][r])/norm;
	  deltaeta = fabs(neweta-eta[u][i][r]);
	  if (deltaeta>maxdelta) maxdelta = deltaeta;
	  eta[u][i][r] = neweta;
	}
      }
    }

    fprintf(stderr,"BP step %i, max change = %g                   \n",
	    steps,maxdelta);
    // if (maxdelta>0.999) break;

  } while ((maxdelta>BP_ACC)&&(++steps<=BP_MAXSTEPS));

  // Free space

  for (u=0; u<G.nvertices; u++) {
    for (i=0; i<G.vertex[u].degree; i++) free(logetaun[u][i]);
    free(logetaun[u]);
  }
  free(logetaun);

  return steps;
}


/* Function to calculate Q's for group r, metadata x */

void getQs(int r, double x, double Q[])
{
  int k;
  double norm;
  double factor,ratio;

  // Avoid divide by zero

  if (x<EPSILON) {
    for (k=0; k<M; k++) Q[k] = 0.0;
    Q[M] = 1.0;
    return;
  }

  // Otherwise...

  factor = pow(x,M);
  ratio = (1-x)/x;
  norm = 0.0;
  for (k=0; k<=M; k++) {
    Q[k] = gmma[r][k]*binom[k]*factor;
    norm += Q[k];
    factor *= ratio;
  }
  for (k=0; k<=M; k++) Q[k] /= norm;
}


/* Function to calculate new parameter values */


double params()
{
  int u,v;
  int i,j,k;
  int r,s;
  double norm,esum,quvrs;
  double deltag,maxdeltag;
  double Q[M+1];
  double gunnorm[N][M+1];
  double oldg[N][M+1];
  double d[N];
  double term[N][N];
  double sum[N][N];
  double L;

  // Loop to calculate the new values of the gammas

  fprintf(stderr,"Calculating parameters...\n");
  do {

    // Zero out the unnormalized estimates

    for (r=0; r<N; r++) {
      for (k=0; k<=M; k++) gunnorm[r][k] = 0.0;
    }

    // Calculate the unnormalized estimates

    for (u=0; u<G.nvertices; u++) {
      for (r=0; r<N; r++) {
	getQs(r,meta[u],Q);
	for (k=0; k<=M; k++) gunnorm[r][k] += q[u][r]*Q[k];
      }
    }

    // Save the old gammas

    for (r=0; r<N; r++) {
      for (k=0; k<=M; k++) oldg[r][k] = gmma[r][k];
    }

    // Calculate the new gammas

    for (k=0; k<=M; k++) {
      norm = 0.0;
      for (r=0; r<N; r++) norm += gunnorm[r][k];
      for (r=0; r<N; r++) gmma[r][k] = gunnorm[r][k]/norm;
    }

    // Find the biggest change

    maxdeltag = 0.0;
    for (r=0; r<N; r++) {
      for (k=0; k<=M; k++) {
	deltag = fabs(oldg[r][k]-gmma[r][k]);
	if (deltag>maxdeltag) maxdeltag = deltag;
      }
    }

    fprintf(stderr,"delta = %g    \n",maxdeltag);

  } while (maxdeltag>GMMA_ACC);

  fprintf(stderr,"\n");

  // Calculate the new values of the omegas

  // Calculate group degrees

  for (r=0; r<N; r++) {
    d[r] = 0.0;
    for (u=0; u<G.nvertices; u++) d[r] += q[u][r]*G.vertex[u].degree;
  }

  // Zero out the sum variables

  for (r=0; r<N; r++) {
    for (s=0; s<N; s++) sum[r][s] = 0.0;
  }

  // Perform the sums

  for (u=0; u<G.nvertices; u++) {
    for (i=0; i<G.vertex[u].degree; i++) {
      v = G.vertex[u].edge[i].target;

      // Find which edge leads back from v to u

      for (j=0; j<G.vertex[v].degree; j++) {
	if (G.vertex[v].edge[j].target==u) break;
      }
      if (j==G.vertex[v].degree) {
	fprintf(stderr,"Error!\n");
	exit(23);
      }

      // Calculate the terms and the normalization factor

      norm = 0.0;
      for (r=0; r<N; r++) {
	for (s=0; s<N; s++) {
	  term[r][s] = omega[r][s]*eta[u][i][r]*eta[v][j][s];
	  norm += term[r][s];
	}
      }

      // Add to the running sums

      for (r=0; r<N; r++) {
    	for (s=0; s<N; s++){
            quvrs = term[r][s]/norm;
            sum[r][s] += quvrs;

            // avoid 0*log(0)
            if(quvrs == 0){
                esum += EPSILON*log(EPSILON);
            } else{
                esum += quvrs*log(quvrs);
            }

          }
      }
    }
    }

  // Calculate the new values of the omega variables

  for (r=0; r<N; r++) {
    for (s=0; s<N; s++) omega[r][s] = sum[r][s]/(d[r]*d[s]);
  }

  // Calculate the expected log-likelihood

  // Internal energy first

  L = 0.0;
  for (r=0; r<N; r++) {
    for (s=0; s<N; s++){
        L += 0.5*sum[r][s]*log(omega[r][s]);
    }
  }

  for (u = 0; u<G.nvertices; u++){
      for (r = 0; r<N; r++){
          L += prior(r, meta[u])*log(prior(r, meta[u]));
      }
  }

  // Now the entropy

  L -= 0.5*esum;
  for (u=0; u<G.nvertices; u++) {
    for (r=0; r<N; r++) {
      if ((q[u][r]>0.0)&&(G.vertex[u].degree>0)) {
	L += (G.vertex[u].degree-1)*q[u][r]*log(q[u][r]);
      }
    }
  }
  return L;
}




int main(int argc, char *argv[])
{
  int u,i,k,r,s;
  int step;
  int bpsteps;
  double c[N][N];
  double oldc[N][N];
  double deltac,maxdelta;
  double L;
  FILE *netfile;
  FILE *metafile;

  // Initialize random number generator

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,time(NULL));

  // Calculate the binomial coefficients

  binoms();

  // Read the network and the metadata from stdin

  netfile = fopen(argv[1],"r");
  metafile = fopen(argv[2],"r");
  fprintf(stderr,"Reading network...\n");
  read_network(&G,netfile);
  for (u=twom=0; u<G.nvertices; u++) twom += G.vertex[u].degree;
  fprintf(stderr,"Read %i nodes and %i edges\n",G.nvertices,twom/2);
  fprintf(stderr,"Reading metadata...\n");
  get_metadata(metafile);

  // Make space for the messages and marginals, and initialize both to
  // random initial values

  eta = malloc(G.nvertices*sizeof(double**));
  q = malloc(G.nvertices*sizeof(double*));

  for (u=0; u<G.nvertices; u++) {

    eta[u] = malloc(G.vertex[u].degree*sizeof(double*));
    for (i=0; i<G.vertex[u].degree; i++) {
      eta[u][i] = malloc(N*sizeof(double));
      for (r=0; r<N; r++) eta[u][i][r] = gsl_rng_uniform(rng);
    }

    q[u] = malloc(N*sizeof(double));
    for (r=0; r<N; r++) q[u][r] = gsl_rng_uniform(rng);

  }

  // Choose random initial values for the parameters, being careful to
  // make sure c[r][s] is symmetric

#define GMIN 0.2
#define GMAX 0.8

  for (r=0; r<N; r++) {
    for (k=0; k<=M; k++) {
      gmma[r][k] = GMIN + (GMAX-GMIN)*gsl_rng_uniform(rng);
    }
  }

  // Choose random values for the omegas, but with a bias toward
  // assortative choices (change if necessary for other networks)

  for (r=0; r<N; r++) {
    for (s=0; s<N; s++) {
        if (r==s) c[r][s] = 1 + gsl_rng_uniform(rng);
        else if (r<s) c[r][s] = gsl_rng_uniform(rng);
        else c[r][s] = c[s][r];
        omega[r][s] = c[r][s]/twom;
    }
  }

  // EM loop

  fprintf(stderr,"Starting EM algorithm...\n");
  step = 1;
  do {

    // Run BP to calculate the messages and one-vertex marginals

    bpsteps = bp();

    // Calculate the new values of the parameters

    L = params();
    if(isnan(L)) return 1;
    for (r=0; r<N; r++) {
      for (s=0; s<N; s++) {
	oldc[r][s] = c[r][s];
	c[r][s] = omega[r][s]*twom;
      }
    }

    // Print out new values of the parameters

    fprintf(stderr,"gamma =\n");
    for (r=0; r<N; r++) {
      for (k=0; k<=M; k++) fprintf(stderr," %g",gmma[r][k]);
      fprintf(stderr,"\n");
    }

    fprintf(stderr,"c =\n");
    for (r=0; r<N; r++) {
      for (s=0; s<N; s++) fprintf(stderr," %g",c[r][s]);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"Log-likelihood = %g\n",L);
    fprintf(stderr,"\n");

    // Find the largest change in any of the c's

    maxdelta = 0.0;
    for (r=0; r<N; r++) {
      for (s=0; s<N; s++) {
	deltac = fabs(c[r][s]-oldc[r][s]);
        if (deltac>maxdelta) maxdelta = deltac;
      }
    }
    if (++step>EM_MAXSTEP) {
          fprintf(stderr,"Solution failed to converge in %i EM steps\n",
              EM_MAXSTEP);
        fail = 1;
        break;
    }

  } while (maxdelta>EM_ACC);

  if (bpsteps>BP_MAXSTEPS) {
    fprintf(stderr,"BP failed converge on final EM step\n");
    fail = 1;
  }

  // Print out results

  for (u=0; u<G.nvertices; u++) {
    printf("%i %g",u,meta[u]);
    for (r=0; r<N; r++) printf(" %g",q[u][r]);
    printf("\n");
  }
  return fail;
}
