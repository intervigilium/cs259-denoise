/*========================================================================
 *
 * RICIANDENOIES3MX.C  3D TV minimization for Rician denoising
 *
 * u = riciandenoise3mx(f,sigma,lambda,Tol) performs denoising on a 3D 
 * volume f with Rician noise with parameter sigma.  The denoised image u
 * is found as the minimizer of 
 *
 *         /                      / [ u^2 + f^2            u f    ]
 *    min  | |grad u| dx + lambda | [ --------- - log I0( ----- ) ] dx.
 *     u   /                      / [ 2 sigma^2          sigma^2  ]
 *
 * Parameter lambda >= 0 determines the strength of the denoising: smaller
 * lambda implies stronger denoising.  Tol specifies the stopping 
 * tolerance, the method stops when ||u^Iter - u^Iter-1||_inf < Tol.
 *
 *======================================================================*/
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "papi.h"

/* Method Parameters */
#define DT          5.0
#define EPSILON     1.0E-20
#define MAXITER     500

/* Macro functions */
#define SQR(x) ((x)*(x))

/* Total variation minimization for rician denoising */
static void riciandenoise3(double *u, const double *f,
			   int M, int N, int P,
			   double sigma, double lambda, double Tol)
{
    double *g;			/* Array storing 1/|grad u| approximation */
    double sigma2, gamma, r, ulast;
    bool Converged;
    int m, n, p;
    int Iter;

    /* Initializations */
    sigma2 = SQR(sigma);
    gamma = lambda / sigma2;
    Converged = false;
    memcpy(u, f, sizeof(double) * M * N * P);	/* Initialize u = f */
    g = calloc(M * N * P, sizeof(double));	/* Allocate temporary work array */

    /*** Main gradient descent loop ***/
    for (Iter = 1; Iter <= MAXITER; Iter++) {
	/* Macros for referring to pixel neighbors */
#define CENTER   (m+M*(n+N*p))
#define RIGHT    (m+M*(n+N*p)+M)
#define LEFT     (m+M*(n+N*p)-M)
#define DOWN     (m+M*(n+N*p)+1)
#define UP       (m+M*(n+N*p)-1)
#define ZOUT     (m+M*(n+N*p+N))
#define ZIN      (m+M*(n+N*p-N))

	/* Approximate g = 1/|grad u| */
	for (p = 1; p < P - 1; p++)
	    for (n = 1; n < N - 1; n++)
		for (m = 1; m < M - 1; m++)
		    g[CENTER] = 1.0 / sqrt(EPSILON
					   + SQR(u[CENTER] - u[RIGHT])
					   + SQR(u[CENTER] - u[LEFT])
					   + SQR(u[CENTER] - u[DOWN])
					   + SQR(u[CENTER] - u[UP])
					   + SQR(u[CENTER] - u[ZOUT])
					   + SQR(u[CENTER] - u[ZIN]));

	/* Update u by a sem-implict step */
	Converged = true;

	for (p = 1; p < P - 1; p++)
	    for (n = 1; n < N - 1; n++)
		for (m = 1; m < M - 1; m++) {
		    /* Evaluate r = I1(u*f/sigma^2) / I0(u*f/sigma^2) with
		       a cubic rational approximation. */
		    r = u[CENTER] * f[CENTER] / sigma2;
		    r = (r * (2.38944 + r * (0.950037 + r)))
			/ (4.65314 + r * (2.57541 + r * (1.48937 + r)));
		    /* Update u */
		    ulast = u[CENTER];
		    u[CENTER] = (u[CENTER] + DT * (u[RIGHT] * g[RIGHT]
						   + u[LEFT] * g[LEFT] +
						   u[DOWN] * g[DOWN] +
						   u[UP] * g[UP]
						   + u[ZOUT] * g[ZOUT] +
						   u[ZIN] * g[ZIN]
						   +
						   gamma * f[CENTER] *
						   r)) / (1.0 +
							  DT * (g[RIGHT] +
								g[LEFT]
								+ g[DOWN] +
								g[UP]
								+ g[ZOUT] +
								g[ZIN] +
								gamma));

		    /* Test for convergence */
		    if (fabs(ulast - u[CENTER]) > Tol)
			Converged = false;
		}

	if (Converged)
	    break;
    }

    /* Done, show exiting message */
    if (Converged)
	printf("Converged in %d iterations with tolerance %g.\n",
	       Iter, Tol);
    else
	printf("Maximum iterations exceeded (MaxIter=%d).\n", MAXITER);

    free(g);			/* Free temporary array */
    return;
}

int main(int argc, char *argv[])
{
    if (argc < 6) {
	printf("riciandnoise3 M N P inputfile outputfile <batch id>\r\n");
	exit(0);
    }
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    int P = atoi(argv[3]);
    int p, n, m;
    FILE *inputfile = fopen(argv[4], "r");
    FILE *outputfile = fopen(argv[5], "w");
    unsigned batch_id = 0;
    if (argc == 7) {
	batch_id = atoi(argv[6]);
    }
    double *f, *u;
    unsigned short *sf;
    f = calloc(M * N * P, sizeof(double));
    sf = calloc(M * N * P, sizeof(unsigned short));
    u = calloc(M * N * P, sizeof(double));
    fread(sf, sizeof(unsigned short), M * N * P, inputfile);
    for (p = 0; p < P; p++) {
	for (n = 0; n < N; n++) {
	    for (m = 0; m < M; m++) {
		f[CENTER] = (double) sf[CENTER];
	    }
	}
    }
    double sigma = 0.05;
    double lamda = 0.065;
    double Tol = 2e-3;
    int Events[5];
    u_long_long papi_values[5];
    util_start_papi(batch_id, Events);
    riciandenoise3(u, f, M, N, P, sigma, lamda, Tol);
    util_stop_papi(batch_id, papi_values);
    util_print_papi(batch_id, papi_values, (batch_id == 0));
    fwrite(u, sizeof(double), M * N * P, outputfile);
    printf("Finished\r\n");
}
