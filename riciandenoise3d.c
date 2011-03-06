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
 * tolerance, the method stops when ||u^i - u^i-1||_inf < Tol.
 *
 *======================================================================*/
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "papi.h"

/* Method Parameters */
#define DT 5.0
#define EPSILON 1.0E-20
#define MAX_ITERATIONS 500

/* Macro functions */
#define SQR(x) ((x)*(x))

/* Macros for referring to pixel neighbors */
#define CENTER   (m+M*(n+N*p))
#define RIGHT    (m+M*(n+N*p)+M)
#define LEFT     (m+M*(n+N*p)-M)
#define DOWN     (m+M*(n+N*p)+1)
#define UP       (m+M*(n+N*p)-1)
#define ZOUT     (m+M*(n+N*p+N))
#define ZIN      (m+M*(n+N*p-N))

/* Total variation minimization for rician denoising */
static void riciandenoise3(double *u, const double *f,
			   int M, int N, int P,
			   double sigma, double lambda, double Tol)
{
	double *g;		/* Array storing 1/|grad u| approximation */
	double sigma2, gamma, r, ulast, numer, denom;
	int converged;
	int i, m, n, p;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;
	converged = 0;
	memcpy(u, f, sizeof(double) * M * N * P);	/* Initialize u = f */
	g = calloc(M * N * P, sizeof(double));	/* Allocate temporary work array */

    /*** Main gradient descent loop ***/
	for (i = 1; i <= MAX_ITERATIONS; i++) {
		/* Approximate g = 1/|grad u| */
		for (p = 1; p < P - 1; p++) {
			for (n = 1; n < N - 1; n++) {
				for (m = 1; m < M - 1; m++) {
					g[CENTER] = 1.0 / sqrt(EPSILON
							       + SQR(u[CENTER] -
								     u[RIGHT])
							       + SQR(u[CENTER] -
								     u[LEFT])
							       + SQR(u[CENTER] -
								     u[DOWN])
							       + SQR(u[CENTER] -
								     u[UP])
							       + SQR(u[CENTER] -
								     u[ZOUT])
							       + SQR(u[CENTER] -
								     u[ZIN]));
				}
			}
		}

		/* Update u by a sem-implict step */
		converged = 1;

		for (p = 1; p < P - 1; p++) {
			for (n = 1; n < N - 1; n++) {
				for (m = 1; m < M - 1; m++) {
					/* Evaluate r = I1(u*f/sigma^2) / I0(u*f/sigma^2) with
					 * a cubic rational approximation. */
					r = u[CENTER] * f[CENTER] / sigma2;
					numer =
					    r * 2.38944 + r * (0.950037 + r);
					denom =
					    4.65314 + r * (2.57541 +
							   r * (1.48937 + r));
					r = numer / denom;

					/* Update u */
					ulast = u[CENTER];
					numer =
					    u[CENTER] +
					    DT * (u[RIGHT] * g[RIGHT] +
						  u[LEFT] * g[LEFT] +
						  u[DOWN] * g[DOWN] +
						  u[UP] * g[UP] +
						  u[ZOUT] * g[ZOUT] +
						  u[ZIN] * g[ZIN] +
						  gamma * f[CENTER] * r);
					denom =
					    1.0 + DT * (g[RIGHT] + g[LEFT] +
							g[DOWN] + g[UP] +
							g[ZOUT] + g[ZIN] +
							gamma);
					u[CENTER] = numer / denom;

					/* Test for convergence */
					if (fabs(ulast - u[CENTER]) > Tol) {
						converged = 0;
					}
				}
			}
		}

		if (converged) {
			break;
		}
	}

	if (converged) {
		printf("Converged in %d iterations with tolerance %g.\n",
		       i, Tol);
	} else {
		printf("Maximum iteraations exceeded (Max=%d).\n",
		       MAX_ITERATIONS);
	}

	free(g);
}

void usage()
{
	printf
	    ("riciandenoise3d [-m|n|p dimensions|-i input|-o output|-b batch id]\n");
}

int main(int argc, char *argv[])
{
	int M, N, P, m, n, p, c;
	FILE *inputfile, *outputfile;
	unsigned batch_id = 0;
	double *f, *u;
	unsigned short *sf;
	double sigma = 0.05;
	double lamda = 0.065;
	double Tol = 2e-3;

	if (argc < 11) {
		usage();
		exit(0);
	}

	while ((c = getopt(argc, argv, "vhm:n:p:i:o:b:")) != -1) {
		switch (c) {
		case 'v':
		case 'h':
		case '?':
		default:
			usage();
			exit(0);
		case 'm':
			M = atoi(optarg);
			break;
		case 'n':
			N = atoi(optarg);
			break;
		case 'p':
			P = atoi(optarg);
			break;
		case 'i':
			inputfile = fopen(optarg, "r");
			break;
		case 'o':
			outputfile = fopen(optarg, "w");
			break;
		case 'b':
			sscanf(optarg, "%u", &batch_id);
			break;
		}
	}

	if (M < 1 || N < 1 || P < 1 || !inputfile || !outputfile) {
		usage();
		exit(0);
	}

	f = calloc(M * N * P, sizeof(double));
	sf = calloc(M * N * P, sizeof(unsigned short));
	u = calloc(M * N * P, sizeof(double));
	fread(sf, sizeof(unsigned short), M * N * P, inputfile);
	for (p = 0; p < P; p++) {
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				f[CENTER] = (double)sf[CENTER];
			}
		}
	}

	int Events[5];
	u_long_long papi_values[5];
	util_start_papi(batch_id, Events);
	riciandenoise3(u, f, M, N, P, sigma, lamda, Tol);
	util_stop_papi(batch_id, papi_values);
	util_print_papi(batch_id, papi_values, (batch_id == 0));
	fwrite(u, sizeof(double), M * N * P, outputfile);
	printf("Finished\r\n");
}
