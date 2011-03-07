/*
 * Rician Denoise 3D kernel for FPGA implementation
 */

#include <math.h>

#define DT 5.0
#define EPSILON 1.0E-20
#define MAX_ITERATIONS 500
#define M 60
#define N 60
#define P 60

void array_copy(double src[M][N][P], double dst[M][N][P])
{
	int i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < P; k++) {
				dst[i][j][k] = src[i][j][k];
			}
		}
	}
}

void riciandenoise3d(double u[M][N][P], const double f[M][N][P], double sigma,
		     double lambda, double tolerance)
{
	/* Array storing 1/|grad u| approximation */
	double g[M][N][P];

	double sigma2, gamma, r, ulast, numer, denom;
	int converged;
	int i, m, n, p;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;
	converged = 0;
	/* Initialize u = f */
	array_copy(f, u);

    /*** Main gradient descent loop ***/
	for (i = 1; i <= MAX_ITERATIONS; i++) {
		/* Approximate g = 1/|grad u| */
		for (p = 1; p < P - 1; p++) {
			for (n = 1; n < N - 1; n++) {
				for (m = 1; m < M - 1; m++) {
					/* stencil */
					denom = sqrt(EPSILON
						     + pow(u[m][n][p] -
							   u[m + 1][n][p], 2)
						     + pow(u[m][n][p] -
							   u[m - 1][n][p], 2)
						     + pow(u[m][n][p] -
							   u[m][n - 1][p], 2)
						     + pow(u[m][n][p] -
							   u[m][n + 1][p], 2)
						     + pow(u[m][n][p] -
							   u[m][n][p - 1], 2)
						     + pow(u[m][n][p] -
							   u[m][n][p + 1], 2));
					g[m][n][p] = 1.0 / denom;
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
					r = u[m][n][p] * f[m][n][p] / sigma2;
					numer =
					    r * (2.38944 + r * (0.950037 + r));
					denom =
					    4.65314 + r * (2.57541 +
							   r * (1.48937 + r));
					r = numer / denom;

					/* Update u */
					ulast = u[m][n][p];
					/* RIGHT = [m+1][n][p]
					 * LEFT = [m-1][n][p]
					 * UP = [m][n+1][p]
					 * DOWN = [m][n-1][p]
					 * IN = [m][n][p+1]
					 * OUT = [m][n][p-1] */
					numer =
					    u[m][n][p] +
					    DT * (u[m + 1][n][p] *
						  g[m + 1][n][p] + u[m -
								     1][n][p] *
						  g[m - 1][n][p] + u[m][n -
									1][p] *
						  g[m][n - 1][p] + u[m][n +
									1][p] *
						  g[m][n + 1][p] + u[m][n][p -
									   1] *
						  g[m][n][p - 1] + u[m][n][p +
									   1] *
						  g[m][n][p + 1] +
						  gamma * f[m][n][p] * r);
					denom =
					    1.0 + DT * (g[m + 1][n][p] +
							g[m - 1][n][p] +
							g[m][n - 1][p] +
							g[m][n + 1][p] +
							g[m][n][p - 1] +
							g[m][n][p + 1] + gamma);
					u[m][n][p] = numer / denom;

					/* Test for convergence */
					if (fabs(ulast - u[m][n][p]) >
					    tolerance) {
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
		       i, tolerance);
	} else {
		printf("Maximum iterations exceeded (Max: %d)\n",
		       MAX_ITERATIONS);
	}
}
