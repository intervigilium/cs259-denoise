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

#define SQR(x) ((x)*(x))

/* macros for array access */
#define U_CENTER u[m][n][p]
#define U_LEFT u[m][n-1][p]
#define U_RIGHT u[m][n+1][p]
#define U_UP u[m-1][n][p]
#define U_DOWN u[m+1][n][p]
#define U_IN u[m][n][p-1]
#define U_OUT u[m][n][p+1]

#define G_CENTER g[m][n][p]
#define G_LEFT g[m][n-1][p]
#define G_RIGHT g[m][n+1][p]
#define G_UP g[m-1][n][p]
#define G_DOWN g[m+1][n][p]
#define G_IN g[m][n][p-1]
#define G_OUT g[m][n][p+1]

inline void array_copy(double src[M][N][P], double dst[M][N][P])
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
	double sigma2, gamma, r;
	double numer, denom;
	double u_last;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	int converged;
	int i, m, n, p;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;
	converged = 0;
	/* Initialize u = f */
	array_copy(f, u);

    /*** Main gradient descent loop ***/
	/* fully pipeline/parallelize this */
	for (i = 1; i <= MAX_ITERATIONS; i++) {
		/* Approximate g = 1/|grad u| */
		for (p = 1; p < P - 1; p++) {
			for (n = 1; n < N - 1; n++) {
				for (m = 1; m < M - 1; m++) {
					/* stencil m-1 = m; m = m+1 */
					if (m == 1) {
						u_stencil_up = U_UP;
						u_stencil_center = U_CENTER;
					} else {
						u_stencil_up = u_stencil_center;
						u_stencil_center =
						    u_stencil_down;
					}
					u_stencil_down = U_DOWN;
					denom =
					    sqrt(EPSILON +
						 SQR(u_stencil_center - U_RIGHT)
						 + SQR(u_stencil_center -
						       U_LEFT) +
						 SQR(u_stencil_center -
						     u_stencil_down) +
						 SQR(u_stencil_center -
						     u_stencil_up) +
						 SQR(u_stencil_center - U_OUT) +
						 SQR(u_stencil_center - U_IN));
					g[m][n][p] = 1.0 / denom;
				}
			}
		}

		/* Update u by a sem-implict step */
		converged = 1;

		/* possible to pipeline? data dependence
		 * due to u[m][n][p] writeback */
		for (p = 1; p < P - 1; p++) {
			for (n = 1; n < N - 1; n++) {
				for (m = 1; m < M - 1; m++) {
					/* stencil m-1 = m; m = m+1 */
					if (m == 1) {
						u_stencil_up = U_UP;
						g_stencil_up = G_UP;
						u_stencil_center = U_CENTER;
						g_stencil_center = G_CENTER;
					} else {
						u_stencil_up = u_stencil_center;
						g_stencil_up = g_stencil_center;
						u_stencil_center =
						    u_stencil_down;
						g_stencil_center =
						    g_stencil_down;
					}
					u_stencil_down = U_DOWN;
					g_stencil_down = G_DOWN;

					/* Update u */
					u_last = u_stencil_center;

					/* Evaluate r = I1(u*f/sigma^2) / I0(u*f/sigma^2) with
					 * a cubic rational approximation. */
					r = ulast * f[m][n][p] / sigma2;
					numer =
					    r * (2.38944 + r * (0.950037 + r));
					denom =
					    4.65314 + r * (2.57541 +
							   r * (1.48937 + r));
					r = numer / denom;

					numer =
					    u_last + DT * (U_RIGHT * G_RIGHT +
							   U_LEFT * G_LEFT +
							   u_stencil_up *
							   g_stencil_up +
							   u_stencil_down *
							   g_stencil_down +
							   U_IN * G_IN +
							   U_OUT * G_OUT +
							   gamma * f[m][n][p] *
							   r);
					denom =
					    1.0 + DT * (G_RIGHT + G_LEFT +
							g_stencil_up +
							g_stencil_down + G_IN +
							G_OUT + gamma);
					/* save modified u_stencil_center value */
					u_stencil_center = numer / denom;
					u[m][n][p] = u_stencil_center;

					/* Test for convergence */
					if (fabs(u_last - u_stencil_center) >
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
