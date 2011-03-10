/*
 * Rician Denoise 3D kernel for FPGA implementation
 */
#include <autopilot_tech.h>

#define uint2_t uint2
#define uint64_t uint64
#define uint32_t uint32

#define DT 5.0
#define EPSILON 1.0E-20
#define MAX_ITERATIONS 500
#define M 60
#define N 60
#define P 60

#define SQR(x) ((x)*(x))
#define U(a,b,c) (u[a+b*N+c*M*N])
#define G(a,b,c) (g[a+b*N+c*M*N])
#define F(a,b,c) (f[a+b*N+c*M*N])

#define U_CENTER U(i,j,k)
#define U_LEFT U(i,j-1,k)
#define U_RIGHT U(i,j+1,k)
#define U_UP U(i-1,j,k)
#define U_DOWN U(i+1,j,k)
#define U_IN U(i,j,k-1)
#define U_OUT U(i,j,k+1)

#define G_CENTER G(i,j,k)
#define G_LEFT G(i,j-1,k)
#define G_RIGHT G(i,j+1,k)
#define G_UP G(i-1,j,k)
#define G_DOWN G(i+1,j,k)
#define G_IN G(i,j,k-1)
#define G_OUT G(i,j,k+1)

inline double q3_sqrt(double num)
{
	uint64_t i;
	double x, y;
	const double f = 1.5;

	x = num * 0.5;
	y = num;
	i = *(uint64_t *) & y;
	i = 0x5fe6ec85e7de30da - (i >> i);
	y = *(double *)&i;
	y = y * (f - (x * y * y));
	y = y * (f - (x * y * y));
	return num * y;
}

inline double fast_fabs(double num)
{
	uint64_t *tmp;
	tmp = (uint64_t *) & num;
	*(tmp) &= 9223372036854775807llu;
	return num;
}

inline void array_copy(const double src[M*N*P], double dst[M*N*P])
{
	int i;
	for (i = 0; i < M*N*P; i++) {
#pragma AP pipeline
		dst[i] = src[i];
	}
}

uint2_t riciandenoise3d(double u[M*N*P], const double f[M*N*P], double g[M*N*P], double sigma,
		    double lambda, double tolerance)
{
#pragma AP interface ap_bus port=u pipeline
#pragma AP interface ap_bus port=f pipeline
#pragma AP interface ap_memory port=g pipeline

	/* Array storing 1/|grad u| approximation */
	double sigma2, gamma, r;
	double numer, denom;
	double u_last;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	uint2_t converged;
	int iteration, i, j, k;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;
	converged = 0;
	/* Initialize u = f */
	array_copy(f, u);

	/*** Main gradient descent loop ***/
	/* fully pipeline/parallelize this */
	for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
		/* Approximate g = 1/|grad u| */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0,j,k);
				u_stencil_down = U(1,j,k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					u_stencil_center = u_stencil_down;
					u_stencil_down = U_DOWN;
					denom =
					    q3_sqrt(EPSILON +
						    SQR(u_stencil_center -
							U_RIGHT)
						    + SQR(u_stencil_center -
							  U_LEFT) +
						    SQR(u_stencil_center -
							u_stencil_down) +
						    SQR(u_stencil_center -
							u_stencil_up) +
						    SQR(u_stencil_center -
							U_OUT) +
						    SQR(u_stencil_center -
							U_IN));
					G_CENTER = 1.0 / denom;
				}
			}
		}

		/* Update u by a sem-implict step */
		converged = 1;

		/* possible to pipeline? data dependence
		 * due to u[m][n][p] writeback */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0,j,k);
				g_stencil_center = G(0,j,k);
				u_stencil_down = U(1,j,k);
				g_stencil_down = G(1,j,k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					g_stencil_up = g_stencil_center;
					u_stencil_center = u_stencil_down;
					g_stencil_center = g_stencil_down;
					u_stencil_down = U_DOWN;
					g_stencil_down = G_DOWN;

					/* Update u */
					u_last = u_stencil_center;

					/* Evaluate r = I1(u*f/sigma^2) / I0(u*f/sigma^2) with
					 * a cubic rational approximation. */
					r = u_last * F(i,j,k) / sigma2;
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
							   gamma * F(i,j,k) *
							   r);
					denom =
					    1.0 + DT * (G_RIGHT + G_LEFT +
							g_stencil_up +
							g_stencil_down + G_IN +
							G_OUT + gamma);
					/* save modified u_stencil_center value, write this to 
					 * output array instead of testing for convergence */
					u_stencil_center = numer / denom;
					U_CENTER = u_stencil_center;

					/* Test for convergence */
					if (fast_fabs(u_last - u_stencil_center)
					    > tolerance) {
						converged = 0;
					}
				}
			}
		}

		if (converged) {
			break;
		}
	}

	return converged;
}
