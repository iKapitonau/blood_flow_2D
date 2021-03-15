#include "config.h"
#include "logger.h"
#include "grid.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static grid_node *grid;

static double L;
static double H;
static size_t N;
static size_t M;

static inline double f(double x, double sigma)
{
	return exp(-(x * x) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
}

static void normalize(double *a, size_t n)
{
	log_write(LOG_DEBUG, "Entering normalize function");
	double sum = 0;
	for (size_t i = 0; i < n; ++i)
		sum += a[i];
	
	sum = 1 / sum;
	for (size_t i = 0; i < n; ++i)
		a[i] *= sum;
	log_write(LOG_DEBUG, "Exiting normalize function");
}

static void grid_generate_dimension(double *a, size_t n, size_t size, double sigma)
{
	log_write(LOG_DEBUG, "Entering grid_generate_dimension function");
	double h = 10.0 / n;
	size_t idx = 0;

	for (double i = -5; fabs(5 - i) > EPS; i += h, idx++)
		a[idx] = f(i + h / 2, sigma) * h;

	normalize(a, n);

	for (size_t i = 0; i < n; ++i)
		a[i] = a[i] * size;
	log_write(LOG_DEBUG, "Exiting grid_generate_dimension function");
}

static void init_grid(double *x, double *z)
{
	log_write(LOG_DEBUG, "Entering init_grid function");
	size_t nodes_number;
	double *u;
	double *w;
	double *p;

	get_velocity_x(&nodes_number, &u);	
	get_velocity_z(&nodes_number, &w);	
	get_pressure(&nodes_number, &p);

	grid = malloc(nodes_number * sizeof(grid_node));

	for (size_t i = 0; i < nodes_number; ++i) {
		grid[i].u = u[i];
		grid[i].w = w[i];
		grid[i].p = p[i];
		grid[i].mu = 0;
	}
	
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < M; ++j) {
			grid[i * M + j].dx = x[j];
			grid[i * M + j].dz = z[i];
		}
	}
	log_write(LOG_DEBUG, "Exiting init_grid function");
}

static void print_gnuplot(double *a, size_t n, double *b, size_t m)
{
	log_write(LOG_DEBUG, "Entering print_gnuplot function");
	double sum = 0;
	for (size_t i = 0; i < n; ++i) {
		printf("set arrow from 0,%.18lf to %lf,%.18lf nohead\n", sum, L, sum);
		sum += a[i];
	}
	printf("set arrow from 0,%.18lf to %lf,%.18lf nohead\n", sum, L, sum);

	sum = 0;
	for (size_t i = 0; i < m; ++i) {
		printf("set arrow from %.18lf,0 to %.18lf,%lf nohead\n", sum, sum, H);
		sum += b[i];
	}
	printf("set arrow from %.18lf,0 to %.18lf,%lf nohead\n", sum, sum, H);
	printf("plot 0\npause -1");
	log_write(LOG_DEBUG, "Exiting print_gnuplot function");
}

void get_grid(grid_node **grid_, size_t *n, size_t *m)
{
	log_write(LOG_DEBUG, "Entering get_grid function");
	*grid_ = grid;
	*n = N + 1;
	*m = M + 1;
	log_write(LOG_DEBUG, "Exiting get_grid function");
}

void grid_generate(double sigma)
{
	log_write(LOG_DEBUG, "Entering grid_generate function");

	L = get_vessel_size_x();
	H = get_vessel_size_z();
	N = get_rows_number();
	M = get_columns_number();

	double *z = malloc(N * sizeof(double));
	double *x = malloc(M * sizeof(double));

	log_write(LOG_INFO, "Generating x dimension...");
	grid_generate_dimension(x, M, L, sigma);
	log_write(LOG_INFO, "Generating z dimension...");
	grid_generate_dimension(z, N, H, sigma);

	log_write(LOG_INFO, "Initializing grid...");
	init_grid(x, z);

	//print_gnuplot(z, N, x, M);

	free(z);
	free(x);
	log_write(LOG_DEBUG, "Exiting grid_generate function");
}

void grid_destroy(void)
{
	log_write(LOG_DEBUG, "Entering grid_destroy function");
	free(grid);
	log_write(LOG_DEBUG, "Exiting grid_destroy function");
}
