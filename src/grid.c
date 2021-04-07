#include "config.h"
#include "logger.h"
#include "grid.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

static const grid_offset ALL_OFFSETS[GRID_OFFSET_NUM] = {
	GRID_OFFSET_U,
	GRID_OFFSET_W,
	GRID_OFFSET_P,
	GRID_OFFSET_MU,
	GRID_OFFSET_DX,
	GRID_OFFSET_DZ
};

static double L;
static double H;
static size_t N;
static size_t M;
static double sigma;
static double *dx;
static double *dz;

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

grid_node *grid_generate(double sigma_, size_t *n, size_t *m)
{
	log_write(LOG_DEBUG, "Entering grid_generate function");

	static bool first_time_in_function = true;

	if (first_time_in_function) {
		log_write(LOG_INFO, "First time in function. Allocating dx, dz...");
		L = get_vessel_size_x();
		H = get_vessel_size_z();
		M = get_rows_number();
		N = get_columns_number();

		dz = malloc(M * sizeof(double));
		dx = malloc(N * sizeof(double));

		first_time_in_function = false;
		log_write(LOG_INFO, "Allocated successfully!");
	}
	
	if (sigma_ != sigma) {
		log_write(LOG_INFO, "New sigma param. Regenerating grid...");
		sigma = sigma_;

		log_write(LOG_INFO, "Generating x dimension...");
		grid_generate_dimension(dx, N, L, sigma);
		log_write(LOG_INFO, "Generating z dimension...");
		grid_generate_dimension(dz, M, H, sigma);
		log_write(LOG_INFO, "Generated successfully!");
	}

	log_write(LOG_INFO, "Allocating grid...");
	grid_node *grid = malloc((N + 1) * (M + 1) * sizeof(grid_node));
	*n = N + 1;
	*m = M + 1;
	log_write(LOG_INFO, "Allocated successfully!");

	//print_gnuplot(z, N, x, M);

	log_write(LOG_DEBUG, "Exiting grid_generate function");

	return grid;
}

void grid_destroy(grid_node *grid)
{
	log_write(LOG_DEBUG, "Entering grid_destroy function");
	free(grid);
	log_write(LOG_DEBUG, "Exiting grid_destroy function");
}

void grid_fill_from_config(grid_node *grid)
{
	log_write(LOG_DEBUG, "Entering grid_fill_from_config function");
	log_write(LOG_INFO, "Initializing grid...");
	size_t nodes_number;
	double *u;
	double *w;
	double *p;

	get_velocity_x(&nodes_number, &u);	
	get_velocity_z(&nodes_number, &w);	
	get_pressure(&nodes_number, &p);

	for (size_t i = 0; i < nodes_number; ++i) {
		grid[i].u = u[i];
		grid[i].w = w[i];
		grid[i].p = p[i];
		grid[i].mu = 0;
	}
	
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < M; ++j) {
			grid[i * M + j].dx = dx[i];
			grid[i * M + j].dz = dz[j];
		}
	}
	log_write(LOG_INFO, "Initialized successfully!");
	log_write(LOG_DEBUG, "Exiting grid_fill_from_config function");
}

void grid_fill(grid_node *grid, double *a, grid_offset offset)
{
	log_write(LOG_DEBUG, "Entering grid_fill function");
	size_t nodes_number = (N + 1) * (M + 1);

	for (size_t i = 0; i < nodes_number; ++i)
		*(double *)((char *)&grid[i] + offset) = a[i];

	log_write(LOG_DEBUG, "Exiting grid_fill function");
}

void grid_print_elem(grid_node *grid, size_t n, size_t m, print_option mode, grid_offset offset)
{
	log_write(LOG_DEBUG, "Entering grid_print_element function");
	if (mode & GRID_PRINT_AS_TABLE) {
		switch (offset) {
		case GRID_OFFSET_U:
			printf("velocity_x:\n");
			break;
		case GRID_OFFSET_W:
			printf("velocity_z:\n");
			break;
		case GRID_OFFSET_P:
			printf("pressure:\n");
			break;
		case GRID_OFFSET_MU:
			printf("mu:\n");
			break;
		case GRID_OFFSET_DX:
			printf("dx:\n");
			break;
		case GRID_OFFSET_DZ:
			printf("dz:\n");
			break;
		default:
			return;
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < m; ++j) {
				printf("%lf\t", *(double *)((char *)&grid[i * m + j] + offset));
			}
			printf("\n");
		}
	} else if (mode & GRID_PRINT_FOR_CONFIG) {
		switch (offset) {
		case GRID_OFFSET_U:
			printf("velocity_x = [ ");
			break;
		case GRID_OFFSET_W:
			printf("velocity_z = [ ");
			break;
		case GRID_OFFSET_P:
			printf("pressure = [ ");
			break;
		case GRID_OFFSET_MU:
			printf("mu = [ ");
			break;
		case GRID_OFFSET_DX:
			printf("dx = [ ");
			break;
		case GRID_OFFSET_DZ:
			printf("dz = [ ");
			break;
		default:
			return;
		}
		for (size_t i = 0; i < n; ++i) 
			for (size_t j = 0; j < m; ++j)
				if (i == n - 1 && j == m - 1)
					printf("%lf ]\n", *(double *)((char *)&grid[i * m + j] + offset));
				else 
					printf("%lf, ", *(double *)((char *)&grid[i * m + j] + offset));
	}

	log_write(LOG_DEBUG, "Exiting grid_print_element function");
}

void grid_print_all(grid_node *grid, size_t n, size_t m, print_option mode)
{
	log_write(LOG_DEBUG, "Entering grid_print_all function");

	if ((mode & GRID_PRINT_AS_TABLE) || (mode & GRID_PRINT_FOR_CONFIG)) {
		for (size_t i = 0; i != GRID_OFFSET_NUM; ++i) {
			grid_print_elem(grid, n, m, mode, ALL_OFFSETS[i]);
		}
	}

	log_write(LOG_DEBUG, "Exiting grid_print_all function");
}
