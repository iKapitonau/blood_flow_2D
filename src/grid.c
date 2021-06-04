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
static FILE *out;

static inline double f(double x, double sigma)
{
	return exp(-(x * x) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
}

static void normalize(double *a, size_t n)
{
	//log_write(LOG_DEBUG, "Entering normalize function");
	double sum = 0;
	for (size_t i = 0; i < n; ++i)
		sum += a[i];
	
	sum = 1 / sum;
	for (size_t i = 0; i < n; ++i)
		a[i] *= sum;
	//log_write(LOG_DEBUG, "Exiting normalize function");
}

static void grid_generate_dimension(double *a, size_t n, double size, double sigma)
{
	//log_write(LOG_DEBUG, "Entering grid_generate_dimension function");
	double h = 10.0 / n;
	size_t idx = 0;

	for (double i = -5; fabs(5 - i) > EPS_GRID; i += h, idx++)
		a[idx] = f(i + h / 2, sigma) * h;

	normalize(a, n);

	for (size_t i = 0; i < n; ++i)
		a[i] = a[i] * size;
	//log_write(LOG_DEBUG, "Exiting grid_generate_dimension function");
}

grid_node *grid_generate(double sigma_, size_t *n, size_t *m)
{
	//log_write(LOG_DEBUG, "Entering grid_generate function");

	static bool first_time_in_function = true;

	if (first_time_in_function) {
		log_write(LOG_INFO, "First time in function. Allocating dx, dz...");
		L = get_vessel_size_x();
		H = get_vessel_size_z();

		M = get_rows_number();
		N = get_columns_number();

		dz = calloc(M + 1, sizeof(double));
		dx = calloc(N + 1, sizeof(double));

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
	grid_node *grid = calloc((N + 1) * (M + 1), sizeof(grid_node));
	*n = N + 1;
	*m = M + 1;
	for (size_t i = 0; i < *n; ++i) {
		for (size_t j = 0; j < *m; ++j) {
			grid[i * *m + j].dx = dx[i];
			grid[i * *m + j].dz = dz[j];
		}
	}
	log_write(LOG_INFO, "Allocated successfully!");

	//log_write(LOG_DEBUG, "Exiting grid_generate function");

	return grid;
}

void grid_copy(grid_node *src, grid_node **dst)
{
	grid_destroy(*dst);

	size_t n, m;
	*dst = grid_generate(sigma, &n, &m);
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			(*dst)[i * m + j] = src[i * m + j];
}

void grid_destroy(grid_node *grid)
{
	//log_write(LOG_DEBUG, "Entering grid_destroy function");
	if (grid)
		free(grid);
	//log_write(LOG_DEBUG, "Exiting grid_destroy function");
}

void grid_fill_from_config(grid_node *grid)
{
	//log_write(LOG_DEBUG, "Entering grid_fill_from_config function");
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

	free(u);
	free(w);
	free(p);
	
	log_write(LOG_INFO, "Initialized successfully!");
	//log_write(LOG_DEBUG, "Exiting grid_fill_from_config function");
}

void grid_fill(grid_node *grid, double *a, grid_offset offset)
{
	//log_write(LOG_DEBUG, "Entering grid_fill function");
	size_t nodes_number = (N + 1) * (M + 1);

	for (size_t i = 0; i < nodes_number; ++i)
		*(double *)((char *)&grid[i] + offset) = a[i];

	//log_write(LOG_DEBUG, "Exiting grid_fill function");
}

void grid_print_elem(grid_node *grid, size_t n, size_t m, print_option mode, grid_offset offset, const char *filename)
{
	//log_write(LOG_DEBUG, "Entering grid_print_element function");
	if (out == NULL) {
		if (filename != NULL) 
			out = fopen(filename, ((mode & GRID_PRINT_BINARY) ? "wb" : "w"));
		else 
			out = stdout;
	}
	if (mode & GRID_PRINT_AS_TABLE) {
		switch (offset) {
		case GRID_OFFSET_U:
			fprintf(out, "velocity_x:\n");
			break;
		case GRID_OFFSET_W:
			fprintf(out, "velocity_z:\n");
			break;
		case GRID_OFFSET_P:
			fprintf(out, "pressure:\n");
			break;
		case GRID_OFFSET_MU:
			fprintf(out, "mu:\n");
			break;
		case GRID_OFFSET_DX:
			fprintf(out, "dx:\n");
			break;
		case GRID_OFFSET_DZ:
			fprintf(out, "dz:\n");
			break;
		default:
			return;
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < m; ++j) {
				fprintf(out, "%lf\t", *(double *)((char *)&grid[i * m + j] + offset));
			}
			fprintf(out, "\n");
		}
	} else if (mode & GRID_PRINT_FOR_CONFIG) {
		switch (offset) {
		case GRID_OFFSET_U:
			fprintf(out, "velocity_x = [ ");
			break;
		case GRID_OFFSET_W:
			fprintf(out, "velocity_z = [ ");
			break;
		case GRID_OFFSET_P:
			fprintf(out, "pressure = [ ");
			break;
		case GRID_OFFSET_MU:
			fprintf(out, "mu = [ ");
			break;
		case GRID_OFFSET_DX:
			fprintf(out, "dx = [ ");
			break;
		case GRID_OFFSET_DZ:
			fprintf(out, "dz = [ ");
			break;
		default:
			return;
		}
		for (size_t i = 0; i < n; ++i) 
			for (size_t j = 0; j < m; ++j)
				if (i == n - 1 && j == m - 1)
					fprintf(out, "%lf ]\n", *(double *)((char *)&grid[i * m + j] + offset));
				else 
					fprintf(out, "%lf, ", *(double *)((char *)&grid[i * m + j] + offset));
	} else if (mode & GRID_PRINT_BINARY) {
		switch (offset) {
		case GRID_OFFSET_U:
			fwrite("velocity_x:\n", sizeof("velocity_x:\n"), 1, out);
			break;
		case GRID_OFFSET_W:
			fwrite("velocity_z:\n", sizeof("velocity_z:\n"), 1, out);
			break;
		case GRID_OFFSET_P:
			fwrite("pressure:\n", sizeof("pressure:\n"), 1, out);
			break;
		case GRID_OFFSET_MU:
			fwrite("mu:\n", sizeof("mu:\n"), 1, out);
			break;
		case GRID_OFFSET_DX:
			fwrite("dx:\n", sizeof("dx:\n"), 1, out);
			break;
		case GRID_OFFSET_DZ:
			fwrite("dz:\n", sizeof("dz:\n"), 1, out);
			break;
		default:
			return;
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < m; ++j) {
				fwrite((double *)((char *)&grid[i * m + j] + offset), sizeof(double), 1, out);
			}
			fwrite("\n", sizeof("\n"), 1, out);
		}
	}

	//log_write(LOG_DEBUG, "Exiting grid_print_element function");
}

void grid_print_all(grid_node *grid, size_t n, size_t m, const char *filename, print_option mode)
{
	//log_write(LOG_DEBUG, "Entering grid_print_all function");

	if ((mode & GRID_PRINT_AS_TABLE) || (mode & GRID_PRINT_FOR_CONFIG) ||
			(mode & GRID_PRINT_BINARY)) {
		for (size_t i = 0; i != GRID_OFFSET_NUM; ++i) {
			grid_print_elem(grid, n, m, mode, ALL_OFFSETS[i], filename);
		}
	}

	//log_write(LOG_DEBUG, "Exiting grid_print_all function");
}

void grid_clear(void)
{
	free(dx);
	free(dz);
	if (out == stdout)
		out = NULL;
	else if (out)
		fclose(out);
}
