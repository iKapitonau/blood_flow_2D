#ifndef GRID_H_
#define GRID_H_

#include <stddef.h>

#define EPS	1e-6

typedef enum print_option_t {
	GRID_PRINT_AS_TABLE =	(1 << 0),
	GRID_PRINT_FOR_CONFIG =	(1 << 1),
	GRID_PRINT_BINARY =		(1 << 2),
} print_option;

typedef struct grid_node_t {
	double u;
	double w;
	double p;
	double mu;
	double dx;
	double dz;
} grid_node;

typedef enum grid_offset_t {
	GRID_OFFSET_U =		offsetof(grid_node, u),
	GRID_OFFSET_W =		offsetof(grid_node, w),
	GRID_OFFSET_P =		offsetof(grid_node, p),
	GRID_OFFSET_MU =	offsetof(grid_node, mu),
	GRID_OFFSET_DX =	offsetof(grid_node, dx),
	GRID_OFFSET_DZ =	offsetof(grid_node, dz),
	GRID_OFFSET_NUM =	6
} grid_offset;

grid_node *grid_generate(double sigma, size_t *n, size_t *m);
void grid_copy(grid_node *src, grid_node **dst);
void grid_destroy(grid_node *grid);

void grid_fill_from_config(grid_node *grid);
void grid_fill(grid_node *grid, double *a, grid_offset offset);

void grid_print_all(grid_node *grid, size_t n, size_t m, const char *filename, print_option mode);
void grid_print_elem(grid_node *grid, size_t n, size_t m, print_option mode, grid_offset offset, const char *filename);

void grid_clear(void);

#endif	// GRID_H_
