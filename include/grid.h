#ifndef GRID_H_
#define GRID_H_

#define EPS	1e-6

typedef struct grid_node_t {
	double u;
	double w;
	double p;
	double mu;
	double dx;
	double dz;
} grid_node;

void get_grid(grid_node **grid, size_t *N, size_t *M);

void grid_generate(double sigma);
void grid_destroy(void);

#endif	// GRID_H_
