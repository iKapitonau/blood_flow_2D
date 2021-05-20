#ifndef CALC_H_
#define CALC_H_

#include "grid.h"

#include <stddef.h>

typedef struct tma_info_t {
	double *a;
	double *b;
	double *c;
	double *d;
	double *y;
	double xi_1;
	double xi_2;
	double mu_1;
	double mu_2;
	size_t n;
} tma_info;

typedef struct prs_info_t {
	size_t n;
	size_t m;
	size_t start_half_1;
	size_t start_half_2;
	tma_info *ti_half_1;
	tma_info *ti_half_2;
} prs_info;

void tma(tma_info *ti);

void prs_half1(prs_info *pi);

void prs_half2(prs_info *pi);

void calculate(double sigma, double t_beg, double t_end);

#endif	// CALC_H_
