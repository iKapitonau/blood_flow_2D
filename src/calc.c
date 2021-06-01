#include "blood_consts.h"
#include "calc.h"
#include "config.h"
#include "logger.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS	1e-9

void tma(tma_info *ti)
{
	size_t n = ti->n;
	double *a = ti->a;
	double *b = ti->b;
	double *c = ti->c;
	double *d = ti->d;
	double *y = ti->y;
	double xi_1 = ti->xi_1;
	double xi_2 = ti->xi_2;
	double mu_1 = ti->mu_1;
	double mu_2 = ti->mu_2;
	double alpha[n + 1], beta[n + 1];

	alpha[1] = xi_1;
	beta[1] = mu_1;

	for (size_t i = 1; i < n; ++i) {
		alpha[i + 1] = -c[i] / (b[i] + a[i] * alpha[i]);
		beta[i + 1] = (d[i] - a[i] * beta[i]) / (b[i] + a[i] * alpha[i]);
	}

	y[n] = (xi_2 * beta[n] + mu_2) / (1 - xi_2 * alpha[n]);
	for (ssize_t i = n - 1; i >= 0; --i)
		y[i] = beta[i + 1] + alpha[i + 1] * y[i + 1];
}

void prs_half_1(prs_info *pi)
{
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t j = pi->start_half_1; j < pi->m; ++j)
		tma(&(pi->ti_half_1[j]));
}

void prs_half_2(prs_info *pi)
{
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = pi->start_half_2; i < pi->n; ++i)
		tma(&(pi->ti_half_2[i]));
}

static inline double w_0(double t)
{
	double f = 0;

	if (t < 0.2598 || t >= 0.7523)
		f = 0.189;
	else if (t >= 0.2598 && t < 0.458) 
		f = -59.5 * pow(t - 0.37, 2) + 0.912;
	else if (t >= 0.458 && t < 0.5901) 
		f = -27.7 * pow(t - 0.5, 2) + 0.5;
	else if (t >= 0.5901 && t < 0.7523)
		f = -23.7 * pow(t - 0.66, 2) + 0.391;
	return f;
	
/*
	if ((t > 0.213793 && t < 0.5) || t > 0.713793)
		return 0.1;
	else if (t <= 0.213793)
		return (-35 * pow(t - 0.1069, 2) + 0.5) * (x * (2 - x / R) / R);
	else
		return (-35 * pow(t - 0.5 - 0.1069, 2) + 0.5) * (x * (2 - x / R) / R);
		*/
}

static inline double w_0_derivative_t(double t)
{
	double f = 0;

	if (t < 0.2598 || t >= 0.7523)
		f = 0;
	else if (t >= 0.2598 && t < 0.458) 
		f = -59.5 * 2 * (t - 0.37);
	else if (t >= 0.458 && t < 0.5901) 
		f = -27.7 * 2 * (t - 0.5);
	else if (t >= 0.5901 && t < 0.7523)
		f = -23.7 * 2 * (t - 0.66);
	return f; 
	
/*
	if ((t > 0.213793 && t < 0.5) || t > 0.713793)
		return 0.1;
	else if (t <= 0.213793)
		return (-35 * pow(t - 0.1069, 2) + 0.5) * (x * (2 - x / R) / R);
	else
		return (-35 * pow(t - 0.5 - 0.1069, 2) + 0.5) * (x * (2 - x / R) / R);
		*/
}

static void calculate_mu(grid_node *grid, size_t n, size_t m)
{
	double u_n, u_s, u_ne, u_se, w_e, w_w, w_ne, w_nw, dz, dx, shear_rate;	
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(u_n, u_s, u_ne, u_se, w_e, w_w, w_ne, w_nw, dz, dx, shear_rate)
#endif
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 0; k < m - 1; ++k) {
			u_n = (k == m - 2 ? 0 : grid[i * m + k + 1].u);
			u_s = (k == 0 ? grid[i * m + k].u : grid[i * m + k - 1].u);
			u_ne = (k == m - 2 ? 0 : grid[(i + 1) * m + k + 1].u);
			u_se = (k == 0 ? grid[(i + 1) * m + k].u : grid[(i + 1) * m + k - 1].u);
			w_e = (i == n - 2 ? -grid[i * m + k].w : grid[(i + 1) * m + k].w);
			w_w = (i == 0 ? -grid[i * m + k].w : grid[(i - 1) * m + k].w);
			w_ne = (i == n - 2 ? -grid[i * m + k + 1].w : grid[(i + 1) * m + k + 1].w);
			w_nw = (i == 0 ? -grid[i * m + k + 1].w : grid[(i - 1) * m + k + 1].w);
			dz = grid[i * m + k].dz;
			dx = grid[i * m + k].dx;

			shear_rate = 
					((u_n - u_s + u_ne - u_se) / dz + 
					(w_e - w_w + w_ne - w_nw) / dx) / 8;
		// CARREAU 
		//	grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF);
			//if (fabs(shear_rate) <= EPS)
				//grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF);
			//else
				//grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF) *
					//pow(1 + pow(B * shear_rate, 2), (N - 1) / 2);
		//	CARREAU-YASUDA
			if (fabs(shear_rate) <= EPS)
				grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF);
			else
				grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF) *
					pow(1 + pow(B * fabs(shear_rate), A), (N - 1) / A);
		//	CROSS
		//	if (shear_rate <= EPS) 
			//	grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF);
		//	else
			//	grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF) /
				//	pow(1 + pow(B * shear_rate, M), A);
		}
	}
}

static inline void calculate_velocity_x_a_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double dx_half, dx_w, c1, d1, mu_w;	
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dx_half, dx_w, c1, d1, mu_w)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 1; ++i) {
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dx_w = grid[(i - 1) * m + k].dx;
			d1 = (grid[i * m + k].u + grid[(i - 1) * m + k].u) / 2;
			c1 = (-d1 >= 0 ? 0 : 1);
			mu_w = grid[(i - 1) * m + k].mu;

			ti[k].a[i] = delta_t * (-2 * mu_w / (density * dx_w * dx_half) - c1 * d1 / dx_half);
		}
	}
}

static inline void calculate_velocity_x_b_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double d1, d2, c1, c2, mu, mu_w, dx_w, dx, dx_half;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(d1, d2, c1, c2, mu, mu_w, dx_w, dx, dx_half)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 1; ++i) {
			d1 = (grid[i * m + k].u + grid[(i - 1) * m + k].u) / 2;
			d2 = (grid[i * m + k].u + grid[(i + 1) * m + k].u) / 2;
			c1 = (-d1 >= 0 ? 0 : 1);
			c2 = (d2 >= 0 ? 0 : 1);
			mu = grid[i * m + k].mu;
			mu_w = grid[(i - 1) * m + k].mu;
			dx_w = grid[(i - 1) * m + k].dx;
			dx = grid[i * m + k].dx;
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;

			ti[k].b[i] = 1 + delta_t * (2 * mu / (density * dx_half * dx) +
				2 * mu_w / (density * dx_half * dx_w) +
				((1 - c2) * d2 - (1 - c1) * d1) / dx_half);
		}
	}
}

static inline void calculate_velocity_x_c_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double dx_half, dx, c2, d2, mu;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dx_half, dx, c2, d2, mu)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 1; ++i) {
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dx = grid[i * m + k].dx;
			d2 = (grid[i * m + k].u + grid[(i + 1) * m + k].u) / 2;
			c2 = (d2 >= 0 ? 0 : 1);
			mu = grid[i * m + k].mu;
			
			ti[k].c[i] =
				delta_t * (-2 * mu / (density * dx_half * dx) + c2 * d2 / dx_half);
		}
	}
}

static inline void calculate_velocity_x_d_half_1(grid_node *grid, grid_node *grid_prev, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double dx_half, dz, dz_half_plus, d3, d4, c3, c4, mu_avg_n, mu_avg, u, u_prev, u_n, w_n, w_nw, w, w_w, p_w, p, dz_half_minus, u_s;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dx_half, dz, dz_half_plus, d3, d4, c3, c4, mu_avg_n, mu_avg, u, u_prev, u_n, w_n, w_nw, w, w_w, p_w, p, dz_half_minus, u_s)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 1; ++i) {
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dz = grid[i * m + k].dz;
			dz_half_plus = (k == m - 2 ? 0 : (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2);
			dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
			d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
			d4 = (grid[(i - 1) * m + k + 1].w + grid[i * m + k + 1].w) / 2;
			c3 = (-d3 >= 0 ? 0 : 1);
			c4 = (d4 >= 0 ? 0 : 1);
			mu_avg_n = (k == m - 2 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
					(grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4);
			mu_avg = (k == 0 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
					(grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4);
			u = grid[i * m + k].u;
			u_prev = grid_prev[i * m + k].u;
			u_n = (k == m - 2 ? 0 : grid[i * m + k + 1].u);
			u_s = (k == 0 ? 0 : grid[i * m + k - 1].u);
			w_n = grid[i * m + k + 1].w;
			w_nw = grid[(i - 1) * m + k + 1].w;
			w = grid[i * m + k].w;
			w_w = grid[(i - 1) * m + k].w;
			p_w = grid[(i - 1) * m + k].p;
			p = grid[i * m + k].p;

			if (k == 0) {
				ti[k].d[i] = u_prev + delta_t * (mu_avg_n * (u_n - u) / (dz * dz_half_plus * density) -
					(c4 * d4 * u_n + ((1 - c4) * d4 - d3) * u) / dz -
					(p - p_w) / (dx_half * density) + (1 / (density * dz * dx_half)) *
					(mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)));
			} else if (k == m - 2) {
				ti[k].d[i] = u_prev + delta_t * ((mu_avg_n * (-u) / dz - mu_avg * (u - u_s) / dz_half_minus) / (dz * density) +
					((1 - c3) * d3 * u + c3 * d3 * u_s) / dz -
					(p - p_w) / (dx_half * density) + (1 / (density * dz * dx_half)) *
					(mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)));
			} else {
				ti[k].d[i] = u_prev + delta_t * ((mu_avg_n * (u_n - u) / dz_half_plus - mu_avg * (u - u_s) / dz_half_minus) / (dz * density) -
					(c4 * d4 * u_n + ((1 - c4) * d4 - (1 - c3) * d3) * u - c3 * d3 * u_s) / dz -
					(p - p_w) / (dx_half * density) + (1 / (density * dz * dx_half)) *
					(mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)));
			}
		}
	}
}

static inline void calculate_velocity_x_a_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double dz_half, dz, d3, c3, mu_avg;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dz_half, dz, d3, c3, mu_avg)
#endif
	for (size_t i = 1; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			dz_half = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			dz = grid[i * m + k].dz;
			d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
			c3 = (-d3 >= 0 ? 0 : 1);
			mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;

			ti[i].a[k] = delta_t * (-mu_avg / (density * dz * dz_half) - c3 * d3 / dz);
		}
	}
}

static inline void calculate_velocity_x_b_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double dz_half_minus, dz_half_plus, dz, d3, d4, c3, c4, mu_avg, mu_avg_n;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dz_half_minus, dz_half_plus, dz, d3, d4, c3, c4, mu_avg, mu_avg_n)
#endif
	for (size_t i = 1; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
			dz = grid[i * m + k].dz;
			d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
			d4 = (grid[(i - 1) * m + k + 1].w + grid[i * m + k + 1].w) / 2;
			c3 = (-d3 >= 0 ? 0 : 1);
			c4 = (d4 >= 0 ? 0 : 1);
			mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4;
			mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;

			ti[i].b[k] = 1 + delta_t * (mu_avg_n / (density * dz * dz_half_plus) +
				mu_avg / (density * dz * dz_half_minus) +
				((1 - c4) * d4 - (1 - c3) * d3) / dz);
		}
	}
}

static inline void calculate_velocity_x_c_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double d4, c4, dz_half_plus, dz, mu_avg_n;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(d4, c4, dz_half_plus, dz, mu_avg_n)
#endif
	for (size_t i = 1; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			d4 = (grid[(i - 1) * m + k + 1].w + grid[i * m + k + 1].w) / 2;
			c4 = (d4 >= 0 ? 0 : 1);
			dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
			dz = grid[i * m + k].dz;
			mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4;

			ti[i].c[k] = delta_t * (-mu_avg_n / (density * dz * dz_half_plus) + c4 * d4 / dz);
		}
	}
}

static inline void calculate_velocity_x_d_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, tma_info *ti_half_2, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double d1, d2, c1, c2, dx_half, dx_w, dx, dz, u, u_w, u_e, u_prev, mu, mu_w, mu_avg_n, mu_avg, p, p_w, w, w_n, w_nw, w_w;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(d1, d2, c1, c2, dx_half, dx_w, dx, dz, u, u_w, u_e, u_prev, mu, mu_w, mu_avg_n, mu_avg, p, p_w, w, w_n, w_nw, w_w)
#endif
	for (size_t i = 1; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			d1 = (grid[i * m + k].u + grid[(i - 1) * m + k].u) / 2;
			d2 = (grid[i * m + k].u + grid[(i + 1) * m + k].u) / 2;
			c1 = (-d1 >= 0 ? 0 : 1);
			c2 = (d2 >= 0 ? 0 : 1);
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dx_w = grid[(i - 1) * m + k].dx;
			dx = grid[i * m + k].dx;
			dz = grid[i * m + k].dz;
			u = ti_half_1[k].y[i];
			u_w = ti_half_1[k].y[i - 1];
			u_e = ti_half_1[k].y[i + 1];
			u_prev = grid_prev[i * m + k].u;
			mu = grid[i * m + k].mu;
			mu_w = grid[(i - 1) * m + k].mu;
			mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4;
			mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;
			p = grid[i * m + k].p;
			p_w = grid[(i - 1) * m + k].p;
			w = grid[i * m + k].w;
			w_n = grid[i * m + k + 1].w;
			w_nw = grid[(i - 1) * m + k + 1].w;
			w_w = grid[(i - 1) * m + k].w;

			ti_half_2[i].d[k] = u_prev + delta_t * ((2 / (density * dx_half)) *
				(mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) -
				(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) /
				dx_half - (p - p_w) / (dx_half * density) +	(1 / (density * dz)) *
				(mu_avg_n * (w_n - w_nw) / dx_half - mu_avg * (w - w_w) / dx_half));
		}
	}
}

static inline double calculate_velocity_x_xi_1_half_2(grid_node *grid, size_t i, size_t m)
{
	// k == 0
	double d3 = (grid[i * m].w + grid[(i - 1) * m].w) / 2;
	double d4 = (grid[(i - 1) * m + 1].w + grid[i * m + 1].w) / 2;
	double c4 = (d4 >= 0 ? 0 : 1);
	double mu_avg_n = (grid[i * m].mu + grid[(i - 1) * m].mu +
			grid[i * m + 1].mu + grid[(i - 1) * m + 1].mu) / 4;
	double dz = grid[i * m].dz;
	double dz_half_plus = (grid[i * m].dz + grid[i * m + 1].dz) / 2;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return delta_t * (mu_avg_n / (density * dz * dz_half_plus) - c4 * d4 / dz) /
			(1 + delta_t * (mu_avg_n / (density * dz * dz_half_plus) +
			 ((1 - c4) * d4 - d3) / dz));
}

static inline double calculate_velocity_x_xi_2_half_2(grid_node *grid, size_t i, ssize_t m)
{
	// k == m - 2
	double d3 = (grid[i * m + m - 2].w + grid[(i - 1) * m + m - 2].w) / 2;
	double c3 = (-d3 >= 0 ? 0 : 1);
	double dz = grid[i * m + m - 2].dz;
	double dz_half_minus = (grid[i * m + m - 2].dz + grid[i * m + m - 3].dz) / 2;
	double mu_avg_n = (grid[i * m + m - 2].mu + grid[(i - 1) * m + m - 2].mu) / 2;
	double mu_avg = (grid[i * m + m - 2].mu + grid[(i - 1) * m + m - 2].mu +
			grid[i * m + m - 3].mu + grid[(i - 1) * m + m - 3].mu) / 4;
	double density = get_blood_density();
	double delta_t = get_delta_t();

	return delta_t * (c3 * d3 / dz + mu_avg / (density * dz_half_minus * dz)) / 
		(1 + delta_t * (mu_avg_n / (density * dz * dz) +
		 mu_avg / (density * dz * dz_half_minus) - ((1 - c3) * d3) / dz));
}

static inline double calculate_velocity_x_mu_1_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, size_t i, size_t m)
{
	// k == 0
	double d1 = (grid[i * m].u + grid[(i - 1) * m].u) / 2;
	double d2 = (grid[i * m].u + grid[(i + 1) * m].u) / 2;
	double d3 = (grid[i * m].w + grid[(i - 1) * m].w) / 2;
	double d4 = (grid[(i - 1) * m + 1].w + grid[i * m + 1].w) / 2;
	double c1 = (-d1 >= 0 ? 0 : 1);
	double c2 = (d2 >= 0 ? 0 : 1);
	double c4 = (d4 >= 0 ? 0 : 1);
	double dx = grid[i * m].dx;
	double dx_w = grid[(i - 1) * m].dx;
	double dx_half = (dx + dx_w) / 2;
	double dz_half_plus = (grid[i * m].dz + grid[i * m + 1].dz) / 2;
	double dz = grid[i * m].dz;
	double u_prev = grid_prev[i * m].u;
	double u_e = ti_half_1[0].y[i + 1];
	double u_w = ti_half_1[0].y[i - 1];
	double u = ti_half_1[0].y[i];
	double mu = grid[i * m].mu;
	double mu_avg = (grid[i * m].mu + grid[(i - 1) * m].mu) / 2;
	double mu_avg_n = (grid[i * m].mu + grid[(i - 1) * m].mu +
			grid[i * m + 1].mu + grid[(i - 1) * m + 1].mu) / 4;
	double mu_w = grid[(i - 1) * m].mu;
	double p = grid[i * m].p;
	double p_w = grid[(i - 1) * m].p;
	double w_n = grid[i * m + 1].w;
	double w_nw = grid[(i - 1) * m + 1].w;
	double w = grid[i * m].w;
	double w_w = grid[(i - 1) * m].w;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return (u_prev + delta_t * ((2 / (density * dx_half)) *
			(mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) -
			(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) / dx_half -
			(1 / (density * dx_half)) * (p - p_w) + (1 / (density * dz)) *
			(mu_avg_n * (w_n - w_nw) / dx_half - mu_avg * (w - w_w) / dx_half))) / 
			(1 + delta_t * (mu_avg_n / (density * dz * dz_half_plus) +
			 ((1 - c4) * d4 - d3) / dz));
}

static inline double calculate_velocity_x_mu_2_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, size_t i, ssize_t m)
{
	// k == m - 2
	double d1 = (grid[i * m + m - 2].u + grid[(i - 1) * m + m - 2].u) / 2;
	double d2 = (grid[i * m + m - 2].u + grid[(i + 1) * m + m - 2].u) / 2;
	double d3 = (grid[i * m + m - 2].w + grid[(i - 1) * m + m - 2].w) / 2;
	double c1 = (-d1 >= 0 ? 0 : 1);
	double c2 = (d2 >= 0 ? 0 : 1);
	double c3 = (-d3 >= 0 ? 0 : 1);
	double dz = grid[i * m + m - 2].dz;
	double dz_half_minus = (grid[i * m + m - 2].dz + grid[i * m + m - 3].dz) / 2;
	double dx_half = (grid[i * m + m - 2].dx + grid[(i - 1) * m + m - 3].dx) / 2;
	double dx_w = grid[(i - 1) * m + m - 2].dx;
	double dx = grid[i * m + m - 2].dx;
	double mu = grid[i * m + m - 2].mu;
	double mu_avg_n = (grid[i * m + m - 2].mu + grid[(i - 1) * m + m - 2].mu) / 2;
	double mu_avg = (grid[i * m + m - 2].mu + grid[(i - 1) * m + m - 2].mu +
			grid[i * m + m - 3].mu + grid[(i - 1) * m + m - 3].mu) / 4;
	double mu_w = grid[(i - 1) * m + m - 2].mu;
	double u_e = ti_half_1[m - 2].y[i + 1];
	double u_w = ti_half_1[m - 2].y[i - 1];
	double u_prev = grid_prev[i * m + m - 2].u;
	double u = ti_half_1[m - 2].y[i];
	double p = grid[i * m + m - 2].p;
	double p_w = grid[(i - 1) * m + m - 2].p;
	double w_n = grid[i * m + m - 1].w;
	double w_nw = grid[(i - 1) * m + m - 1].w;
	double w = grid[i * m + m - 2].w;
	double w_w = grid[(i - 1) * m + m - 2].w;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return (u_prev + delta_t * ((2 / (density * dx_half)) * (mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) -
			(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) / dx_half -
			(p - p_w) / (density * dx_half) + 1 / (density * dz) * 
			(mu_avg_n * (w_n - w_nw) / dx_half - mu_avg * (w - w_w) / dx_half))) /
			(1 + delta_t * (mu_avg_n / (density * dz * dz) +
			 mu_avg / (density * dz * dz_half_minus) -
			 (1 - c3) * d3 / dz));
}

static void calculate_velocity_x(grid_node *grid_prev, grid_node *grid_cur, grid_node *grid_last_layer, size_t n, size_t m)
{
	prs_info pi;
	pi.n = n - 1;
	pi.m = m - 1;
	pi.start_half_1 = 0;
	pi.start_half_2 = 1;
	pi.ti_half_1 = calloc(pi.m, sizeof(tma_info));

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t k = 0; k < pi.m; ++k) {
		pi.ti_half_1[k].n = pi.n;
		pi.ti_half_1[k].a = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].b = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].c = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].d = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].y = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].xi_1 = 0;
		pi.ti_half_1[k].xi_2 = 0;
		pi.ti_half_1[k].mu_1 = 0;
		pi.ti_half_1[k].mu_2 = 0;
	}

	calculate_velocity_x_a_half_1(grid_prev, pi.ti_half_1, n, m);
	/*fprintf(stderr, "VELOCITY_X_A_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_X_B_HALF_1:\n");
	*/
	calculate_velocity_x_b_half_1(grid_prev, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_X_C_HALF_1:\n");
	*/
	calculate_velocity_x_c_half_1(grid_prev, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_X_D_HALF_1:\n");
	*/
	calculate_velocity_x_d_half_1(grid_prev, grid_last_layer, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].d[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	prs_half_1(&pi);
	/*
	fprintf(stderr, "VELOCITY_X_Y_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].y[i]);
		}
		fprintf(stderr, "\n");
	}
	*/

	pi.ti_half_2 = calloc(pi.n, sizeof(tma_info));

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 1; i < pi.n; ++i) {
		pi.ti_half_2[i].n = pi.m - 1;
		pi.ti_half_2[i].a = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].b = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].c = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].d = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].y = calloc(m, sizeof(double));
		pi.ti_half_2[i].xi_1 = calculate_velocity_x_xi_1_half_2(grid_prev, i, m);
		pi.ti_half_2[i].xi_2 = calculate_velocity_x_xi_2_half_2(grid_prev, i, m);
		pi.ti_half_2[i].mu_1 = calculate_velocity_x_mu_1_half_2(grid_prev, grid_last_layer, pi.ti_half_1, i, m);
		pi.ti_half_2[i].mu_2 = calculate_velocity_x_mu_2_half_2(grid_prev, grid_last_layer, pi.ti_half_1, i, m);
		//fprintf(stderr, "XI_1 = %lf, XI_2 = %lf, MU_1 = %lf, MU_2 = %lf\n", pi.ti_half_2[i].xi_1, pi.ti_half_2[i].xi_2, pi.ti_half_2[i].mu_1, pi.ti_half_2[i].mu_2);
	}

	calculate_velocity_x_a_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_X_A_HALF_2:\n");
	for (size_t k = 1; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	calculate_velocity_x_b_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_X_B_HALF_2:\n");
	for (size_t k = 1; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	calculate_velocity_x_c_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_X_C_HALF_2:\n");
	for (size_t k = 1; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	calculate_velocity_x_d_half_2(grid_prev, grid_last_layer, pi.ti_half_1, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_X_D_HALF_2:\n");
	for (size_t k = 1; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].d[i]);
		}
		fprintf(stderr, "\n");
	}
	*/

	prs_half_2(&pi);

	/*
	fprintf(stderr, "VELOCITY_X_Y_HALF_2:\n");
	for (size_t k = 1; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].y[i]);
		}
		fprintf(stderr, "\n");
	}
	*/

	double *tmp = calloc(n * m, sizeof(double));
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 1; i < pi.n; ++i) {
		memcpy(&tmp[i * m], pi.ti_half_2[i].y, m * sizeof(double));
	}

	grid_fill(grid_cur, tmp, GRID_OFFSET_U);

	free(tmp);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 1; i < pi.n; ++i) {
		free(pi.ti_half_2[i].a);
		free(pi.ti_half_2[i].b);
		free(pi.ti_half_2[i].c);
		free(pi.ti_half_2[i].d);
		free(pi.ti_half_2[i].y);
	}
	free(pi.ti_half_2);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.m; ++i) {
		free(pi.ti_half_1[i].a);
		free(pi.ti_half_1[i].b);
		free(pi.ti_half_1[i].c);
		free(pi.ti_half_1[i].d);
		free(pi.ti_half_1[i].y);
	}
	free(pi.ti_half_1);
}

static inline void calculate_velocity_z_a_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double mu_avg, dx_half, dx, d1, c1;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(mu_avg, dx_half, dx, d1, c1)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			mu_avg = (k == 0 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
									(grid[i * m + k].mu + grid[(i - 1) * m + k].mu + 
									 grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4);
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dx = grid[i * m + k].dx;
			d1 = (k == 0 ? -grid[i * m + k].u :
					-(grid[i * m + k].u + grid[i * m + k - 1].u) / 2);
			c1 = (d1 >= 0 ? 0 : 1);

			ti[k].a[i] = delta_t * (-mu_avg / (density * dx_half * dx) + c1 * d1 / dx);
		}
	}
}

static inline void calculate_velocity_z_b_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double d1, d2, c1, c2, mu_avg, mu_avg_e, dx, dx_half, dx_half_plus;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(d1, d2, c1, c2, mu_avg, mu_avg_e, dx, dx_half, dx_half_plus)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			d1 = (k == 0 ? -grid[i * m + k].u :
					-(grid[i * m + k].u + grid[i * m + k - 1].u) / 2);
			d2 = (k == 0 ? grid[(i + 1) * m + k].u :
					(grid[(i + 1) * m + k].u + grid[(i + 1) * m + k - 1].u) / 2);
			c1 = (d1 >= 0 ? 0 : 1);
			c2 = (d2 >= 0 ? 0 : 1);
			mu_avg = (k == 0 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
									(grid[i * m + k].mu + grid[(i - 1) * m + k].mu + 
									 grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4);
			mu_avg_e = (k == 0 ? (grid[i * m + k].mu + grid[(i + 1) * m + k].mu) / 2 :
									(grid[i * m + k].mu + grid[(i + 1) * m + k].mu + 
									 grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4);
			dx = grid[i * m + k].dx;
			dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;

			ti[k].b[i] = 1 + delta_t * (mu_avg_e / (density * dx * dx_half_plus) + 
				mu_avg / (density * dx * dx_half) + ((1 - c2) * d2 + (1 - c1) * d1) / dx);
		}
	}
}

static inline void calculate_velocity_z_c_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double mu_avg_e, dx_half_plus, dx, d2, c2;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(mu_avg_e, dx_half_plus, dx, d2, c2)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			mu_avg_e = (k == 0 ? (grid[i * m + k].mu + grid[(i + 1) * m + k].mu) / 2 :
									(grid[i * m + k].mu + grid[(i + 1) * m + k].mu + 
									 grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4);
			dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
			dx = grid[i * m + k].dx;
			d2 = (k == 0 ? grid[(i + 1) * m + k].u :
					(grid[(i + 1) * m + k].u + grid[(i + 1) * m + k - 1].u) / 2);
			c2 = (d2 >= 0 ? 0 : 1);
			
			ti[k].c[i] = delta_t * (-mu_avg_e / (density * dx * dx_half_plus) + c2 * d2 / dx);
		}
	}
}

static inline void calculate_velocity_z_d_half_1(grid_node *grid, grid_node *grid_prev, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double w_prev, mu, w, w_n, dz, d3, c3, dz_half_minus, dz_s, mu_s, w_s, d4, c4, p, p_s, dx, u, u_e, u_s, u_se, mu_avg, mu_avg_e;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(w_prev, mu, w, w_n, dz, d3, c3, dz_half_minus, dz_s, mu_s, w_s, d4, c4, p, p_s, dx, u, u_e, u_s, u_se, mu_avg, mu_avg_e)
#endif
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			w_prev = grid_prev[i * m + k].w;
			mu = grid[i * m + k].mu;
			w = grid[i * m + k].w;
			w_n = grid[i * m + k + 1].w;
			dz = grid[i * m + k].dz;
			d3 = (k == 0 ? -(grid[i * m + k].w + grid[i * m + k + 1].w) / 2 :
					-(grid[i * m + k].w + grid[i * m + k - 1].w) / 2);
			c3 = (d3 >= 0 ? 0 : 1);

			if (k == 0) {
				ti[k].d[i] = w_prev + delta_t * (4 * mu * (w_n - w) / (density * dz * dz) -
					(2 * c3 - 1) * d3 * (w_n - w) / dz);
			} else {
				dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
				dz_s = grid[i * m + k - 1].dz;
				mu_s = grid[i * m + k - 1].mu;
				w_s = grid[i * m + k - 1].w;
				d4 = (grid[i * m + k].w + grid[i * m + k + 1].w) / 2;
				c4 = (d4 >= 0 ? 0 : 1);
				p = grid[i * m + k].p;
				p_s = grid[i * m + k - 1].p;
				dx = grid[i * m + k].dx;
				u = grid[i * m + k].u;
				u_e = grid[(i + 1) * m + k].u;
				u_s = grid[i * m + k - 1].u;
				u_se = grid[(i + 1) * m + k - 1].u;
				mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu + 
								grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;
				mu_avg_e = (grid[i * m + k].mu + grid[(i + 1) * m + k].mu + 
								grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4;
				//if (i == 1 || i == n - 3) {
					//fprintf(stderr, "%zu: w_prev = %.18lf, mu = %.18lf, w = %.18lf, w_n = %.18lf, dz = %.18lf, d3 = %.18lf, c3 = %.18lf, dz_half_minus = %.18lf, dz_s = %.18lf, mu_s = %.18lf, w_s = %.18lf, d4 = %.18lf, c4 = %.18lf, p = %.18lf, p_s = %.18lf, dx = %.18lf, u = %.18lf, u_e = %.18lf, u_s = %.18lf, u_se = %.18lf, mu_avg = %.18lf, mu_avg_e = %.18lf\n", i, w_prev, mu, w, w_n, dz, d3, c3, dz_half_minus, dz_s, mu_s, w_s, d4, c4, p, p_s, dx, u, u_e, u_s, u_se, mu_avg, mu_avg_e);
				//}

				ti[k].d[i] = w_prev + delta_t * ((2 * mu * (w_n - w) / dz - mu_s * (w - w_s) / dz_s) / (density * dz_half_minus) -
					(c4 * d4 * w_n + ((1 - c4) * d4 + (1 - c3) * d3) * w + c3 * d3 * w_s) / dz_half_minus -
					(p - p_s) / (density * dz_half_minus) + (mu_avg_e * (u_e - u_se) - mu_avg * (u - u_s)) / (density * dx * dz_half_minus));
			}
		}
	}
}

static inline double calculate_velocity_z_xi_1_half_1(grid_node *grid, size_t k, size_t m)
{
	size_t i = 0;
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double mu_avg_e = (k == 0 ? (grid[i * m + k].mu + grid[(i + 1) * m + k].mu) / 2 :
							(grid[i * m + k].mu + grid[(i + 1) * m + k].mu + 
							 grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4);
	double mu_avg = (k == 0 ? grid[i * m + k].mu :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2);
	double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
	double dx = grid[i * m + k].dx;
	double d2 = (k == 0 ? grid[(i + 1) * m + k].u :
			(grid[(i + 1) * m + k].u + grid[(i + 1) * m + k - 1].u) / 2);
	double c2 = (d2 >= 0 ? 0 : 1);

	return delta_t * (mu_avg_e / (density * dx_half_plus * dx) - c2 * d2 / dx) / 
		(1 + delta_t * (mu_avg_e / (density * dx_half_plus * dx) + 2 * mu_avg / (density * dx * dx) +
		(1 - c2) * d2 / dx));
}

static inline double calculate_velocity_z_xi_2_half_1(grid_node *grid, size_t k, ssize_t n, ssize_t m)
{
	ssize_t i = n - 2;
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double mu_avg_e = (k == 0 ? grid[i * m + k].mu :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2);
	double mu_avg = (k == 0 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu +
							 grid[(i - 1) * m + k].mu + grid[(i - 1) * m + k - 1].mu) / 4);
	double dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
	double dx = grid[i * m + k].dx;
	double d1 = (k == 0 ? -grid[i * m + k].u :
			-(grid[i * m + k].u + grid[i * m + k - 1].u) / 2);
	double c1 = (d1 >= 0 ? 0 : 1);

	return delta_t * (mu_avg / (density * dx * dx_half) - c1 * d1 / dx) /
			(1 + delta_t * (2 * mu_avg_e / (density * dx * dx) + mu_avg / (density * dx_half * dx) +
			(1 - c1) * d1 / dx));
}

static inline double calculate_velocity_z_mu_1_half_1(grid_node *grid, grid_node *grid_prev, size_t k, size_t m)
{
	size_t i = 0;
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double w_prev = grid_prev[i * m + k].w;
	double mu = grid[i * m + k].mu;
	double mu_s = (k == 0 ? 0 : grid[i * m + k - 1].mu);
	double mu_avg_e = (k == 0 ? (grid[i * m + k].mu + grid[(i + 1) * m + k].mu) / 2 :
							(grid[i * m + k].mu + grid[(i + 1) * m + k].mu + 
							 grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4);
	double mu_avg = (k == 0 ? grid[i * m + k].mu :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2);
	double w = grid[i * m + k].w;
	double w_n = grid[i * m + k + 1].w;
	double w_s = (k == 0 ? 0 : grid[i * m + k - 1].w);
	double dz = grid[i * m + k].dz;
	double dz_s = (k == 0 ? 0 : grid[i * m + k - 1].dz);
	double dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
	double d2 = (k == 0 ? grid[(i + 1) * m + k].u :
			(grid[(i + 1) * m + k].u + grid[(i + 1) * m + k - 1].u) / 2);
	double d3 = (k == 0 ? -(grid[i * m + k].w + grid[i * m + k + 1].w) / 2 :
			-(grid[i * m + k].w + grid[i * m + k - 1].w) / 2);
	double d4 = (k == 0 ? -d3 : (grid[i * m + k].w + grid[i * m + k + 1].w) / 2);
	double c2 = (d2 >= 0 ? 0 : 1);
	double c3 = (d3 >= 0 ? 0 : 1);
	double c4 = (d4 >= 0 ? 0 : 1);
	double p = grid[i * m + k].p;
	double p_s = (k == 0 ? 0 : grid[i * m + k - 1].p);
	double dx = grid[i * m + k].dx;
	double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
	double u_e = grid[(i + 1) * m + k].u;
	double u_se = (k == 0 ? 0 : grid[(i + 1) * m + k - 1].u);

	return (w_prev + delta_t * (k == 0 ? (4 * mu_avg * (w_n - w) / (density * dz * dz) -
			(2 * c3 - 1) * (w_n - w) * d3 / dz) :
			((2 / (density * dz_half_minus)) * (mu * (w_n - w) / dz - mu_s * (w - w_s) / dz_s) -
			(c4 * d4 * w_n + ((1 - c4) * d4 + (1 - c3) * d3) * w + c3 * d3 * w_s) / dz_half_minus -
			(1 / (density * dz_half_minus)) * (p - p_s) + (1 / (density * dx * dz_half_minus)) * mu_avg_e * (u_e - u_se)))) / 
			(1 + delta_t * (mu_avg_e / (density * dx_half_plus * dx) + 2 * mu_avg / (density * dx * dx) +
			(1 - c2) * d2 / dx));
}

static inline double calculate_velocity_z_mu_2_half_1(grid_node *grid, grid_node *grid_prev, size_t k, ssize_t n, ssize_t m)
{
	ssize_t i = n - 2;
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double mu_avg_e = (k == 0 ? grid[i * m + k].mu :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2);
	double mu_avg = (k == 0 ? (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2 :
							(grid[i * m + k].mu + grid[i * m + k - 1].mu +
							 grid[(i - 1) * m + k].mu + grid[(i - 1) * m + k - 1].mu) / 4);
	double dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
	double dx = grid[i * m + k].dx;
	double d1 = (k == 0 ? -grid[i * m + k].u :
			-(grid[i * m + k].u + grid[i * m + k - 1].u) / 2);
	double d3 = (k == 0 ? -(grid[i * m + k].w + grid[i * m + k + 1].w) / 2 :
			-(grid[i * m + k].w + grid[i * m + k - 1].w) / 2);
	double d4 = (k == 0 ? -d3 : (grid[i * m + k].w + grid[i * m + k + 1].w) / 2);
	double c1 = (d1 >= 0 ? 0 : 1);
	double c3 = (d3 >= 0 ? 0 : 1);
	double c4 = (d4 >= 0 ? 0 : 1);
	double w_prev = grid_prev[i * m + k].w;
	double w = grid[i * m + k].w;
	double w_s = (k == 0 ? 0 : grid[i * m + k - 1].w);
	double w_n = grid[i * m + k + 1].w;
	double dz = grid[i * m + k].dz;
	double dz_s = (k == 0 ? 0 : grid[i * m + k - 1].dz);
	double dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
	double mu = grid[i * m + k].mu;
	double mu_s = (k == 0 ? 0 : grid[i * m + k - 1].mu);
	double p = grid[i * m + k].p;
	double p_s = (k == 0 ? 0 : grid[i * m + k - 1].p);
	double u = grid[i * m + k].u;
	double u_s = (k == 0 ? 0 : grid[i * m + k - 1].u);

	return (k == 0 ? (w_prev + delta_t * (4 * mu * (w_n - w) / (density * dz * dz) -
				(2 * c3 - 1) * (w_n - w) * d3 / dz)) :
			(w_prev + delta_t * ((2 / (density * dz_half_minus)) * (mu * (w_n - w) / dz - mu_s * (w - w_s) / dz_s) -
			  (c4 * d4 * w_n + ((1 - c4) * d4 + (1 - c3) * d3) * w + c3 * d3 * w_s) / dz_half_minus - 
			  (p - p_s) / (dz_half_minus * density) - mu_avg * (u - u_s) / (density * dz_half_minus * dx)))) /
			(1 + delta_t * (2 * mu_avg_e / (density * dx * dx) + mu_avg / (density * dx_half * dx) +
			(1 - c1) * d1 / dx));
}

static inline void calculate_velocity_z_a_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double dz_half, dz_s, d3, c3, mu_s;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dz_half, dz_s, d3, c3, mu_s)
#endif
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 1; ++k) {
			dz_half = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			dz_s = grid[i * m + k - 1].dz;
			d3 = -(grid[i * m + k].w + grid[i * m + k - 1].w) / 2;
			c3 = (d3 >= 0 ? 0 : 1);
			mu_s = grid[i * m + k - 1].mu;

			ti[i].a[k] = delta_t * (-2 * mu_s / (density * dz_half * dz_s) + c3 * d3 / dz_half);
		}
	}
}

static inline void calculate_velocity_z_b_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double dz_half_minus, dz, dz_s, d3, d4, c3, c4, mu_s, mu;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(dz_half_minus, dz, dz_s, d3, d4, c3, c4, mu_s, mu)
#endif
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 1; ++k) {
			dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			dz = grid[i * m + k].dz;
			dz_s = grid[i * m + k - 1].dz;
			d3 = -(grid[i * m + k].w + grid[i * m + k - 1].w) / 2;
			d4 = (grid[i * m + k].w + grid[i * m + k + 1].w) / 2;
			c3 = (d3 >= 0 ? 0 : 1);
			c4 = (d4 >= 0 ? 0 : 1);
			mu_s = grid[i * m + k - 1].mu;
			mu = grid[i * m + k].mu;

			ti[i].b[k] = 1 + delta_t * (2 * mu / (density * dz * dz_half_minus) +
				2 * mu_s / (density * dz_s * dz_half_minus) +
				((1 - c4) * d4 + (1 - c3) * d3) / dz_half_minus);
		}
	}
}

static inline void calculate_velocity_z_c_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double d4, c4, dz_half_minus, dz, mu;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(d4, c4, dz_half_minus, dz, mu)
#endif
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 1; ++k) {
			d4 = (grid[i * m + k].w + grid[i * m + k + 1].w) / 2;
			c4 = (d4 >= 0 ? 0 : 1);
			dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			dz = grid[i * m + k].dz;
			mu = grid[i * m + k].mu;

			ti[i].c[k] = delta_t * (-2 * mu / (density * dz * dz_half_minus) + c4 * d4 / dz_half_minus);
		}
	}
}

static inline void calculate_velocity_z_d_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, tma_info *ti_half_2, size_t n, size_t m)
{
	double delta_t = get_delta_t();
	double density = get_blood_density();
	double w_prev, mu_avg, mu_avg_e, w, w_w, w_e, dx, dx_half_plus, dx_half, d1, d2, c1, c2, u, u_s, u_e, u_se, dz_half_minus, p, p_s;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) private(w_prev, mu_avg, mu_avg_e, w, w_w, w_e, dx, dx_half_plus, dx_half, d1, d2, c1, c2, u, u_s, u_e, u_se, dz_half_minus, p, p_s)
#endif
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 1; ++k) {
			w_prev = grid_prev[i * m + k].w;
			mu_avg = (i == 0 ? (grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2 :
					(grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4);
			mu_avg_e = (i == n - 2 ? (grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2 :
					(grid[i * m + k].mu + grid[(i + 1) * m + k].mu +
					grid[i * m + k - 1].mu + grid[(i + 1) * m + k - 1].mu) / 4);
			w = ti_half_1[k].y[i];
			w_w = (i == 0 ? 0 : ti_half_1[k].y[i - 1]);
			w_e = (i == n - 2 ? 0 : ti_half_1[k].y[i + 1]);
			dx = grid[i * m + k].dx;
			dx_half_plus = (i == n - 2 ? 0 : (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2);
			dx_half = (i == 0 ? 0 : (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2);
			d1 =	-(grid[i * m + k].u + grid[i * m + k - 1].u) / 2;
			d2 = (grid[(i + 1) * m + k].u + grid[(i + 1) * m + k - 1].u) / 2;
			c1 = (d1 >= 0 ? 0 : 1);
			c2 = (d2 >= 0 ? 0 : 1);
			u = grid[i * m + k].u;
			u_s = grid[i * m + k - 1].u;
			u_e = grid[(i + 1) * m + k].u;
			u_se = grid[(i + 1) * m + k - 1].u;
			dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			p = grid[i * m + k].p;
			p_s = grid[i * m + k - 1].p;

			if (i == 0) {
				ti_half_2[i].d[k] = w_prev + delta_t * ((1 / (density * dx)) * 
					(mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (2 * w) / dx) -
					(c2 * d2 * w_e + (1 - c2) * d2 * w) / dx - (p - p_s) / (density * dz_half_minus) +
					(1 / (density * dx * dz_half_minus)) * mu_avg_e * (u_e - u_se));
			} else if (i == n - 2) {
				mu_avg_e = (grid[i * m + k].mu + grid[i * m + k - 1].mu) / 2;

				ti_half_2[i].d[k] = w_prev + delta_t * ((1 / (density * dx)) * 
					(-mu_avg_e * (2 * w) / dx - mu_avg * (w - w_w) / dx_half) -
					(c1 * d1 * w_w + (1 - c1) * d1 * w) / dx -
					(1 / (density * dz_half_minus)) * (p - p_s) - (mu_avg / (density * dx * dz_half_minus)) *
					(u - u_s));
			} else {
				ti_half_2[i].d[k] = w_prev + delta_t * ((1 / (density * dx)) *
					(mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (w - w_w) / dx_half) -
					(c2 * d2 * w_e + ((1 - c2) * d2 + (1 - c1) * d1) * w + c1 * d1 * w_w) / dx +
					(1 / (density * dx * dz_half_minus)) * (mu_avg_e * (u_e - u_se) - mu_avg * (u - u_s)) -
					(1 / density) * (p - p_s) / dz_half_minus);
			}
		}
	}
}

static inline double calculate_velocity_z_xi_1_half_2(grid_node *grid, size_t i, size_t m)
{
	size_t k = 0;
	double d3 = -(grid[i * m + k].w + grid[i * m + k + 1].w) / 2;
	double c3 = (d3 >= 0 ? 0 : 1);
	double mu = grid[i * m + k].mu;
	double dz = grid[i * m + k].dz;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz) /
			(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
}

static inline double calculate_velocity_z_mu_1_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, ssize_t i, ssize_t n, size_t m)
{
	ssize_t k = 0;
	double density = get_blood_density();
	double delta_t = get_delta_t();
	double w_prev = grid_prev[i * m + k].w;
	double dx = grid[i * m + k].dx;
	double dx_half_plus = (i == n - 2 ? 0 : (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2);
	double dx_half_minus = (i == 0 ? 0 : (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2);
	double mu_avg_e = (i == n - 2 ? grid[i * m + k].mu : (grid[i * m + k].mu + grid[(i + 1) * m + k].mu) / 2);
	double mu_avg = (i == 0 ? grid[i * m + k].mu : (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2);
	double w = ti_half_1[k].y[i];
	double w_e = (i == n - 2 ? 0 : ti_half_1[k].y[i + 1]);
	double w_w = (i == 0 ? 0 : ti_half_1[k].y[i - 1]);
	double d1 = -grid[i * m + k].u;
	double c1 = (d1 >= 0 ? 0 : 1);
	double d2 = grid[(i + 1) * m + k].u;
	double d3 = -(grid[i * m + k].w + grid[i * m + k + 1].w) / 2;
	double c2 = (d2 >= 0 ? 0 : 1);
	double c3 = (d3 >= 0 ? 0 : 1);
	double mu = grid[i * m + k].mu;
	double dz = grid[i * m + k].dz;

	if (i == 0) {
		return  (w_prev + delta_t * ((1 / (density * dx)) * (mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (2 * w) / dx) -
				(c2 * w_e + (1 - c2) * w) * d2 / dx)) /
				(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
	} else if (i == n - 2) {
		return  (w_prev + delta_t * ((1 / (density * dx)) * (-mu_avg_e * (2 * w) / dx - mu_avg * (w - w_w) / dx_half_minus) - 
				(c1 * d1 * w_w + (1 - c1) * d1 * w) / dx)) /
				(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
	} else {
		return (delta_t * ((1 / (density * dx)) * (mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (w - w_w) / dx_half_minus) -
				(c2 * d2 * w_e + ((1 - c2) * d2 + (1 - c1) * d1) * w + c1 * d1 * w_w) / dx) + w_prev) / 
				(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
	}
}

static inline double calculate_velocity_z_mu_2_half_2(grid_node *grid, size_t i, ssize_t m)
{
	size_t k = m - 1;

	return grid[i * m + k].w;
}

static void calculate_velocity_z(grid_node *grid_prev, grid_node *grid_cur, grid_node *grid_last_layer, size_t n, size_t m)
{
	prs_info pi;
	pi.n = n - 2;
	pi.m = m - 1;
	pi.start_half_1 = 0;
	pi.start_half_2 = 0;
	pi.ti_half_1 = calloc(pi.m, sizeof(tma_info));

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t k = 0; k < pi.m; ++k) {
		pi.ti_half_1[k].n = pi.n;
		pi.ti_half_1[k].a = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].b = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].c = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].d = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].y = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].xi_1 = calculate_velocity_z_xi_1_half_1(grid_prev, k, m);
		pi.ti_half_1[k].xi_2 = calculate_velocity_z_xi_2_half_1(grid_prev, k, n, m);
		pi.ti_half_1[k].mu_1 = calculate_velocity_z_mu_1_half_1(grid_prev, grid_last_layer, k, m);
		pi.ti_half_1[k].mu_2 = calculate_velocity_z_mu_2_half_1(grid_prev, grid_last_layer, k, n, m);
		//fprintf(stderr, "XI_1 = %lf, XI_2 = %lf, MU_1 = %lf, MU_2 = %lf\n", pi.ti_half_1[k].xi_1, pi.ti_half_1[k].xi_2, pi.ti_half_1[k].mu_1, pi.ti_half_1[k].mu_2);
		if (fabs(pi.ti_half_1[k].xi_1 - pi.ti_half_1[k].xi_2) > EPS || fabs(pi.ti_half_1[k].mu_1 - pi.ti_half_1[k].mu_2) > EPS) {
			fprintf(stderr, "ASYMMETRY XI_HALF1!\n");
		}
	}

	calculate_velocity_z_a_half_1(grid_prev, pi.ti_half_1, n, m);
	/*
	fprintf(stderr, "VELOCITY_Z_A_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_Z_B_HALF_1:\n");
	*/
	calculate_velocity_z_b_half_1(grid_prev, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_Z_C_HALF_1:\n");
	*/
	calculate_velocity_z_c_half_1(grid_prev, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "VELOCITY_Z_D_HALF_1:\n");
	*/
	calculate_velocity_z_d_half_1(grid_prev, grid_last_layer, pi.ti_half_1, n, m);
	/*
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].d[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			if (fabs(pi.ti_half_1[k].d[i] - pi.ti_half_1[k].d[pi.ti_half_1[k].n - i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZDHALF1 %zu %zu!\n", k, i);
			}
		}
	}

	prs_half_1(&pi);
	/*
	fprintf(stderr, "VELOCITY_Z_Y_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].y[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			if (fabs(pi.ti_half_1[k].y[i] - pi.ti_half_1[k].y[pi.ti_half_1[k].n - i]) > EPS) {
				fprintf(stderr, "ASYMMETRY Z_Y_HALF1!\n");
			}
		}
	}

	pi.n = n - 1;
	pi.ti_half_2 = calloc(pi.n, sizeof(tma_info));

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.n; ++i) {
		pi.ti_half_2[i].n = pi.m;
		pi.ti_half_2[i].a = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].b = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].c = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].d = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].y = calloc(m, sizeof(double));
		pi.ti_half_2[i].xi_1 = calculate_velocity_z_xi_1_half_2(grid_prev, i, m);
		pi.ti_half_2[i].xi_2 = 0; 
		pi.ti_half_2[i].mu_1 = calculate_velocity_z_mu_1_half_2(grid_prev, grid_last_layer, pi.ti_half_1, i, n, m);
		pi.ti_half_2[i].mu_2 = calculate_velocity_z_mu_2_half_2(grid_prev, i, m);
		//fprintf(stderr, "XI_1 = %lf, XI_2 = %lf, MU_1 = %lf, MU_2 = %lf\n", pi.ti_half_2[i].xi_1, pi.ti_half_2[i].xi_2, pi.ti_half_2[i].mu_1, pi.ti_half_2[i].mu_2);
	}

	for (size_t i = 0; i < pi.n; ++i) {
		if (fabs(pi.ti_half_2[i].xi_1 - pi.ti_half_2[pi.n - 1 - i].xi_1) > EPS || fabs(pi.ti_half_2[i].mu_1 - pi.ti_half_2[pi.n - 1 - i].mu_1) > EPS
			|| fabs(pi.ti_half_2[i].xi_2 - pi.ti_half_2[pi.n - 1 - i].xi_2) > EPS || fabs(pi.ti_half_2[i].mu_2 - pi.ti_half_2[pi.n - 1 - i].mu_2) > EPS) {
			fprintf(stderr, "ASYMMETRY XI_HALF2!\n");
		}
	}

	calculate_velocity_z_a_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_Z_A_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			if (fabs(pi.ti_half_2[k].a[i] - pi.ti_half_2[pi.n - 1 - k].a[i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZAHALF2!\n");
			}
		}
	}
	calculate_velocity_z_b_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_Z_B_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			if (fabs(pi.ti_half_2[k].b[i] - pi.ti_half_2[pi.n - 1 - k].b[i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZBHALF2!\n");
			}
		}
	}
	calculate_velocity_z_c_half_2(grid_prev, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_Z_C_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			if (fabs(pi.ti_half_2[k].c[i] - pi.ti_half_2[pi.n - 1 - k].c[i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZCHALF2!\n");
			}
		}
	}
	calculate_velocity_z_d_half_2(grid_prev, grid_last_layer, pi.ti_half_1, pi.ti_half_2, n, m);
	/*
	fprintf(stderr, "VELOCITY_Z_D_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].d[i]);
		}
		fprintf(stderr, "\n");
	}
	*/

	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			if (fabs(pi.ti_half_2[k].d[i] - pi.ti_half_2[pi.n - 1 - k].d[i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZDHALF2!\n");
			}
		}
	}
	prs_half_2(&pi);

	/*
	fprintf(stderr, "VELOCITY_Z_Y_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].y[i]);
		}
		fprintf(stderr, "\n");
	}
	*/
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			if (fabs(pi.ti_half_2[k].y[i] - pi.ti_half_2[pi.n - 1 - k].y[i]) > EPS) {
				fprintf(stderr, "ASYMMETRY ZYHALF2!\n");
			}
		}
	}

	double *tmp = calloc(n * m, sizeof(double));
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.n; ++i) {
		memcpy(&tmp[i * m], pi.ti_half_2[i].y, m * sizeof(double));
	}

	grid_fill(grid_cur, tmp, GRID_OFFSET_W);

	free(tmp);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.n; ++i) {
		free(pi.ti_half_2[i].a);
		free(pi.ti_half_2[i].b);
		free(pi.ti_half_2[i].c);
		free(pi.ti_half_2[i].d);
		free(pi.ti_half_2[i].y);
	}
	free(pi.ti_half_2);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.m; ++i) {
		free(pi.ti_half_1[i].a);
		free(pi.ti_half_1[i].b);
		free(pi.ti_half_1[i].c);
		free(pi.ti_half_1[i].d);
		free(pi.ti_half_1[i].y);
	}
	free(pi.ti_half_1);
}

static inline void calculate_p_a_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			double dx = grid[i * m + k].dx;
			double dx_half_minus = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;

			ti[k].a[i] = -1 / (dx * dx_half_minus);
		}
	}
}

static inline void calculate_p_b_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m, double tau)
{
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			double dx = grid[i * m + k].dx;
			double dx_half_minus = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
			double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
			
			ti[k].b[i] = 1 / (dx * dx_half_minus) + 1 / (dx * dx_half_plus) + 1 / tau;
		}
	}
}

static inline void calculate_p_c_half_1(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	for (size_t k = 0; k < m - 1; ++k) {
		for (size_t i = 1; i < n - 2; ++i) {
			double dx = grid[i * m + k].dx;
			double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;

			ti[k].c[i] = -1 / (dx * dx_half_plus);
		}
	}
}

static inline void calculate_p_d_half_1(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, tma_info *ti, size_t n, size_t m, double tau, double t)
{
	double density = get_blood_density();
	for (size_t k = 0; k < m - 1; ++k) {
		double x = 0;
		for (size_t i = 1; i < n - 2; ++i) {
			double p_prev = grid_last_layer[i * m + k].p;
			double p = grid_prev[i * m + k].p;
			double p_n = (k == m - 2 ? 0 : grid_prev[i * m + k + 1].p);
			double p_s = (k == 0 ? 0 : grid_prev[i * m + k - 1].p);
			double dz = grid[i * m + k].dz;
			double dz_half_plus = (k == m - 2 ? 0 : (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2);
			double dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
			double dx = grid[i * m + k].dx;
			double u = grid[i * m + k].u;
			double u_n = (k == m - 2 ? 0 : grid[i * m + k + 1].u);
			double u_s = (k == 0 ? 0 : grid[i * m + k - 1].u);
			double u_e = grid[(i + 1) * m + k].u;
			double u_ne = (k == m - 2 ? 0 : grid[(i + 1) * m + k + 1].u);
			double u_se = (k == 0 ? 0 : grid[(i + 1) * m + k - 1].u);
			double w = grid[i * m + k].w;
			double w_n = grid[i * m + k + 1].w;
			double w_w = grid[(i - 1) * m + k].w;
			double w_e = grid[(i + 1) * m + k].w;
			double w_ne = grid[(i + 1) * m + k + 1].w;
			double w_nw = grid[(i - 1) * m + k + 1].w;
			double mu = grid_prev[i * m + k].mu;
		
			if (k == 0) {
				ti[k].d[i] = p_prev / tau + (p_n - p) / (dz * dz_half_plus) -
					2 * density * ((u_e - u) * (w_n - w) / (dz * dx) - 
					(u_n - u + u_ne - u_e) * (w_e - w_w + w_ne - w_nw) / (16 * dz * dx));
			} else if (k == m - 2) {
				double R = get_vessel_size_x() / 2;
				double x_ = x + dx / 2;
				
				ti[k].d[i] = p_prev / tau + 1 / dz * ((density * w_0_derivative_t(t) * x_ * (2 - x_ / R) / R +
							mu * 2 * w_0(t) / (R * R)) - (p - p_s) / dz_half_minus) -
					2 * density * ((u_e - u) * (w_n - w) / (dz * dx) - 
					(-u_s - u_se) * (w_e - w_w + w_ne - w_nw) / (16 * dz * dx));
			} else {
				ti[k].d[i] = p_prev / tau + 1 / dz * ((p_n - p) / dz_half_plus - (p - p_s) / dz_half_minus) -
					2 * density * ((u_e - u) * (w_n - w) / (dz * dx) - 
					(u_n - u_s + u_ne - u_se) * (w_e - w_w + w_ne - w_nw) / (16 * dz * dx));
			}
			x += dx;
		}
	}
}

static inline double calculate_p_xi_1_half_1(grid_node *grid, size_t k, size_t m, double tau)
{
	size_t i = 0;
	double dx = grid[i * m + k].dx;
	double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
	
	return tau / (tau + dx * dx_half_plus);
}

static inline double calculate_p_mu_1_half_1(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, size_t k, size_t m, double t, double tau)
{
	ssize_t i = 0;
	double density = get_blood_density();
	double mu = grid_prev[i * m + k].mu;
	double dz = grid[i * m + k].dz;
	double dz_half_plus = (k == m - 2 ? 0 : (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2);
	double dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
	double dx = grid[i * m + k].dx;
	double dx_e = grid[(i + 1) * m + k].dx;
	double dx_half_plus = (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2;
	double u_e = grid[(i + 1) * m + k].u;
	double u_ee = grid[(i + 2) * m + k].u;
	double u_se = (k == 0 ? 0 : grid[(i + 1) * m + k - 1].u);
	double u_ne = (k == m - 2 ? 0 : grid[(i + 1) * m + k + 1].u);
	double w = grid[i * m + k].w;
	double w_e = grid[(i + 1) * m + k].w;
	double w_n = grid[i * m + k + 1].w;
	double w_ne = grid[(i + 1) * m + k + 1].w;
	double p = grid_prev[i * m + k].p;
	double p_prev = grid_last_layer[i * m + k].p;
	double p_n = (k == m - 2 ? 0 : grid_prev[i * m + k + 1].p);
	double p_s = (k == 0 ? 0 : grid_prev[i * m + k - 1].p);

	if (k == 0) {
		return (p_prev / tau - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) / (dx * dx * dx_e * (dx + dx_e)) +
			(p_n - p) / (dz * dz_half_plus) -
			2 * density * (u_e * (w_n - w) -
			(u_ne - u_e) * (w_e + w_ne + w + w_n) / 16) / (dz * dx)) /
			(1 / tau + 1 / (dx * dx_half_plus));
	} else if (k == m - 2) {
		double x = grid[i * m + k].dx / 2;
		double R = get_vessel_size_x() / 2;
		return (p_prev / tau - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) / (dx * dx * dx_e * (dx + dx_e)) +
			(density * w_0_derivative_t(t) * x * (2 - x / R) / R + mu * 2 * w_0(t) / (R * R) - (p - p_s) / dz_half_minus) / dz -
			2 * density * (u_e * (w_n - w) +
					(u_se * (w_e + w_ne + w_n + w) / 16)) / (dx * dz)) /
			(1 / tau + 1 / (dx * dx_half_plus));
	} else {
		return (p_prev / tau - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) / (dx * dx * dx_e * (dx + dx_e)) +
			((p_n - p) / dz_half_plus - (p - p_s) / dz_half_minus) / dz -
			2 * density * (u_e * (w_n - w) -
			(u_ne - u_se) * (w_e + w_ne + w + w_n) / 16) / (dz * dx)) /
			(1 / tau + 1 / (dx * dx_half_plus));
	}
}

static inline double calculate_p_xi_2_half_1(grid_node *grid, size_t k, ssize_t n, size_t m, double tau)
{
	ssize_t i = n - 2;
	double dx = grid[i * m + k].dx;
	double dx_half_minus = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
	
	return tau / (tau + dx * dx_half_minus);
}

static inline double calculate_p_mu_2_half_1(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, size_t k, ssize_t n, size_t m, double t, double tau)
{
	ssize_t i = n - 2;
	double density = get_blood_density();
	double mu = grid_prev[i * m + k].mu;
	double dz = grid[i * m + k].dz;
	double dz_half_plus = (k == m - 2 ? 0 : (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2);
	double dz_half_minus = (k == 0 ? 0 : (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2);
	double dx = grid[i * m + k].dx;
	double dx_w = grid[(i - 1) * m + k].dx;
	double dx_half_minus = (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2;
	double u_w = grid[(i - 1) * m + k].u;
	double u_n = (k == m - 2 ? 0 : grid[i * m + k + 1].u);
	double u_s = (k == 0 ? 0 : grid[i * m + k - 1].u);
	double u = grid[i * m + k].u;
	double w = grid[i * m + k].w;
	double w_w = grid[(i - 1) * m + k].w;
	double w_n = grid[i * m + k + 1].w;
	double w_nw = grid[(i - 1) * m + k + 1].w;
	double p = grid_prev[i * m + k].p;
	double p_prev = grid_last_layer[i * m + k].p;
	double p_n = (k == m - 2 ? 0 : grid_prev[i * m + k + 1].p);
	double p_s = (k == 0 ? 0 : grid_prev[i * m + k - 1].p);

	if (k == 0) {
		double u = grid[i * m + k].u;
		return (p_prev / tau + 4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx * dx_w * (dx + dx_w)) +
			(p_n - p) / (dz * dz_half_plus) -
			2 * density * (-u * (w_n - w) -
			(u_n - u) * (-w_w - w_nw - w - w_n) / 16) / (dz * dx)) /
			(1 / tau + 1 / (dx * dx_half_minus));
	} else if (k == m - 2) {
		double R = get_vessel_size_x() / 2;
		double x = 2 * R - dx / 2; 
		return (p_prev / tau + 4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx * dx_w * (dx + dx_w)) +
			(density * w_0_derivative_t(t) * 
			x * (2 - x / R) / R + mu * 2 * w_0(t) / (R * R) - (p - p_s) / dz_half_minus) / dz -
			2 * density * (-u * (w_n - w) +
					(u_s * (-w_w - w_nw - w - w_n) / 16)) / (dx * dz)) /
			(1 / tau + 1 / (dx * dx_half_minus));
	} else {
		return (p_prev / tau + 4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx * dx_w * (dx + dx_w)) +
			((p_n - p) / dz_half_plus - (p - p_s) / dz_half_minus) / dz -
			2 * density * (-u * (w_n - w) -
			(u_n - u_s) * (-w_w - w_nw - w - w_n) / 16) / (dz * dx)) / 
			(1 / tau + 1 / (dx * dx_half_minus));
	}
}

static inline void calculate_p_a_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			double dz = grid[i * m + k].dz;
			double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;

			ti[i].a[k] = 1 / (dz * dz_half_minus);
		}
	}
}

static inline void calculate_p_b_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m, double tau)
{
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			double dz = grid[i * m + k].dz;
			double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
			double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;

			ti[i].b[k] = -1 / (dz * dz_half_minus) - 1 / (dz * dz_half_plus) + 1 / tau;
		}
	}
}

static inline void calculate_p_c_half_2(grid_node *grid, tma_info *ti, size_t n, size_t m)
{
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			double dz = grid[i * m + k].dz;
			double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;

			ti[i].c[k] = 1 / (dz * dz_half_plus);
		}
	}
}

static inline void calculate_p_d_half_2(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, tma_info *ti_half_1, tma_info *ti_half_2, size_t n, size_t m, double tau)
{
	double density = get_blood_density();
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t k = 1; k < m - 2; ++k) {
			double p_prev = grid_last_layer[i * m + k].p;
			double p = ti_half_1[k].y[i];
			double p_w = (i == 0 ? 0 : ti_half_1[k].y[i - 1]);
			double p_e = (i == n - 2 ? 0 : ti_half_1[k].y[i + 1]);
			double mu = grid_prev[i * m + k].mu;
			double dz = grid[i * m + k].dz;
			double dx = grid[i * m + k].dx;
			double dx_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].dx);
			double dx_w = (i == 0 ? 0 : grid[(i - 1) * m + k].dx);
			double dx_half_plus = (i == n - 2 ? 0 : (grid[(i + 1) * m + k].dx + grid[i * m + k].dx) / 2);
			double dx_half_minus = (i == 0 ? 0 : (grid[(i - 1) * m + k].dx + grid[i * m + k].dx) / 2);
			double u = grid[i * m + k].u;
			double u_n = grid[i * m + k + 1].u;
			double u_s = grid[i * m + k - 1].u;
			double u_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].u);
			double u_w = (i == 0 ? 0 : grid[(i - 1) * m + k].u);
			double u_ee = (i == n - 2 ? 0 : grid[(i + 2) * m + k].u);
			double u_ne = grid[(i + 1) * m + k + 1].u;
			double u_se = grid[(i + 1) * m + k - 1].u;
			double w = grid[i * m + k].w;
			double w_n = grid[i * m + k + 1].w;
			double w_w = (i == 0 ? 0 : grid[(i - 1) * m + k].w);
			double w_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].w);
			double w_ne = (i == n - 2 ? 0 : grid[(i + 1) * m + k + 1].w);
			double w_nw = (i == 0 ? 0 : grid[(i - 1) * m + k + 1].w);

			if (i == 0) {
				ti_half_2[i].d[k] = p_prev / tau + ((p_e - p) / dx_half_plus - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) / (dx * dx_e * (dx + dx_e))) / dx -
					2 * density * (u_e * (w_n - w) / (dx * dz) - 
					(u_ne - u_se) * (w_e + w + w_ne + w_n) / (16 * dz * dx));
			} else if (i == n - 2) {
				ti_half_2[i].d[k] = p_prev / tau + (4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx_w * (dx + dx_w)) - (p - p_w) / dx_half_minus) / dx -
					2 * density * ((-u) * (w_n - w) / (dx * dz) - 
					(u_n - u_s) * (-w - w_w - w_n - w_nw) / (16 * dz * dx));
			} else {
				ti_half_2[i].d[k] = p_prev / tau + ((p_e - p) / dx_half_plus - (p - p_w) / dx_half_minus) / dx -
					2 * density * ((u_e - u) * (w_n - w) / (dx * dz) - 
					(u_n - u_s + u_ne - u_se) * (w_e - w_w + w_ne - w_nw) / (16 * dz * dx));
			}
		}
	}
}

static inline double calculate_p_xi_1_half_2(grid_node *grid, size_t i, size_t m, double tau)
{
	ssize_t k = 0;
	double dz = grid[i * m + k].dz;
	double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
	
	return tau / (tau + dz * dz_half_plus);
}

static inline double calculate_p_mu_1_half_2(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, tma_info *ti_half_1, size_t i, size_t n, size_t m, double tau)
{
	size_t k = 0;
	double mu = grid_prev[i * m + k].mu;
	double density = get_blood_density();
	double dz = grid[i * m + k].dz;
	double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
	double dx = grid[i * m + k].dx;
	double dx_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].dx);
	double dx_w = (i == 0 ? 0 : grid[(i - 1) * m + k].dx);
	double dx_half_plus = (i == n - 2 ? 0 : (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2);
	double dx_half_minus = (i == 0 ? 0 : (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2);
	double p_prev = grid_last_layer[i * m + k].p;
	double p = ti_half_1[k].y[i];
	double p_w = (i == 0 ? 0 : ti_half_1[k].y[i - 1]);
	double p_e = (i == n - 2 ? 0 : ti_half_1[k].y[i + 1]);
	double u = grid[i * m + k].u;
	double u_n = grid[i * m + k + 1].u;
	double u_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].u);
	double u_w = (i == 0 ? 0 : grid[(i - 1) * m + k].u);
	double u_ee = (i == n - 2 ? 0 : grid[(i + 2) * m + k].u);
	double u_ne = (i == n - 2 ? 0 : grid[(i + 1) * m + k + 1].u);
	double w = grid[i * m + k].w;
	double w_n = grid[i * m + k + 1].w;
	double w_w = (i == 0 ? 0 : grid[(i - 1) * m + k].w);
	double w_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].w);
	double w_ne = (i == n - 2 ? 0 : grid[(i + 1) * m + k + 1].w);
	double w_nw = (i == 0 ? 0 : grid[(i - 1) * m + k + 1].w);

	if (i == 0) {
		return (p_prev / tau + ((p_e - p) / dx_half_plus - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) /
				(dx * dx_e * (dx + dx_e))) / dx -
			2 * density * (u_e * (w_n - w) - (u_ne - u_e) *
			(w_e + w_ne + w_n + w) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_plus));
	} else if (i == n - 2) {
		return (p_prev / tau + (4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx_w * (dx + dx_w))
				- (p - p_w) / dx_half_minus) / dx -
			2 * density * (-u * (w_n - w) - (u_n - u) *
			(-w - w_n - w_w - w_nw) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_plus));
	} else {
		return (p_prev / tau + ((p_e - p) / dx_half_plus - (p - p_w) / dx_half_minus) / dx -
			2 * density * ((u_e - u) * (w_n - w) - (u_n + u_ne - u - u_e) *
			(w_e + w_ne - w_w - w_nw) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_plus));
	}
}

static inline double calculate_p_xi_2_half_2(grid_node *grid, size_t i, size_t m, double tau)
{
	ssize_t k = m - 2;
	double dz = grid[i * m + k].dz;
	double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
	
	return tau / (tau + dz * dz_half_minus);
}

static inline double calculate_p_mu_2_half_2(grid_node *grid, grid_node *grid_prev, grid_node *grid_last_layer, tma_info *ti_half_1, size_t i, size_t n, size_t m, double t, double tau)
{
	ssize_t k = m - 2;
	double mu = grid_prev[i * m + k].mu;
	double density = get_blood_density();
	double dz = grid[i * m + k].dz;
	double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
	double dx = grid[i * m + k].dx;
	double dx_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].dx);
	double dx_w = (i == 0 ? 0 : grid[(i - 1) * m + k].dx);
	double dx_half_plus = (i == n - 2 ? 0 : (grid[i * m + k].dx + grid[(i + 1) * m + k].dx) / 2);
	double dx_half_minus = (i == 0 ? 0 : (grid[i * m + k].dx + grid[(i - 1) * m + k].dx) / 2);
	double p_prev = grid_last_layer[i * m + k].p;
	double p = ti_half_1[k].y[i];
	double p_w = (i == 0 ? 0 : ti_half_1[k].y[i - 1]);
	double p_e = (i == n - 2 ? 0 : ti_half_1[k].y[i + 1]);
	double u = grid[i * m + k].u;
	double u_s = grid[i * m + k - 1].u;
	double u_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].u);
	double u_w = (i == 0 ? 0 : grid[(i - 1) * m + k].u);
	double u_ee = (i == n - 2 ? 0 : grid[(i + 2) * m + k].u);
	double u_se = (i == n - 2 ? 0 : grid[(i + 1) * m + k - 1].u);
	double w = grid[i * m + k].w;
	double w_n = grid[i * m + k + 1].w;
	double w_w = (i == 0 ? 0 : grid[(i - 1) * m + k].w);
	double w_e = (i == n - 2 ? 0 : grid[(i + 1) * m + k].w);
	double w_ne = (i == n - 2 ? 0 : grid[(i + 1) * m + k + 1].w);
	double w_nw = (i == 0 ? 0 : grid[(i - 1) * m + k + 1].w);
	double R = get_vessel_size_x() / 2;

	if (i == 0) {
		double x = dx / 2;
		return (p_prev / tau + ((p_e - p) / dx_half_plus - 4 * mu * (dx * u_ee - (dx + dx_e) * u_e) / (dx * dx_e * (dx + dx_e))) / dx +
			(density * w_0_derivative_t(t) * x * (2 - x / R) / R + mu * 2 * w_0(t) / (R * R)) / dz -
			2 * density * (u_e * (w_n - w) - (-u_se) *
			(w_e + w_ne + w + w_n) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_minus));
	} else if (i == n - 2) {
		double x = 2 * R - dx / 2;
		return (p_prev / tau + (4 * mu * (dx * u_w - (dx + dx_w) * u) / (dx * dx_w * (dx + dx_w)) - (p - p_w) / dx_half_minus) / dx +
			(density * w_0_derivative_t(t) * x * (2 - x / R) / R + mu * 2 * w_0(t) / (R * R)) / dz -
			2 * density * ((-u) * (w_n - w) - (-u_s) *
			(-w - w_n - w_w - w_nw) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_minus));
	} else {
		double x = 0;
		for (size_t l = 0; l < i; ++l) {
			x += grid[l * m + k].dx;
		}
		x += dx / 2;
		return (p_prev / tau + ((p_e - p) / dx_half_plus - (p - p_w) / dx_half_minus) / dx +
			(density * w_0_derivative_t(t) * x * (2 - x / R) / R + mu * 2 * w_0(t) / (R * R)) / dz -
			2 * density * ((u_e - u) * (w_n - w) - (-u_s - u_se) *
			(w_e + w_ne - w_w - w_nw) / 16) / (dx * dz)) /
			(1 / tau + 1 / (dz * dz_half_minus));
	}
}

static void calculate_p(grid_node *grid_prev, grid_node *grid_cur, grid_node *grid_last_layer, size_t n, size_t m, double tau, double t)
{
	prs_info pi;
	pi.n = n - 2;
	pi.m = m - 1;
	pi.start_half_1 = 0;
	pi.start_half_2 = 0;

	pi.ti_half_1 = calloc(pi.m, sizeof(tma_info));

	for (size_t k = 0; k < pi.m; ++k) {
		pi.ti_half_1[k].n = pi.n;
		pi.ti_half_1[k].a = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].b = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].c = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].d = calloc((pi.ti_half_1[k].n + 1), sizeof(double));
		pi.ti_half_1[k].y = calloc(n, sizeof(double));
		pi.ti_half_1[k].xi_1 = calculate_p_xi_1_half_1(grid_cur, k, m, tau);
		pi.ti_half_1[k].xi_2 = calculate_p_xi_2_half_1(grid_cur, k, n, m, tau);
		pi.ti_half_1[k].mu_1 = calculate_p_mu_1_half_1(grid_cur, grid_prev, grid_last_layer, k, m, t, tau);
		pi.ti_half_1[k].mu_2 = calculate_p_mu_2_half_1(grid_cur, grid_prev, grid_last_layer, k, n, m, t, tau);
		fprintf(stderr, "XI_1 = %lf, XI_2 = %lf, MU_1 = %lf, MU_2 = %lf\n", pi.ti_half_1[k].xi_1, pi.ti_half_1[k].xi_2, pi.ti_half_1[k].mu_1, pi.ti_half_1[k].mu_2);
	}

	calculate_p_a_half_1(grid_cur, pi.ti_half_1, n, m);
	fprintf(stderr, "PRESSURE_A_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "PRESSURE_B_HALF_1:\n");
	calculate_p_b_half_1(grid_cur, pi.ti_half_1, n, m, tau);
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "PRESSURE_C_HALF_1:\n");
	calculate_p_c_half_1(grid_cur, pi.ti_half_1, n, m);
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "PRESSURE_D_HALF_1:\n");
	calculate_p_d_half_1(grid_cur, grid_prev, grid_last_layer, pi.ti_half_1, n, m, tau, t);
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].d[i]);
		}
		fprintf(stderr, "\n");
	}

	prs_half_1(&pi);
	fprintf(stderr, "PRESSURE_Y_HALF_1:\n");
	for (size_t k = 0; k < pi.m; ++k) {
		for (size_t i = 0; i <= pi.ti_half_1[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_1[k].y[i]);
		}
		fprintf(stderr, "\n");
	}

	pi.n = n - 1;

	pi.ti_half_2 = calloc(pi.n, sizeof(tma_info));
	for (size_t i = 0; i < pi.n; ++i) {
		pi.ti_half_2[i].n = pi.m - 1;
		pi.ti_half_2[i].a = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].b = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].c = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].d = calloc((pi.ti_half_2[i].n + 1), sizeof(double));
		pi.ti_half_2[i].y = calloc(m, sizeof(double));
		pi.ti_half_2[i].xi_1 = calculate_p_xi_1_half_2(grid_cur, i, m, tau);
		pi.ti_half_2[i].xi_2 = calculate_p_xi_2_half_2(grid_cur, i, m, tau);
		pi.ti_half_2[i].mu_1 = calculate_p_mu_1_half_2(grid_cur, grid_prev, grid_last_layer, pi.ti_half_1, i, n, m, tau);
		pi.ti_half_2[i].mu_2 = calculate_p_mu_2_half_2(grid_cur, grid_prev, grid_last_layer, pi.ti_half_1, i, n, m, t, tau);
		fprintf(stderr, "XI_1 = %lf, XI_2 = %lf, MU_1 = %lf, MU_2 = %lf\n", pi.ti_half_2[i].xi_1, pi.ti_half_2[i].xi_2, pi.ti_half_2[i].mu_1, pi.ti_half_2[i].mu_2);
	}

	calculate_p_a_half_2(grid_cur, pi.ti_half_2, n, m);
	fprintf(stderr, "PRESSURE_A_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].a[i]);
		}
		fprintf(stderr, "\n");
	}
	calculate_p_b_half_2(grid_cur, pi.ti_half_2, n, m, tau);
	fprintf(stderr, "PRESSURE_B_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].b[i]);
		}
		fprintf(stderr, "\n");
	}
	calculate_p_c_half_2(grid_cur, pi.ti_half_2, n, m);
	fprintf(stderr, "PRESSURE_C_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].c[i]);
		}
		fprintf(stderr, "\n");
	}
	calculate_p_d_half_2(grid_cur, grid_prev, grid_last_layer, pi.ti_half_1, pi.ti_half_2, n, m, tau);
	fprintf(stderr, "PRESSURE_D_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].d[i]);
		}
		fprintf(stderr, "\n");
	}

	prs_half_2(&pi);

	fprintf(stderr, "PRESSURE_Y_HALF_2:\n");
	for (size_t k = 0; k < pi.n; ++k) {
		for (size_t i = 0; i <= pi.ti_half_2[k].n; ++i) {
			fprintf(stderr, "%lf ", pi.ti_half_2[k].y[i]);
		}
		fprintf(stderr, "\n");
	}

	memcpy(pi.ti_half_2[0].y, pi.ti_half_2[1].y, m * sizeof(double));
	memcpy(pi.ti_half_2[pi.n - 1].y, pi.ti_half_2[pi.n - 2].y, m * sizeof(double));


	double *tmp = calloc(n * m, sizeof(double));
	for (size_t i = 0; i < pi.n; ++i) {
		memcpy(&tmp[i * m], pi.ti_half_2[i].y, m * sizeof(double));
	}

	grid_fill(grid_cur, tmp, GRID_OFFSET_P);

	for (size_t i = 0; i < n; ++i) {
		grid_cur[i * m].p = grid_cur[i * m + 1].p;
		grid_cur[i * m + m - 2].p = grid_cur[i * m + m - 3].p;
	}

	free(tmp);
	for (size_t i = 0; i < pi.n; ++i) {
		free(pi.ti_half_2[i].a);
		free(pi.ti_half_2[i].b);
		free(pi.ti_half_2[i].c);
		free(pi.ti_half_2[i].d);
		free(pi.ti_half_2[i].y);
	}
	free(pi.ti_half_2);

	for (size_t i = 0; i < pi.m; ++i) {
		free(pi.ti_half_1[i].a);
		free(pi.ti_half_1[i].b);
		free(pi.ti_half_1[i].c);
		free(pi.ti_half_1[i].d);
		free(pi.ti_half_1[i].y);
	}
	free(pi.ti_half_1);
}

static inline double init_w(double t, double x)
{
	double R = get_vessel_size_x() / 2.0;

	return w_0(t) * (x * (2.0 - x / R)) / R; 
}

static void calculate_boundary_conditions(grid_node *grid, size_t n, size_t m, double t)
{
#ifdef WITH_OPENMP_
#pragma omp parallel 
{
#pragma omp for nowait
#endif
	for (size_t j = 0; j < m; ++j) {
		grid[j].u = 0;
		//grid[j].w = 0;
		grid[(n - 1) * m + j].u = 0;
		//grid[(n - 1) * m + j].w = 0;
	}

	//fprintf(stderr, "W_INIT:\n");
	double cur_x = 0.0;
	for (size_t i = 0; i < n; ++i) {
		cur_x += grid[i * m + m - 1].dx / 2;
		grid[i * m + m - 1].u = 0;
		grid[i * m + m - 1].w = -init_w(t, cur_x);
	//	fprintf(stderr, "%.18lf ", grid[i * m + m - 1].w);
		cur_x += grid[i * m + m - 1].dx / 2;
	}
	//fprintf(stderr, "\n");

	// DEBUG
#ifdef WITH_OPENMP_
#pragma omp for collapse(2) nowait
#endif
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			grid[i * m + j].p = 13332;
#ifdef WITH_OPENMP_
}
#endif
}

static bool check_convergence(grid_node *prev, grid_node *cur, size_t n, size_t m, double eps)
{
	int res = 0;
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2) reduction(+:res)
#endif
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			if (fabs(prev[i * m + j].u - cur[i * m + j].u) > eps ||
				fabs(prev[i * m + j].w - cur[i * m + j].w) > eps ||
				fabs(prev[i * m + j].p - cur[i * m + j].p) > eps) {
	//			fprintf(stderr, "DIVERGE: U = %lf, W = %lf, P = %lf\n", fabs(prev[i * m + j].u - cur[i * m + j].u),
	//					fabs(prev[i * m + j].w - cur[i * m + j].w), fabs(prev[i * m + j].p - cur[i * m + j].p));
				res += 1;
			}
		}
	}
	//fprintf(stderr, "RES %d\n", res);
	return (res == 0);
}

void calculate(double sigma, double t_beg, double t_end)
{
	double delta_t = get_delta_t();
	grid_node *grid_cur_layer = NULL;
	grid_node *grid_next_layer = NULL;

	size_t layers;
	layers = ceil((t_end - t_beg) / delta_t);
	
	size_t n, m;
	grid_cur_layer = grid_generate(sigma, &n, &m);
	grid_fill_from_config(grid_cur_layer);

	grid_node *grid_prev = NULL;
	grid_node *grid_cur = NULL;

	grid_cur = grid_generate(sigma, &n, &m);

	for (size_t cur_layer = 0; cur_layer + 1 < layers; ++cur_layer) {
		grid_copy(grid_cur_layer, &grid_prev);

		calculate_boundary_conditions(grid_cur, n, m, cur_layer * delta_t);
		calculate_boundary_conditions(grid_prev, n, m, cur_layer * delta_t);
		calculate_boundary_conditions(grid_cur_layer, n, m, cur_layer * delta_t);

		fprintf(stderr, "LAYER #%zu\n", cur_layer);
		grid_print_all(grid_cur_layer, n, m, "text", GRID_PRINT_AS_TABLE);

		size_t iter = 0;
		while (true) {
			fprintf(stderr, "ITER #%zu\n", iter);
			calculate_mu(grid_prev, n, m);
			calculate_velocity_z(grid_prev, grid_cur, grid_cur_layer, n, m);
			calculate_boundary_conditions(grid_cur, n, m, cur_layer * delta_t);
			calculate_boundary_conditions(grid_prev, n, m, cur_layer * delta_t);
			calculate_velocity_x(grid_prev, grid_cur, grid_cur_layer, n, m);
			calculate_boundary_conditions(grid_cur, n, m, cur_layer * delta_t);
			calculate_boundary_conditions(grid_prev, n, m, cur_layer * delta_t);
//			calculate_p(grid_prev, grid_cur, grid_cur_layer, n, m, 0.1, (cur_layer + 1) * delta_t);
			//fprintf(stderr, "GRID_PREV\n");
			//grid_print_all(grid_prev, n, m, GRID_PRINT_AS_TABLE);
			//fprintf(stderr, "GRID_CUR\n");
			//grid_print_all(grid_cur, n, m, GRID_PRINT_AS_TABLE);
			//grid_print_all(grid_cur, n, m, GRID_PRINT_AS_TABLE);
			//exit(0);
			if (check_convergence(grid_prev, grid_cur, n, m, 1e-9)) {
				grid_copy(grid_cur, &grid_next_layer);
				break;
			} else {
				grid_node *tmp = grid_prev;
				grid_prev = grid_cur;
				grid_cur = tmp;
			}
			++iter;
		}
		grid_node *tmp = grid_cur_layer;
		grid_cur_layer = grid_next_layer;
		grid_next_layer = tmp;
	}
	fprintf(stderr, "LAYER #last\n");
	grid_print_all(grid_cur_layer, n, m, "text", GRID_PRINT_AS_TABLE);
	grid_destroy(grid_prev);
	grid_destroy(grid_cur);
}

int main(void)
{
	log_open(NULL, LOG_STDERR);

	read_config("config.ex1");
//	setvbuf(stdout, NULL, _IONBF, 0);

	double t1 = omp_get_wtime();
	calculate(5, 0, 1);
	double t2 = omp_get_wtime();

	fprintf(stderr, "TOTAL=%lf\n", t2 - t1);

	clear_config();
	grid_clear();
}
