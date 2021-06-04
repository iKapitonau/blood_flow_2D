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

	/*
	if (t < 0.2598 || t >= 0.7523)
		f = 0.189;
	else if (t >= 0.2598 && t < 0.458) 
		f = -59.5 * pow(t - 0.37, 2) + 0.912;
	else if (t >= 0.458 && t < 0.5901) 
		f = -27.7 * pow(t - 0.5, 2) + 0.5;
	else if (t >= 0.5901 && t < 0.7523)
		f = -23.7 * pow(t - 0.66, 2) + 0.391;
	*/
	
	if ((t > 0.213793 && t < 0.5) || t > 0.713793)
		f = 0.1;
	else if (t <= 0.213793)
		f = -35 * pow(t - 0.1069, 2) + 0.5;
	else
		f = -35 * pow(t - 0.5 - 0.1069, 2) + 0.5;

	return f;
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
		return 0;
	else if (t <= 0.213793)
		return -35 * 2 * (t - 0.1069);
	else
		return -35 * 2 * (t - 0.5 - 0.1069);
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
			
		//	CARREAU-YASUDA
			if (fabs(shear_rate) <= EPS)
				grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF);
			else
				grid[i * m + k].mu = MU_INF + (MU_0 - MU_INF) *
					pow(1 + pow(B * fabs(shear_rate), A), (N - 1) / A);
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

			ti[k].a[i] = -delta_t * (2 * mu_w / (density * dx_w) + c1 * d1) / dx_half;
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

			ti[k].b[i] = 1 + delta_t * (2 * mu / (density * dx) +
				2 * mu_w / (density * dx_w) + ((1 - c2) * d2 - (1 - c1) * d1)) / dx_half;
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
			
			ti[k].c[i] = delta_t * (-2 * mu / (density * dx) + c2 * d2) / dx_half;
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
					(c4 * d4 * u_n + ((1 - c4) * d4 - d3) * u) / dz - (p - p_w) / (dx_half * density) +
					(mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)) / (density * dz * dx_half));
			} else if (k == m - 2) {
				ti[k].d[i] = u_prev + delta_t * ((-mu_avg_n * u / dz - mu_avg * (u - u_s) / dz_half_minus) / (dz * density) +
					((1 - c3) * d3 * u + c3 * d3 * u_s) / dz - (p - p_w) / (dx_half * density) +
					(mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)) / (density * dz * dx_half));
			} else {
				ti[k].d[i] = u_prev + delta_t * ((mu_avg_n * (u_n - u) / dz_half_plus - mu_avg * (u - u_s) / dz_half_minus) / (dz * density) -
					(c4 * d4 * u_n + ((1 - c4) * d4 - (1 - c3) * d3) * u - c3 * d3 * u_s) / dz - (p - p_w) / (dx_half * density) + (mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)) / (density * dz * dx_half));
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

			ti[i].a[k] = delta_t * (-mu_avg / (density * dz_half) - c3 * d3) / dz;
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

			ti[i].b[k] = 1 + delta_t * (mu_avg_n / (density * dz_half_plus) +
				mu_avg / (density * dz_half_minus) + ((1 - c4) * d4 - (1 - c3) * d3)) / dz;
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

			ti[i].c[k] = delta_t * (-mu_avg_n / (density * dz_half_plus) + c4 * d4) / dz;
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

			ti_half_2[i].d[k] = u_prev + delta_t * (2 * (mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) / (density * dx_half) -
				(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) / dx_half -
				(p - p_w) / (dx_half * density) + (mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)) / (density * dx_half * dz)); 
		}
	}
}

static inline double calculate_velocity_x_xi_1_half_2(grid_node *grid, size_t i, size_t m)
{
	size_t k = 0;
	double d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
	double d4 = (grid[(i - 1) * m + k + 1].w + grid[i * m + k + 1].w) / 2;
	double c4 = (d4 >= 0 ? 0 : 1);
	double mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
			grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4;
	double dz = grid[i * m + k].dz;
	double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return delta_t * (mu_avg_n / (density * dz_half_plus) - c4 * d4) /
			(dz + delta_t * (mu_avg_n / (density * dz_half_plus) + ((1 - c4) * d4 - d3)));
}

static inline double calculate_velocity_x_xi_2_half_2(grid_node *grid, size_t i, ssize_t m)
{
	ssize_t k = m - 2;
	double d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
	double c3 = (-d3 >= 0 ? 0 : 1);
	double dz = grid[i * m + k].dz;
	double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
	double mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2;
	double mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
			grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;
	double density = get_blood_density();
	double delta_t = get_delta_t();

	return delta_t * (c3 * d3 + mu_avg / (density * dz_half_minus)) / 
		(dz + delta_t * (mu_avg_n / (density * dz) + mu_avg / (density * dz_half_minus) - (1 - c3) * d3));
}

static inline double calculate_velocity_x_mu_1_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, size_t i, size_t m)
{
	size_t k = 0;
	double d1 = (grid[i * m + k].u + grid[(i - 1) * m + k].u) / 2;
	double d2 = (grid[i * m + k].u + grid[(i + 1) * m + k].u) / 2;
	double d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
	double d4 = (grid[(i - 1) * m + k + 1].w + grid[i * m + k + 1].w) / 2;
	double c1 = (-d1 >= 0 ? 0 : 1);
	double c2 = (d2 >= 0 ? 0 : 1);
	double c4 = (d4 >= 0 ? 0 : 1);
	double dx = grid[i * m + k].dx;
	double dx_w = grid[(i - 1) * m + k].dx;
	double dx_half = (dx + dx_w) / 2;
	double dz_half_plus = (grid[i * m + k].dz + grid[i * m + k + 1].dz) / 2;
	double dz = grid[i * m + k].dz;
	double u_prev = grid_prev[i * m + k].u;
	double u_e = ti_half_1[k].y[i + 1];
	double u_w = ti_half_1[k].y[i - 1];
	double u = ti_half_1[k].y[i];
	double mu = grid[i * m + k].mu;
	double mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2;
	double mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
			grid[i * m + k + 1].mu + grid[(i - 1) * m + k + 1].mu) / 4;
	double mu_w = grid[(i - 1) * m + k].mu;
	double p = grid[i * m + k].p;
	double p_w = grid[(i - 1) * m + k].p;
	double w_n = grid[i * m + k + 1].w;
	double w_nw = grid[(i - 1) * m + k + 1].w;
	double w = grid[i * m + k].w;
	double w_w = grid[(i - 1) * m + k].w;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return (u_prev + delta_t * (2 * (mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) / (density * dx_half) -
			(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) / dx_half -
			(p - p_w) / (density * dx_half) + (mu_avg_n * (w_n - w_nw) / dx_half - mu_avg * (w - w_w) / dx_half) / (density * dz))) /
			(1 + delta_t * (mu_avg_n / (density * dz_half_plus) + ((1 - c4) * d4 - d3)) / dz);
}

static inline double calculate_velocity_x_mu_2_half_2(grid_node *grid, grid_node *grid_prev, tma_info *ti_half_1, size_t i, ssize_t m)
{
	ssize_t k = m - 2;
	double d1 = (grid[i * m + k].u + grid[(i - 1) * m + k].u) / 2;
	double d2 = (grid[i * m + k].u + grid[(i + 1) * m + k].u) / 2;
	double d3 = (grid[i * m + k].w + grid[(i - 1) * m + k].w) / 2;
	double c1 = (-d1 >= 0 ? 0 : 1);
	double c2 = (d2 >= 0 ? 0 : 1);
	double c3 = (-d3 >= 0 ? 0 : 1);
	double dz = grid[i * m + k].dz;
	double dz_half_minus = (grid[i * m + k].dz + grid[i * m + k - 1].dz) / 2;
	double dx_half = (grid[i * m + k].dx + grid[(i - 1) * m + k - 1].dx) / 2;
	double dx_w = grid[(i - 1) * m + k].dx;
	double dx = grid[i * m + k].dx;
	double mu = grid[i * m + k].mu;
	double mu_avg_n = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu) / 2;
	double mu_avg = (grid[i * m + k].mu + grid[(i - 1) * m + k].mu +
			grid[i * m + k - 1].mu + grid[(i - 1) * m + k - 1].mu) / 4;
	double mu_w = grid[(i - 1) * m + k].mu;
	double u_e = ti_half_1[k].y[i + 1];
	double u_w = ti_half_1[k].y[i - 1];
	double u_prev = grid_prev[i * m + k].u;
	double u = ti_half_1[k].y[i];
	double p = grid[i * m + k].p;
	double p_w = grid[(i - 1) * m + k].p;
	double w_n = grid[i * m + k + 1].w;
	double w_nw = grid[(i - 1) * m + k + 1].w;
	double w = grid[i * m + k].w;
	double w_w = grid[(i - 1) * m + k].w;
	double delta_t = get_delta_t();
	double density = get_blood_density();

	return (u_prev + delta_t * (2 * (mu * (u_e - u) / dx - mu_w * (u - u_w) / dx_w) / (density * dx_half) -
			(c2 * d2 * u_e + ((1 - c2) * d2 - (1 - c1) * d1) * u - c1 * d1 * u_w) / dx_half -
			(p - p_w) / (density * dx_half) + (mu_avg_n * (w_n - w_nw) - mu_avg * (w - w_w)) / (density * dx_half * dz))) /
			(1 + delta_t * (mu_avg_n / (density * dz) + mu_avg / (density * dz_half_minus) - (1 - c3) * d3) / dz);
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
	calculate_velocity_x_b_half_1(grid_prev, pi.ti_half_1, n, m);
	calculate_velocity_x_c_half_1(grid_prev, pi.ti_half_1, n, m);
	calculate_velocity_x_d_half_1(grid_prev, grid_last_layer, pi.ti_half_1, n, m);
	prs_half_1(&pi);

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
	}

	calculate_velocity_x_a_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_x_b_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_x_c_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_x_d_half_2(grid_prev, grid_last_layer, pi.ti_half_1, pi.ti_half_2, n, m);
	prs_half_2(&pi);

	double *tmp = calloc(n * m, sizeof(double));
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 1; i < pi.n; ++i)
		memcpy(&tmp[i * m], pi.ti_half_2[i].y, m * sizeof(double));

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

			ti[k].a[i] = delta_t * (-mu_avg / (density * dx_half) + c1 * d1) / dx;
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

			ti[k].b[i] = 1 + delta_t * (mu_avg_e / (density * dx_half_plus) + 
				mu_avg / (density * dx_half) + ((1 - c2) * d2 + (1 - c1) * d1)) / dx;
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
			
			ti[k].c[i] = delta_t * (-mu_avg_e / (density * dx_half_plus) + c2 * d2) / dx;
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

	return delta_t * (mu_avg_e / (density * dx_half_plus) - c2 * d2) / 
		(dx + delta_t * (mu_avg_e / (density * dx_half_plus) + 2 * mu_avg / (density * dx) + (1 - c2) * d2));
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

	return delta_t * (mu_avg / (density * dx_half) - c1 * d1) /
			(dx + delta_t * (2 * mu_avg_e / (density * dx) + mu_avg / (density * dx_half) + (1 - c1) * d1));
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
			(2 * (mu * (w_n - w) / dz - mu_s * (w - w_s) / dz_s) / (density * dz_half_minus) -
			(c4 * d4 * w_n + ((1 - c4) * d4 + (1 - c3) * d3) * w + c3 * d3 * w_s) / dz_half_minus -
			(p - p_s) / (density * dz_half_minus) + mu_avg_e * (u_e - u_se) / (density * dx * dz_half_minus)))) / 
			(1 + delta_t * (mu_avg_e / (density * dx_half_plus * dx) + 2 * mu_avg / (density * dx * dx) + (1 - c2) * d2 / dx));
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

	return (w_prev + delta_t * (k == 0 ? (4 * mu * (w_n - w) / (density * dz * dz) - (2 * c3 - 1) * (w_n - w) * d3 / dz) :
			(2 * (mu * (w_n - w) / dz - mu_s * (w - w_s) / dz_s) / (density * dz_half_minus) -
			(c4 * d4 * w_n + ((1 - c4) * d4 + (1 - c3) * d3) * w + c3 * d3 * w_s) / dz_half_minus - 
			(p - p_s) / (dz_half_minus * density) - mu_avg * (u - u_s) / (density * dz_half_minus * dx)))) /
			(1 + delta_t * (2 * mu_avg_e / (density * dx) + mu_avg / (density * dx_half) + (1 - c1) * d1) / dx);
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

			ti[i].a[k] = delta_t * (-2 * mu_s / (density * dz_s) + c3 * d3) / dz_half;
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

			ti[i].b[k] = 1 + delta_t * (2 * mu / (density * dz) + 2 * mu_s / (density * dz_s) + 
					((1 - c4) * d4 + (1 - c3) * d3)) / dz_half_minus;
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

			ti[i].c[k] = delta_t * (-2 * mu / (density * dz) + c4 * d4) / dz_half_minus;
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
			d1 = -(grid[i * m + k].u + grid[i * m + k - 1].u) / 2;
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
				ti_half_2[i].d[k] = w_prev + delta_t * ((mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (2 * w) / dx) / (density * dx) - 
					(c2 * d2 * w_e + (1 - c2) * d2 * w) / dx - (p - p_s) / (density * dz_half_minus) +
					mu_avg_e * (u_e - u_se) / (density * dx * dz_half_minus));
			} else if (i == n - 2) {
				ti_half_2[i].d[k] = w_prev + delta_t * ((-mu_avg_e * (2 * w) / dx - mu_avg * (w - w_w) / dx_half) / (density * dx) -
					(c1 * d1 * w_w + (1 - c1) * d1 * w) / dx -
					(p - p_s) / (density * dz_half_minus) - mu_avg * (u - u_s) / (density * dx * dz_half_minus));
			} else {
				ti_half_2[i].d[k] = w_prev + delta_t * ((mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (w - w_w) / dx_half) / (density * dx) -
					(c2 * d2 * w_e + ((1 - c2) * d2 + (1 - c1) * d1) * w + c1 * d1 * w_w) / dx +
					(mu_avg_e * (u_e - u_se) - mu_avg * (u - u_s)) / (density * dx * dz_half_minus) -
					(p - p_s) / (density * dz_half_minus));
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

	return delta_t * (4 * mu / (density * dz) - (2 * c3 - 1) * d3) /
			(dz + delta_t * (4 * mu / (density * dz) - (2 * c3 - 1) * d3));
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
		return  (w_prev + delta_t * ((mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (2 * w) / dx) / (density * dx) -
				(c2 * w_e + (1 - c2) * w) * d2 / dx)) /
				(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
	} else if (i == n - 2) {
		return  (w_prev + delta_t * ((-mu_avg_e * (2 * w) / dx - mu_avg * (w - w_w) / dx_half_minus) / (density * dx) - 
				(c1 * d1 * w_w + (1 - c1) * d1 * w) / dx)) /
				(1 + delta_t * (4 * mu / (density * dz * dz) - (2 * c3 - 1) * d3 / dz));
	} else {
		return (w_prev + delta_t * ((mu_avg_e * (w_e - w) / dx_half_plus - mu_avg * (w - w_w) / dx_half_minus) / (density * dx) -
				(c2 * d2 * w_e + ((1 - c2) * d2 + (1 - c1) * d1) * w + c1 * d1 * w_w) / dx)) / 
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
	}

	calculate_velocity_z_a_half_1(grid_prev, pi.ti_half_1, n, m);
	calculate_velocity_z_b_half_1(grid_prev, pi.ti_half_1, n, m);
	calculate_velocity_z_c_half_1(grid_prev, pi.ti_half_1, n, m);
	calculate_velocity_z_d_half_1(grid_prev, grid_last_layer, pi.ti_half_1, n, m);
	prs_half_1(&pi);

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
	}

	calculate_velocity_z_a_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_z_b_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_z_c_half_2(grid_prev, pi.ti_half_2, n, m);
	calculate_velocity_z_d_half_2(grid_prev, grid_last_layer, pi.ti_half_1, pi.ti_half_2, n, m);
	prs_half_2(&pi);

	double *tmp = calloc(n * m, sizeof(double));
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif
	for (size_t i = 0; i < pi.n; ++i)
		memcpy(&tmp[i * m], pi.ti_half_2[i].y, m * sizeof(double));

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
		grid[(n - 1) * m + j].u = 0;
	}

	double cur_x = 0.0;
	for (size_t i = 0; i < n; ++i) {
		cur_x += grid[i * m + m - 1].dx / 2;
		grid[i * m + m - 1].u = 0;
		grid[i * m + m - 1].w = -init_w(t, cur_x);
		cur_x += grid[i * m + m - 1].dx / 2;
	}

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
				res += 1;
			}
		}
	}
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
			//calculate_boundary_conditions(grid_cur, n, m, cur_layer * delta_t);
			//calculate_boundary_conditions(grid_prev, n, m, cur_layer * delta_t);
			if (check_convergence(grid_prev, grid_cur, n, m, EPS)) {
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
