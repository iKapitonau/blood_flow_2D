#include "logger.h"
#include "config.h"
#include "calc.h"

#include <stdio.h>
#include <omp.h>

int main(void)
{
	//log_open(NULL, LOG_STDERR);

	read_config("config.ex1");

	double t1 = omp_get_wtime();
	calculate(5, 0, 1);
	double t2 = omp_get_wtime();

	fprintf(stderr, "TOTAL=%lf\n", t2 - t1);

	clear_config();
}
