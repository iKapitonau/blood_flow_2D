#include "logger.h"
#include "config.h"
#include "calc.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

static const char *std_config_filename = "config_example";

void print_usage(void)
{
	fprintf(stderr, "Usage: blood_flow_2d [config_filename]\n");
}

int main(int argc, char *argv[])
{
	static const char *config_filename;
	if (argc > 2) {
		print_usage();
		return EXIT_SUCCESS;
	} else if (argc == 1) {
		config_filename = std_config_filename;
	} else {
		config_filename = argv[1];
	}

	log_open(NULL, LOG_STDERR);

	int ret = read_config(config_filename);
	if (ret != EXIT_SUCCESS) {
		log_write(LOG_CRIT, "Error while reading config file '%s'. Exiting...", config_filename);
		exit(EXIT_FAILURE);
	}

	/* default params */
	double sigma = 5;
	double t_from = 0;
	double t_to = 1;

	log_write(LOG_INFO, "Calculating with SIGMA = %lf, TIME_RANGE = [%lf, %lf]", sigma, t_from, t_to);

	double t1 = omp_get_wtime();
	calculate(sigma, t_from, t_to);
	double t2 = omp_get_wtime();

	log_write(LOG_INFO, "TIME=%lf", t2 - t1);

	clear_config();

	return EXIT_SUCCESS;
}
