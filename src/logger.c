/**
 * @file logger.c
 * @author Ilya Kapitonau <ilya.th.kaptionov@gmail.com>
 * @brief Implementation of logging functions.
 */

#include "logger.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF_SIZE		256		///< Buffer size for log message.
#define TIMESTAMP_SIZE	25		///< Buffer size for timestamp.

static FILE *log_file;
static Log_option log_mode;

int log_open(const char *log_name, Log_option mode)
{
	log_mode = mode;

	if (log_mode & LOG_FILE) {
		log_file = fopen(log_name, "a");

		if (!log_file)
			perror("fopen");
			return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

void log_write(Log_lvl lvl, const char *fmt, ...)
{
	/**
	 * Get time of log writing.
	 */
	time_t rawtime = time(NULL);

	if (rawtime == (time_t)-1) {
		perror("time");
		return;
	}

	struct tm *timestamp = localtime(&rawtime);

	if (!timestamp) {
		perror("localtime");
		return;
	}

	char strtime[TIMESTAMP_SIZE];	///< String representation of timestamp.

	size_t time_len = strftime(strtime, TIMESTAMP_SIZE, "%d-%m-%y %r", timestamp);

	char buf[BUF_SIZE];
	sprintf(buf, "[%s][%s]: %s\n", strtime, strloglvl[lvl], fmt);

	va_list vl;
	if (log_mode & LOG_FILE) {
		va_start(vl, buf);
		vfprintf(log_file, buf, vl);
		va_end(vl);
	}
	if (log_mode & LOG_STDERR) {
		va_start(vl, buf);
		vfprintf(stderr, buf, vl);
		va_end(vl);
	}
}

int log_close(void)
{
	if (log_mode & LOG_FILE) {
		int ret = fclose(log_file);

		if (ret == EOF)
			return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
