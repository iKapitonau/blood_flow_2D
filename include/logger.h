#ifndef LOGGER_H_
#define LOGGER_H_

// Log levels

typedef enum Log_lvl_t {
	LOG_EMERG = 0,	// 0
	LOG_ALERT,		// 1
	LOG_CRIT,		// 2
	LOG_ERR,		// 3
	LOG_WARNING,	// 4
	LOG_NOTICE,		// 5
	LOG_INFO,		// 6
	LOG_DEBUG,		// 7
} Log_lvl;

static const char *strloglvl[] = {
	"EMERG",
	"ALERT",
	"CRIT",
	"ERROR",
	"WARN",
	"NOTICE",
	"INFO",
	"DEBUG"
};

// Log options

typedef enum Log_option_t {
	LOG_FILE	= 1,		// write to log_file
	LOG_STDERR	= 2,		// write to stderr
} Log_option;

int log_open(const char *log_name, Log_option mode);
void log_write(Log_lvl lvl, const char *fmt, ...);
int log_close();

#endif	// LOGGER_H_
