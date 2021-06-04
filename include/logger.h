/**
 * @file logger.h
 * @author Ilya Kapitonau <ilya.th.kapitonov@gmail.com>
 * @brief Provides basic interface for logging.
 */

#ifndef LOGGER_H_
#define LOGGER_H_

/**
 * @brief Logging levels
 */
typedef enum Log_lvl_t {
	LOG_EMERG = 0,	///< 0
	LOG_ALERT,		///< 1
	LOG_CRIT,		///< 2
	LOG_ERR,		///< 3
	LOG_WARNING,	///< 4
	LOG_NOTICE,		///< 5
	LOG_INFO,		///< 6
	LOG_DEBUG,		///< 7
} Log_lvl;

/**
 * @brief Logger options enum
 *
 * Logger options specify where logger will write:
 * to log_file or/and to stderr.
 */
typedef enum Log_option_t {
	LOG_FILE	= 1,		///< write to log_file
	LOG_STDERR	= 2,		///< write to stderr
} Log_option;

/**
 * @brief Opens log file.
 *
 * This function is used to open log file with name @p log_name if @c LOG_FILE
 * option is set in @p mode argument. It MUST be called before the first
 * logging.
 *
 * @param log_name Name of log file or @c NULL if @c LOG_FILE is not set in @p mode.
 * @param mode Specifies whether logs will be written to file or/and to stderr.
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE if error has occured.
 * @warning If @p mode is @c NULL then logs will not be written anywhere.
 */
int log_open(const char *log_name, Log_option mode);

/**
 * @brief Writes to opened log file.
 *
 * This function writes log message with level @p lvl in @c printf format to
 * the log file opened by @c log_open function.
 * 
 * @param lvl Logging level of the message.
 * @param fmt Format string in printf style.
 */
void log_write(Log_lvl lvl, const char *fmt, ...);

/**
 * @brief Closes open log file.
 * 
 * This function closes log file if any is opened by @c log_open function.
 *
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE if error has occured.
 */
int log_close(void);

#endif	// LOGGER_H_
