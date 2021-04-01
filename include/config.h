/**
 * @file config.h
 * @author Ilya Kapitonau <ilya.th.kapitonov@gmail.com>
 * @brief Provides interface for working with config file and data in it.
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <stddef.h>

/**
 * @brief Info about stenosis inside vessel.
 */
typedef struct stenosis_info_t {
	double pos_x1;		///< X-position of bottom left stenosis corner.
	double pos_z1;		///< Z-position of bottom left stenosis corner.
	double pos_x2;		///< X-position of top right stenosis corner.
	double pos_z2;		///< Z-position of top right stenosis corner.
} stenosis_info;

/**
 * @brief Struct that contains all data stored in config file.
 */
typedef struct config_data_t {
	// Vessel info
	double vessel_size_x;
	double vessel_size_z;
	double blood_density;
	
	// Stenoses info
	size_t stenoses_number;
	stenosis_info *stenoses;
	
	// Grid params
	int rows_number;
	int columns_number;
	double delta_t;

	// Initial
	size_t nodes_number;
	double *velocity_x;
	double *velocity_z;
	double *pressure;
} config_data;

/**
 * @brief Writes data to config file.
 *
 * This function writes data from @p data structure to config file with name 
 * @p config_name.
 *
 * @param config_name Config file name.
 * @param data Data that will be written to config file.
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE if error has occured.
 */
int write_config(const char *config_name, const config_data *data);

/**
 * @brief Reads data from config file.
 *
 * This function reads data from config file with name @p config_name and
 * stores data from it inside static @c config_data structure.
 *
 * @param config_name Config file name.
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE if error has occured.
 */
int read_config(const char *config_name);

/**
 * @brief Frees memory allocated for static @c config_data structure.
 */
void clear_config(void);

/**
 * @brief Returns vessel size along X-axis.
 * @return Value of vessel size along X-axis.
 */
double get_vessel_size_x(void);

/**
 * @brief Returns vessel size along Z-axis.
 * @return Value of vessel size along Z-axis.
 */
double get_vessel_size_z(void);

/**
 * @brief Returns blood pressure in vessel.
 * @return Value of blood pressure in vessel.
 */
double get_blood_density(void);

/**
 * @brief Returns info about stenoses in vessel.
 * @param stenoses_number Pointer to variable where number of stenoses will be written.
 * @param stenoses Pointer to @c stenosis_info struct where info about stenoses will be allocated and written.
 * @note Memory for @p stenoses is allocated inside the function. It is up to user to free the memory!
 */
void get_stenoses_info(size_t *stenoses_number, stenosis_info **stenoses);

/**
 * @brief Returns number of rows in computational grid.
 * @return Value of rows number.
 */
int get_rows_number(void);

/**
 * @brief Returns number of columns in computational grid.
 * @return Value of columns number.
 */
int get_columns_number(void);

/**
 * @brief Returns time step used in numerical algorithms.
 * @return Value of time step.
 */
double get_delta_t(void);

/**
 * @brief Returns info about blood velocity along X-axis.
 * @param nodes_number Pointer to variable where length of @p velocity_x will be written.
 * @param velocity_x Pointer where blood velocity along X-axis will be allocated and written.
 * @note Memory for @p velocity_x is allocated inside the function. It is up to user to free the memory!
 * @note Info in @p velocity_x is actually a 2D array squashed into 1D. So to get dimensions of the
 * original array, you should call @c get_rows_number and @c get_columns_number functions.
 */
void get_velocity_x(size_t *nodes_number, double **velocity_x);

/**
 * @brief Returns info about blood velocity along Z-axis.
 * @param nodes_number Pointer to variable where length of @p velocity_z will be written.
 * @param velocity_z Pointer where blood velocity along Z-axis will be allocated and written.
 * @note Memory for @p velocity_z is allocated inside the function. It is up to user to free the memory!
 * @note Info in @p velocity_z is actually a 2D array squashed into 1D. So to get dimensions of the
 * original array, you should call @c get_rows_number and @c get_columns_number functions.
 */
void get_velocity_z(size_t *nodes_number, double **velocity_z);

/**
 * @brief Returns info about blood pressure.
 * @param nodes_number Pointer to variable where length of @p pressure will be written.
 * @param velocity_x Pointer where blood pressure info will be allocated and written.
 * @note Memory for @p pressure is allocated inside the function. It is up to user to free the memory!
 * @note Info in @p pressure is actually a 2D array squashed into 1D. So to get dimensions of the
 * original array, you should call @c get_rows_number and @c get_columns_number functions.
 */
void get_pressure(size_t *nodes_number, double **pressure);

#endif	// CONFIG_H_
