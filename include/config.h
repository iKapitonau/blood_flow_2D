#ifndef CONFIG_H_
#define CONFIG_H_

#include <stddef.h>

typedef struct stenosis_info_t {
	double pos_x;
	double pos_z;
	double size_x;
	double size_z;
} stenosis_info;

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

void write_config(const char *config_name, const config_data *data);
void read_config(const char *config_name);
void clear_config(void);

// Vessel info
double get_vessel_size_x(void);
double get_vessel_size_z(void);
double get_blood_density(void);

// Stenoses info
void get_stenoses_info(size_t *stenoses_number, stenosis_info **stenoses);

// Grid params
int get_rows_number(void);
int get_columns_number(void);
double get_delta_t(void);

// Initial
void get_velocity_x(size_t *nodes_number, double **velocity_x);
void get_velocity_z(size_t *nodes_number, double **velocity_z);
void get_pressure(size_t *nodes_number, double **pressure);

#endif	// CONFIG_H_
