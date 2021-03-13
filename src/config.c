#include "config.h"
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h>

// Vessel info
static double vessel_size_x;
static double vessel_size_z;
static double blood_density;

// Stenoses info
static size_t stenoses_number;
static stenosis_info *stenoses;

// Grid params
static int rows_number;
static int columns_number;
static double delta_t;

// Initial
static double *velocity_x;
static double *velocity_z;
static double *pressure;

void write_config(const char *config_name, const config_data *data)
{
	log_write(LOG_DEBUG, "Entering write_config function");

	config_t cfg;
	config_setting_t *root, *group, *setting, *array, *subgroup, *list;

	config_init(&cfg);
	root = config_root_setting(&cfg);

	// Vessel info
	group = config_setting_add(root, "vessel_params", CONFIG_TYPE_GROUP);

	setting = config_setting_add(group, "vessel_size_x", CONFIG_TYPE_FLOAT);
	config_setting_set_float(setting, data->vessel_size_x);

	setting = config_setting_add(group, "vessel_size_z", CONFIG_TYPE_FLOAT);
	config_setting_set_float(setting, data->vessel_size_z);

	setting = config_setting_add(group, "blood_density", CONFIG_TYPE_FLOAT);
	config_setting_set_float(setting, data->blood_density);

	// Stenoses info
	list = config_setting_add(group, "stenoses_params", CONFIG_TYPE_LIST);
	for (size_t i = 0; i < data->stenoses_number; ++i) {
		subgroup = config_setting_add(list, NULL, CONFIG_TYPE_GROUP);
		setting = config_setting_add(subgroup, "stenosis_pos_x", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].pos_x);
		setting = config_setting_add(subgroup, "stenosis_pos_z", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].pos_z);
		setting = config_setting_add(subgroup, "stenosis_size_x", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].size_x);
		setting = config_setting_add(subgroup, "stenosis_size_z", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].size_z);
	}

	// Grid params
	group = config_setting_add(root, "grid_params", CONFIG_TYPE_GROUP);

	setting = config_setting_add(group, "rows_number", CONFIG_TYPE_INT);
	config_setting_set_int(setting, data->rows_number);

	setting = config_setting_add(group, "columns_number", CONFIG_TYPE_INT);
	config_setting_set_int(setting, data->columns_number);

	setting = config_setting_add(group, "delta_t", CONFIG_TYPE_FLOAT);
	config_setting_set_float(setting, data->delta_t);

	// Initial params
	group = config_setting_add(root, "initial", CONFIG_TYPE_GROUP);

	array = config_setting_add(group, "velocity_x", CONFIG_TYPE_ARRAY);
	for (size_t i = 0; i < data->nodes_number; ++i) {
		setting = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->velocity_x[i]);
	}

	array = config_setting_add(group, "velocity_z", CONFIG_TYPE_ARRAY);
	for (size_t i = 0; i < data->nodes_number; ++i) {
		setting = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->velocity_z[i]);
	}

	array = config_setting_add(group, "pressure", CONFIG_TYPE_ARRAY);
	for (size_t i = 0; i < data->nodes_number; ++i) {
		setting = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->pressure[i]);
	}

	if (config_write_file(&cfg, config_name) == CONFIG_FALSE) {
		log_write(LOG_ERR, "Config: Error while writing file.");
		config_destroy(&cfg);
		return;
	}

	config_destroy(&cfg);

	log_write(LOG_DEBUG, "Returning from write_config function");
}

void read_config(const char *config_name)
{
	log_write(LOG_DEBUG, "Entering read_config function");

	config_t cfg;
	config_setting_t *setting;

	config_init(&cfg);

	if (config_read_file(&cfg, config_name) == CONFIG_FALSE) {
		log_write(LOG_ERR, "Config: %s:%d %s", config_error_file(&cfg),
			config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return;
	}

	// Vessel info
	config_lookup_float(&cfg, "vessel_params.vessel_size_x", &vessel_size_x);
	config_lookup_float(&cfg, "vessel_params.vessel_size_z", &vessel_size_z);
	config_lookup_float(&cfg, "vessel_params.blood_density", &blood_density);
	
	// Stenoses info
	setting = config_lookup(&cfg, "vessel_params.stenoses_params");
	stenoses_number = config_setting_length(setting);
	stenoses = malloc(stenoses_number * sizeof(stenosis_info));
	for (size_t i = 0; i < stenoses_number; ++i) {
		config_setting_t *cur_elem = config_setting_get_elem(setting, i);
		config_setting_lookup_float(cur_elem, "stenosis_pos_x", &stenoses[i].pos_x);
		config_setting_lookup_float(cur_elem, "stenosis_pos_z", &stenoses[i].pos_z);
		config_setting_lookup_float(cur_elem, "stenosis_size_x", &stenoses[i].size_x);
		config_setting_lookup_float(cur_elem, "stenosis_size_z", &stenoses[i].size_z);
	}

	// Grid params
	config_lookup_int(&cfg, "grid_params.rows_number", &rows_number);
	config_lookup_int(&cfg, "grid_params.columns_number", &columns_number);
	config_lookup_float(&cfg, "grid_params.delta_t", &delta_t);

	// Initial
	size_t nodes_number = rows_number * columns_number;
	velocity_x = malloc(nodes_number * sizeof(double));
	velocity_z = malloc(nodes_number * sizeof(double));
	pressure = malloc(nodes_number * sizeof(double));

	config_setting_t *setting_velocity_x = config_lookup(&cfg, "initial.velocity_x");
	config_setting_t *setting_velocity_z = config_lookup(&cfg, "initial.velocity_z");
	config_setting_t *setting_pressure = config_lookup(&cfg, "initial.pressure");
	for (size_t i = 0; i < nodes_number; ++i) {
		velocity_x[i] = config_setting_get_float_elem(setting_velocity_x, i);
		velocity_z[i] = config_setting_get_float_elem(setting_velocity_z, i);
		pressure[i] = config_setting_get_float_elem(setting_pressure, i);
	}

	config_destroy(&cfg);

	log_write(LOG_DEBUG, "Returning from read_config function");
}

void clear_config(void)
{
	log_write(LOG_DEBUG, "Entering clear_config function");

	free(velocity_x);
	free(velocity_z);
	free(pressure);
	free(stenoses);

	log_write(LOG_DEBUG, "Returning from clear_config function");
}

double get_vessel_size_x(void)
{
	return vessel_size_x;
}

double get_vessel_size_z(void)
{
	return vessel_size_z;
}

double get_blood_density(void)
{
	return blood_density;
}

void get_stenoses_info(size_t *stenoses_n, stenosis_info **stenoses_arr)
{
	*stenoses_n = stenoses_number;
	*stenoses_arr = stenoses;
}

int get_rows_number(void)
{
	return rows_number;
}

int get_columns_number(void)
{
	return columns_number;
}

double get_delta_t(void)
{
	return delta_t;
}

void get_velocity_x(size_t *nodes_number, double **vel_x)
{
	*nodes_number = rows_number * columns_number;
	*vel_x = velocity_x;
}

void get_velocity_z(size_t *nodes_number, double **vel_z)
{
	*nodes_number = rows_number * columns_number;
	*vel_z = velocity_z;
}

void get_pressure(size_t *nodes_number, double **pres)
{
	*nodes_number = rows_number * columns_number;
	*pres = pressure;
}
