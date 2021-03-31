#include "config.h"
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>

static config_data data_;

int write_config(const char *config_name, const config_data *data)
{
	log_write(LOG_INFO, "Config: writing to config '%s'...", config_name);

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
		config_setting_set_float(setting, data->stenoses[i].pos_x1);
		setting = config_setting_add(subgroup, "stenosis_pos_z", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].pos_z1);
		setting = config_setting_add(subgroup, "stenosis_size_x", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].pos_x2);
		setting = config_setting_add(subgroup, "stenosis_size_z", CONFIG_TYPE_FLOAT);
		config_setting_set_float(setting, data->stenoses[i].pos_z2);
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
		return EXIT_FAILURE;
	}

	config_destroy(&cfg);

	log_write(LOG_INFO, "Config: writing done!");

	return EXIT_SUCCESS;
}

int read_config(const char *config_name)
{
	log_write(LOG_INFO, "Config: Reading config '%s'...", config_name);

	config_t cfg;
	config_setting_t *setting;
	int ret;

	config_init(&cfg);

	if (config_read_file(&cfg, config_name) == CONFIG_FALSE) {
		log_write(LOG_ERR, "Config: %s:%d %s", config_error_file(&cfg),
			config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return EXIT_FAILURE;
	}

	// Vessel info
	ret = config_lookup_float(&cfg, "vessel_params.vessel_size_x", &data_.vessel_size_x);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'vessel_params.vessel_size_x' was not read. "
				"Missing or with wrong type (must be double)");

	ret = config_lookup_float(&cfg, "vessel_params.vessel_size_z", &data_.vessel_size_z);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'vessel_params.vessel_size_z' was not read. "
				"Missing or with wrong type (must be double)");

	config_lookup_float(&cfg, "vessel_params.blood_density", &data_.blood_density);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'vessel_params.blood_density' was not read. "
				"Missing or with wrong type (must be double)");
	
	// Stenoses info
	setting = config_lookup(&cfg, "vessel_params.stenoses_params");
	if (!setting) {
		log_write(LOG_WARNING, "Config: group 'vessel_params.stenoses_params' was not read (missing?).");
	} else {
		data_.stenoses_number = config_setting_length(setting);
		if (!data_.stenoses_number)
			log_write(LOG_WARNING, "Config: group 'vessel_params.stenoses_params' is empty.");
		else 
			data_.stenoses = malloc(data_.stenoses_number * sizeof(stenosis_info));

		for (size_t i = 0; i < data_.stenoses_number; ++i) {
			config_setting_t *cur_elem = config_setting_get_elem(setting, i);
			ret = config_setting_lookup_float(cur_elem, "stenosis_pos_x1", &data_.stenoses[i].pos_x1);
			if (ret == CONFIG_FALSE)
				log_write(LOG_WARNING, "Config: parameter 'vessel_params.stenoses[%zu].params.stenosis_pos_x1' was not read. "
						"Missing or with wrong type (must be double)", i);

			ret = config_setting_lookup_float(cur_elem, "stenosis_pos_z1", &data_.stenoses[i].pos_z1);
			if (ret == CONFIG_FALSE)
				log_write(LOG_WARNING, "Config: parameter 'vessel_params.stenoses[%zu].params.stenosis_pos_z1' was not read. "
						"Missing or with wrong type (must be double)", i);

			ret = config_setting_lookup_float(cur_elem, "stenosis_pos_x2", &data_.stenoses[i].pos_x2);
			if (ret == CONFIG_FALSE)
				log_write(LOG_WARNING, "Config: parameter 'vessel_params.stenoses[%zu].params.stenosis_pos_x2' was not read. "
						"Missing or with wrong type (must be double)", i);

			ret = config_setting_lookup_float(cur_elem, "stenosis_pos_z2", &data_.stenoses[i].pos_z2);
			if (ret == CONFIG_FALSE)
				log_write(LOG_WARNING, "Config: parameter 'vessel_params.stenoses[%zu].params.stenosis_pos_z2' was not read. "
						"Missing or with wrong type (must be double)", i);
		}
	}

	// Grid params
	ret = config_lookup_int(&cfg, "grid_params.rows_number", &data_.rows_number);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'grid_params.rows_number' was not read. "
				"Missing or with wrong type (must be int)");

	ret = config_lookup_int(&cfg, "grid_params.columns_number", &data_.columns_number);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'grid_params.columns_number' was not read. "
				"Missing or with wrong type (must be int)");

	ret = config_lookup_float(&cfg, "grid_params.delta_t", &data_.delta_t);
	if (ret == CONFIG_FALSE)
		log_write(LOG_WARNING, "Config: parameter 'grid_params.delta_t' was not read. "
				"Missing or with wrong type (must be double)");

	// Initial
	data_.nodes_number = (data_.rows_number + 1) * (data_.columns_number + 1);
	if (!data_.rows_number || !data_.columns_number) {
		data_.nodes_number = 0;
		log_write(LOG_WARNING, "Config: rows_number or columns_number is zero.");
	} else {
		data_.velocity_x = malloc(data_.nodes_number * sizeof(double));
		data_.velocity_z = malloc(data_.nodes_number * sizeof(double));
		data_.pressure = malloc(data_.nodes_number * sizeof(double));
	}

	config_setting_t *setting_velocity_x = config_lookup(&cfg, "initial.velocity_x");
	if (!setting_velocity_x) {
		log_write(LOG_WARNING, "Config: parameter 'initial.velocity_x' was not read. "
				"Missing or with wrong type (must be array)");
	} else {
		for (size_t i = 0; i < data_.nodes_number; ++i)
			data_.velocity_x[i] = config_setting_get_float_elem(setting_velocity_x, i);
	}

	config_setting_t *setting_velocity_z = config_lookup(&cfg, "initial.velocity_z");
	if (!setting_velocity_z) {
		log_write(LOG_WARNING, "Config: parameter 'initial.velocity_z' was not read. "
				"Missing or with wrong type (must be array)");
	} else {
		for (size_t i = 0; i < data_.nodes_number; ++i)
			data_.velocity_z[i] = config_setting_get_float_elem(setting_velocity_z, i);
	}

	config_setting_t *setting_pressure = config_lookup(&cfg, "initial.pressure");
	if (!setting_pressure) {
		log_write(LOG_WARNING, "Config: parameter 'initial.pressure' was not read. "
				"Missing or with wrong type (must be array)");
	} else {
		for (size_t i = 0; i < data_.nodes_number; ++i)
			data_.pressure[i] = config_setting_get_float_elem(setting_pressure, i);
	}

	config_destroy(&cfg);

	log_write(LOG_INFO, "Config: Reading done!");

	return EXIT_SUCCESS;
}

void clear_config(void)
{
	free(data_.velocity_x);
	free(data_.velocity_z);
	free(data_.pressure);
	free(data_.stenoses);
}

// Vessel info
double get_vessel_size_x(void)
{
	return data_.vessel_size_x;
}

double get_vessel_size_z(void)
{
	return data_.vessel_size_z;
}

double get_blood_density(void)
{
	return data_.blood_density;
}

// Stenoses info
void get_stenoses_info(size_t *stenoses_number, stenosis_info **stenoses)
{
	*stenoses_number = data_.stenoses_number;
	*stenoses = malloc(*stenoses_number * sizeof(stenosis_info));
	memcpy(*stenoses, data_.stenoses, *stenoses_number * sizeof(stenosis_info));
}

// Grid params
int get_rows_number(void)
{
	return data_.rows_number;
}

int get_columns_number(void)
{
	return data_.columns_number;
}

double get_delta_t(void)
{
	return data_.delta_t;
}

// Initial
void get_velocity_x(size_t *nodes_number, double **velocity_x)
{
	*nodes_number = data_.nodes_number;
	*velocity_x = malloc(*nodes_number * sizeof(double));
	memcpy(*velocity_x, data_.velocity_x, *nodes_number * sizeof(double));
}

void get_velocity_z(size_t *nodes_number, double **velocity_z)
{
	*nodes_number = data_.nodes_number;
	*velocity_z = malloc(*nodes_number * sizeof(double));
	memcpy(*velocity_z, data_.velocity_z, *nodes_number * sizeof(double));
}

void get_pressure(size_t *nodes_number, double **pressure)
{
	*nodes_number = data_.nodes_number;
	*pressure = malloc(*nodes_number * sizeof(double));
	memcpy(*pressure, data_.pressure, *nodes_number * sizeof(double));
}
