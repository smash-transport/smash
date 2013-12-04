/*
 *
 *    Copyright (c) 2012
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

static void assign_params(std::list<Parameters> *configuration) {
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("Looking for match %s %s\n", key, value);
        
        /* integer values */
        if (strcmp(key, "STEPS") == 0) {
            parameters->set_steps(abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "RANDOMSEED") == 0) {
            /* negative seed means random startup value */
            if (atol(value) > 0)
                parameters->set_seed(atol(value));
            else
                parameters->set_seed(time(NULL));
            match = true;
        }
        if (strcmp(key, "UPDATE") == 0) {
            parameters->set_output_interval(abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "TESTPARTICLES") == 0) {
            parameters->set_testparticles(abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "MODUS") == 0) {
            parameters->set_modus(abs(atoi(value)));
            match = true;
        }
        
        
        /* double or float values */
        if (strcmp(key, "EPS") == 0) {
            parameters->set_eps(fabs(atof(value)));
            match = true;
        }
        if (strcmp(key, "SIGMA") == 0) {
            parameters->set_cross_section(fabs(atof(value)));
            match = true;
        }
        
        /* remove processed entry */
        if (match) {
            printd("Erasing %s %s\n", key, value);
            i = configuration->erase(i);
            match = false;
        } else {
            ++i;
        }
    }
}

static void assign_params(std::list<Parameters> *configuration, Box *box) {
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("%s %s\n", key, value);
        
        /* double or float values */
        if (strcmp(key, "LENGTH") == 0) {
            box->set_length(fabs(atof(value)));
            match = true;
        }
        if (strcmp(key, "TEMPERATURE") == 0) {
            box->set_temperature(fabs(atof(value)));
            match = true;
        }
        /* int values */
        if (strcmp(key, "INITIAL_CONDITION") == 0) {
            box->set_initial_condition(abs(atoi(value)));
            match = true;
        }
        /* remove processed entry */
        if (match) {
            i = configuration->erase(i);
            match = false;
        } else {
            ++i;
        }
    }
}

static void assign_params(std::list<Parameters> *configuration, Sphere *ball) {
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("%s %s\n", key, value);
        
        /* double or float values */
        if (strcmp(key, "RADIUS") == 0) {
            ball->set_radius(fabs(atof(value)));
            match = true;
        }
        
        /* remove processed entry */
        if (match) {
            i = configuration->erase(i);
            match = false;
        } else {
            ++i;
        }
    }
}



/* process_box_config - configuration handling */
void Box::process_config(char *path) {
    std::list<Parameters> configuration;
    size_t len = strlen("./config_box.txt") + strlen(path) + 1;
    char *config_path = reinterpret_cast<char *>(malloc(len));
    snprintf(config_path, len, "%s/config_box.txt", path);
    process_params(config_path, &configuration);
    assign_params(&configuration, this);
    warn_wrong_params(&configuration);
    free(config_path);
}

/* process_sphere_config - configuration handling */
void Sphere::process_config(char *path) {
    std::list<Parameters> configuration;
    size_t len = strlen("./config_sphere.txt") + strlen(path) + 1;
    char *config_path = reinterpret_cast<char *>(malloc(len));
    snprintf(config_path, len, "%s/config_sphere.txt", path);
    process_params(config_path, &configuration);
    assign_params(&configuration, this);
    warn_wrong_params(&configuration);
    free(config_path);
}
