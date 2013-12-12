/*
 *
 *    Copyright (c) 2012
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include "include/BoundaryConditions.h"
#include "include/outputroutines.h"

void BoundaryConditions::assign_params_general(std::list<Parameters> *configuration) {
    bool match = false;
    std::list<Parameters>::iterator i = configuration->begin();
    while (i != configuration->end()) {
        char *key = i->key();
        char *value = i->value();
        printd("Looking for match %s %s\n", key, value);
        
        /* integer values */
        if (strcmp(key, "STEPS") == 0) {
            steps = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "RANDOMSEED") == 0) {
            /* negative seed means random startup value */
            if (atol(value) > 0)
                seed = (atol(value));
            else
                seed = (time(NULL));
            match = true;
        }
        if (strcmp(key, "UPDATE") == 0) {
            output_interval = (abs(atoi(value)));
            match = true;
        }
        if (strcmp(key, "TESTPARTICLES") == 0) {
            testparticles = (abs(atoi(value)));
            match = true;
        }
        
        /* double or float values */
        if (strcmp(key, "EPS") == 0) {
            eps = (fabs(atof(value)));
            match = true;
        }
        if (strcmp(key, "SIGMA") == 0) {
            cross_section = (fabs(atof(value)));
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

/* print_startup - console output on startup of general parameters */
void BoundaryConditions::print_startup() {
    printf("Elastic cross section: %g [mb]\n", cross_section);
    printf("Using temporal stepsize: %g [GeV]\n", eps);
    printf("Maximum number of steps: %i \n", steps);
    printf("Random number seed: %lli \n", seed);
}




