/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include <fstream>

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/constants.h"
#include "include/configuration.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/listmodus.h"
#include "include/threevector.h"

namespace Smash {

    /*!\Userguide
     * \page input_modi_list_ List
     *
     * \key File_Directory (string, required): \n
     * The directory of the files that store external particle list
     *
     * \key File_Prefix (string, required): \n
     * The prefix of the files that store external particle list
     *
     * \key Start_Time (float, required):\n
     * Starting time of List calculation.
     *
     * Example data format of particle list in the input file:
     * \verbatim
     #!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID
     # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
     0.1 6.42036 1.66473 9.38499 0.138 0.232871 0.116953 -0.115553 0.0903033 111 0
     \endverbatim
     * It means that one pi^{0} at spatial coordinates (t,x,y,z)=(0.1, 6.42036, 1.66473, 9.38499) fm 
     * and 4-momenta (p0, px. py, pz)=(0.232871, 0.116953, -0.115553, 0.0903033 ) GeV, 
     * with mass=0.138, pdg=111 and id=0 will be initialized.
     */


    ListModus::ListModus(Configuration modus_config,
            const ExperimentParameters &)
        :particle_list_file_directory_(modus_config.take({"List", "File_Directory"})),
        particle_list_file_prefix_(modus_config.take({"List", "File_Prefix"})),
        start_time_(modus_config.take({"List", "Start_Time"}))
        {
        }

    /* console output on startup of List specific parameters */
    std::ostream &operator<<(std::ostream &out, const ListModus &m) {
        //out << "-- List Modus:\nPath of the partile list file: " << m.fpath_particle_list_
        out << "\nStarting time for List calculation: " << m.start_time_ << '\n';
        out << "\nFile path for List calculation: " << m.particle_list_file_directory_<< '\n';
        return out;
    }


    /* initial_conditions - sets particle data for @particles */
    float ListModus::initial_conditions(Particles *particles,
            const ExperimentParameters &parameters) {
        const auto &log = logger<LogArea::List>();

        /* Readin PARTICLES from file */
        std::ifstream fin("particle_lists.oscar");

        std::string line, pdg;
        float t, x, y, z, mass, E, px, py, pz, id;

        while( fin.good() ){
            std::getline(fin, line);
            if( fin.eof() ){
                break;
            }
            if( line.length() != 0 && line[0] !='#' ){
                /** \todo test the data format of particle list*/
                std::stringstream(line)>>t>>x>>y>>z>>mass>>E>>px>>py>>pz>>pdg>>id;
                log.debug() << "Particle " << pdg << " (mass,x,y,z)= " << "( ";
                log.debug() << mass <<","<<x<<","<<y<<","<<z<<")\n";

                ParticleData particle(ParticleType::find(PdgCode(pdg)));
                particle.set_4momentum(FourVector(E,px,py,pz));
                particle.set_4position(FourVector(start_time_,x,y,z));
                particles->add_data( particle );
            }
        }
        fin.close();

        return start_time_;
    }
}  // namespace Smash
