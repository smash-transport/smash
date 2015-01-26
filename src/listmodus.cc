/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>


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
     * \key Start_Time (float, required):\n
     * Starting time of List calculation.
     *
     * \key Shift_Id (int, required):\n
     * Starting id for event_id_
     *
     * Example input in particle_lists_in/event{id}
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
        :start_time_(modus_config.take({"List", "Start_Time"})),
        shift_id_(modus_config.take({"List", "Shift_Id"})) {
            event_id_ = shift_id_;
        }


    /* console output on startup of List specific parameters */
    std::ostream &operator<<(std::ostream &out, const ListModus &m) {
        out << "\nStarting time for List calculation: "
            << m.start_time_ << '\n';
        return out;
    }

    /* initial_conditions - sets particle data for @particles */
    float ListModus::initial_conditions(Particles *particles,
            const ExperimentParameters &parameters) {
        const auto &log = logger<LogArea::List>();

        /* Readin PARTICLES from file */

        std::stringstream fname;
        fname << "event" << event_id_;

        const bf::path default_path = bf::absolute("particle_lists_in");
        log.debug() << default_path.filename().native() << '\n';

        const bf::path fpath = default_path/fname.str();
        if ( !bf::exists(fpath) ) {
            log.fatal() << fpath.filename().native()
                << " does not exist! \n"
                <<"\n Usage of smash with external particle lists:\n"
                <<"1. Create particle_lists_in/ directory in build/\n"
                <<"2. Put the external particle lists in file \n"
                <<"particle_lists_in/event{id} where {id} traversal [0,Nevent-1]\n"
                <<"3. Each row should contains t x y z mass p0 px py pz pdg ID \n"
                <<"in units of: fm fm fm fm GeV GeV GeV GeV GeV none none \n";

            exit(EXIT_FAILURE);
        }

        bf::ifstream fin(fpath);

        std::string line, pdg;
        float t, x, y, z, mass, E, px, py, pz, id;


        while ( fin.good() ) {
            std::getline(fin, line);
            if ( fin.eof() ) {
                break;
            }
            if ( line.length() != 0 && line[0] !='#' ) {
                /** \todo test the data format of particle list*/
                std::stringstream(line) >> t >> x >> y >> z
                >> mass >> E >> px >> py >> pz >> pdg >> id;
                log.debug() << "Particle " << pdg << " (mass,x,y,z)= " << "( ";
                log.debug() << mass << "," << x << "," << y << "," << z << ")\n";

                ParticleData particle(ParticleType::find(PdgCode(pdg)));
                particle.set_4momentum(FourVector(E, px, py, pz));
                particle.set_4position(FourVector(start_time_, x, y, z));
                particles->add_data(particle);
            }
        }
        fin.close();

        event_id_++;

        return start_time_;
    }
}  // namespace Smash
