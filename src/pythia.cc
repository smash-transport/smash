/*
 *
 *    Copyright (c) 2015-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *    If Pythia cite 
 *    T. Sjöstrand, S. Mrenna and P. Skands, JHEP05 (2006) 026,
 *                          Comput. Phys. Comm. 178 (2008) 852.
 *    
 */

#include "include/pythia.h"

#include "include/forwarddeclarations.h"
#include "include/logging.h"
#include "include/particledata.h"

///#ifdef Pythia_FOUND
#include "Pythia8/Pythia.h"
///#include "Pythia8/LHAPDFInterface.h"
///#endif

namespace Smash {  
  /// This function will generate outgoing particles in CM frame from a hard process
  ParticleList string_excitation(const ParticleList &incoming_particles_) {
	const auto &log = logger<LogArea::ScatterAction>();  
///    #ifdef Pythia_FOUND
	  /// set all necessary parameters for Pythia call
	  /// Call Pythia
      std::string xmlpath = PYTHIA_XML_DIR;
      Pythia8::Pythia pythia( xmlpath, false ); 
///    #else
///      std::string errMsg = "Pythia 8 not available for string excitation";
///      throw std::runtime_error( errMsg );
///    #endif     
    ParticleList outgoing_particles_; 
    return outgoing_particles_; 
  }	
}
