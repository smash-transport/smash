/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_DENSITY_H_
#define SRC_INCLUDE_DENSITY_H_

#include <vector>

#include "fourvector.h"
#include "threevector.h"
#include "particledata.h"
#include "pdgcode.h"
#include "forwarddeclarations.h"

namespace Smash {

  /** Allows to choose wich kind of density to calculate.
   *  For Fermi momenta and symmetry potential one needs
   *  to know proton and neutron densities. Baryon density
   *  is necessary for Skyrme potential.
   */
  enum Density_type {baryon, proton, neutron};

  /** A small check if particle PDG code belongs to a given type.
   *  Currently checks for protons, neutrons and baryons.
   *
   *  \param[in] pdg PDG code of particle to be tested
   *  \param[in] dens_type density type: baryon, proton, neutron
   */
  bool particle_in_denstype(const PdgCode pdg, Density_type dens_type);

  /** Calculates 4-current in the computational frame.
   *  \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N C_i u^{\mu}_i
   *  exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
   *  \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
   *  \f[ \rho^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
   *  Here \f$ C_i \f$ is a corresponding value of "charge". If baryon
   *  current option is selected then \f$ C_i \f$ is 1 for baryons,
   *  -1 for antibaryons and 0 otherwize. For proton/neutron
   *  current \f$ C_i = 1\f$ for proton/neutron and 0 otherwize.
   *
   * \param[in] r Arbitrary space point where 4-current is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   * \param[in] dens_type type of four-currect to be calculated:
   *            baryon, proton or neutron options are currently available
   * \param[in] ntest Number of test-particles
   */
  FourVector four_current(const ThreeVector &r, const ParticleList &plist,
                          double gs_sigma, Density_type dens_type, int ntest);

  /** Calculates the gradient of Eckart rest frame density with
   *  respect to computational frame coordinates using analytical formula.
   *  \f[ \frac{d\rho_{Eck}}{d \vec r} = \frac{\frac{dj^{\mu}}{d \vec r}
   *  j_{\mu}}{\sqrt{j^{\mu}j_{\mu}}} \f]
   *  The formula \f$ j^{\mu}(\vec r) \f$ itself is given for the
   *  four_current function. Input parameters are the same that for
   *  four_current function.
   */
  std::pair<double, ThreeVector> rho_eckart_gradient(const ThreeVector &r,
                               const ParticleList &plist, double gs_sigma,
                               Density_type dens_type, int ntest);

  /** Prints density along the specified line. Useful to make 1D plots of
    * density profiles.
   */
  void density_along_line(const char * file_name, const ParticleList &plist,
                        double gs_sigma, Density_type dens_type, int ntest,
                        const ThreeVector &line_start,
                        const ThreeVector &line_end, int n_points);
}  // namespace Smash

#endif  // SRC_INCLUDE_DENSITY_H_
