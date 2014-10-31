/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <fstream>
#include <iomanip>
#include "include/constants.h"
#include "include/density.h"

namespace Smash {

bool particle_in_denstype(const PdgCode pdg, Density_type dens_type) {
  if ( (dens_type == baryon  && (pdg.baryon_number() != 0) ) ||
       (dens_type == proton  && pdg == 0x2212) ||
       (dens_type == neutron && pdg == 0x2112) ) {
    return true;
  } else {
    return false;
  }
}

FourVector four_current(const ThreeVector &r, const ParticleList &plist,
                        double gs_sigma, Density_type dens_type) {
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  double tmp;

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (6*gs_sigma) * (6*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    const double inv_gammai = p.inverse_gamma();

    // Get distance between particle and r in the particle rest frame
    tmp = ((r - ri) * betai) / (inv_gammai * (1. + inv_gammai));
    const ThreeVector dr_rest = r - ri + betai * tmp;

    tmp = std::exp(- 0.5 * dr_rest.sqr() / (gs_sigma * gs_sigma)) / inv_gammai;
    if (dens_type == baryon) {
      tmp *= p.pdgcode().baryon_number();
    }
    jmu += FourVector(1., betai) * tmp;
  }

  const double norm = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jmu /= norm;

  // j^0 = jmu.x0() is computational frame density
  // jmu.abs() = sqrt(j^mu j_mu) is Eckart rest frame density
  return jmu;
}

std::pair<double, ThreeVector> rho_eckart_gradient(const ThreeVector &r,
                                const ParticleList &plist, double gs_sigma,
                                Density_type dens_type) {
  // baryon four-current in computational frame
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  // derivatives of baryon four-current in computational frame
  FourVector djmu_dx(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dy(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dz(0.0, 0.0, 0.0, 0.0);
  double tmp1, tmp2;

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (6*gs_sigma) * (6*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    const double inv_gammai = p.inverse_gamma();

    // Get distance between particle and r in the particle rest frame
    tmp1 = inv_gammai * (1. + inv_gammai);
    const ThreeVector dr_rest = r - ri + betai * (((r - ri) * betai) / tmp1);

    tmp2 = std::exp(- 0.5 * dr_rest.sqr() / (gs_sigma*gs_sigma)) / inv_gammai;
    if (dens_type == baryon) {
      tmp2 *= p.pdgcode().baryon_number();
    }
    jmu += FourVector(1., betai) * tmp2;

    // Calculate the gradient: d(0.5 * r_rest^2)/d \vec{r}
    const ThreeVector drest2_dr = dr_rest + betai * ((dr_rest * betai) / tmp1);

    djmu_dx += FourVector(1., betai) * (tmp2 * drest2_dr.x1());
    djmu_dy += FourVector(1., betai) * (tmp2 * drest2_dr.x2());
    djmu_dz += FourVector(1., betai) * (tmp2 * drest2_dr.x3());
  }
  const double norm1 = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jmu /= norm1;
  const double norm2 = - norm1 * gs_sigma*gs_sigma;
  djmu_dx /= norm2;
  djmu_dy /= norm2;
  djmu_dz /= norm2;

  // Eckart rest frame density
  const double rho = jmu.abs();

  // Eckart rest frame density and its gradient
  return std::make_pair(rho, ThreeVector(jmu.Dot(djmu_dx),
                                         jmu.Dot(djmu_dy),
                                         jmu.Dot(djmu_dz)) / rho);
}

void vtk_density_map(const char * file_name, const ParticleList &plist,
                     double gs_sigma, Density_type dens_type,
                     int nx, int ny, int nz, double dx, double dy, double dz) {
  ThreeVector r;
  double rho_eck;
  std::ofstream a_file;
  a_file.open(file_name, std::ios::out);
  a_file << "# vtk DataFile Version 2.0\n" <<
            "density\n" <<
            "ASCII\n" <<
            "DATASET STRUCTURED_POINTS\n" <<
            "DIMENSIONS " << 2*nx+1 << " " << 2*ny+1 << " " << 2*nz+1 <<"\n" <<
            "SPACING 1 1 1\n"
            "ORIGIN " << -nx << " " << -ny << " " << -nz << "\n" <<
            "POINT_DATA " << (2*nx+1)*(2*ny+1)*(2*nz+1) << "\n" <<
            "SCALARS density float 1\n" <<
            "LOOKUP_TABLE default\n";

  a_file << std::setprecision(8);
  a_file << std::fixed;
  for (auto iz = -nz; iz <= nz; iz++) {
    for (auto iy = -ny; iy <= ny; iy++) {
      for (auto ix = -nx; ix <= nx; ix++) {
        r = ThreeVector(ix*dx, iy*dy, iz*dz);
        rho_eck = four_current(r, plist, gs_sigma, dens_type).abs();
        a_file << rho_eck << " ";
      }
      a_file << "\n";
    }
  }
  a_file.close();
}

void density_along_line(const char * file_name, const ParticleList &plist,
                        double gs_sigma, Density_type dens_type,
                        const ThreeVector &line_start,
                        const ThreeVector &line_end, int n_points) {
  ThreeVector r;
  double rho_eck;
  std::ofstream a_file;
  a_file.open(file_name, std::ios::out);

  for (int i = 0; i <= n_points; i++) {
    r = line_start + (line_end - line_start) * (1.0 * i / n_points);
    rho_eck = four_current(r, plist, gs_sigma, dens_type).abs();
    a_file << r.x1() << " " <<
              r.x2() << " " <<
              r.x3() << " " << rho_eck << "\n";
  }
  a_file.close();
}

}  // namespace Smash
