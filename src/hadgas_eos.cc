/*
 *
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <gsl/gsl_sf_bessel.h>

#include <iostream>

#include "include/hadgas_eos.h"

namespace Smash {

HadgasEos::HadgasEos() :
  x_(gsl_vector_alloc(n_equations_)) {
  const gsl_multiroot_fsolver_type *solver_type;
  solver_type = gsl_multiroot_fsolver_hybrids;
  solver_ = gsl_multiroot_fsolver_alloc(solver_type, n_equations_);
}

HadgasEos::~HadgasEos() {
  gsl_multiroot_fsolver_free(solver_);
  gsl_vector_free(x_);
}

double HadgasEos::hadgas_partial_density(const ParticleType& ptype,
                                         double beta, double mub, double mus) {
  const double z = ptype.mass()*beta;
  const double mu = ptype.baryon_number()*mub + ptype.strangeness()*mus;
  const unsigned int g = ptype.spin() + 1;
  return z*z * g * std::exp(mu*beta) * gsl_sf_bessel_Kn(2, z);
}

double HadgasEos::hadgas_energy_density(double T, double mub, double mus) {
  const double beta = 1.0/T;
  double e = 0.0;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_hadron()) {
      continue;
    }
    const double z = ptype.mass()*beta;
    const double mu = mub*ptype.baryon_number() + mus*ptype.strangeness();
    const unsigned int g = ptype.spin() + 1;
    e += z*z * g * std::exp(mu*beta) * (3.0*gsl_sf_bessel_Kn(2, z) +
                                       gsl_sf_bessel_K1(z));
  }
  e *= prefactor_ * T*T*T*T;
  return e;
}

double HadgasEos::hadgas_density(double T, double mub, double mus) {
  const double beta = 1.0/T;
  double rho = 0.0;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_hadron()) {
      continue;
    }
    rho += hadgas_partial_density(ptype, beta, mub, mus);
  }
  rho *= prefactor_ * T*T*T;
  return rho;
}

double HadgasEos::hadgas_net_baryon_density(double T, double mub, double mus) {
  const double beta = 1.0/T;
  double rho = 0.0;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_baryon()) {
      continue;
    }
    rho += hadgas_partial_density(ptype, beta, mub, mus) *
           ptype.baryon_number();
  }
  rho *= prefactor_ * T*T*T;
  return rho;
}

double HadgasEos::hadgas_net_strange_density(double T, double mub, double mus) {
  const double beta = 1.0/T;
  double rho = 0.0;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (ptype.strangeness() == 0) {
      continue;
    }
    rho += hadgas_partial_density(ptype, beta, mub, mus) *
           ptype.strangeness();
  }
  rho *= prefactor_ * T*T*T;
  return rho;
}

double HadgasEos::mus_net_strangeness0(double T, double mub) {
  // Binary search
  double mus_u = mub + T;
  double mus_l = 0.0;
  double mus, rhos;
  int iteration = 0;
  // 30 iterations should give precision 2^-30 ~ 10^-9
  const int max_iteration = 30;
  do {
    mus = 0.5 * (mus_u + mus_l);
    rhos = hadgas_net_strange_density(T, mub, mus);
    if (rhos > 0.0) {
      mus_u = mus;
    } else {
      mus_l = mus;
    }
    iteration++;
  } while (std::abs(rhos) > tolerance_ && iteration < max_iteration);
  if (iteration == max_iteration) {
    throw std::runtime_error("Solving rho_s = 0: too many iterations.");
  }
  return mus;
}

int HadgasEos::hadgas_eos_equations(const gsl_vector* x,
                                    void *params, gsl_vector* f) {
  double e =  reinterpret_cast<struct rparams*>(params)->e;
  double nb = reinterpret_cast<struct rparams*>(params)->nb;
  double ns = reinterpret_cast<struct rparams*>(params)->ns;

  const double T   = gsl_vector_get(x, 0);
  const double mub = gsl_vector_get(x, 1);
  const double mus = gsl_vector_get(x, 2);

  gsl_vector_set(f, 0, hadgas_energy_density(T, mub, mus) - e);
  gsl_vector_set(f, 1, hadgas_net_baryon_density(T, mub, mus) - nb);
  gsl_vector_set(f, 2, hadgas_net_strange_density(T, mub, mus) - ns);

  return GSL_SUCCESS;
}

std::array<double, 3> HadgasEos::solve_hadgas_eos(double e,
                                                  double nb, double ns) {
  int status;
  size_t iter = 0;

  struct rparams p = {e, nb, ns};
  gsl_multiroot_function f = {&HadgasEos::hadgas_eos_equations,
                              n_equations_, &p};
  // Initial approximation
  gsl_vector_set(x_, 0, 0.15);
  gsl_vector_set(x_, 1, 0.2);
  gsl_vector_set(x_, 2, 0.05);

  gsl_multiroot_fsolver_set(solver_, &f, x_);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(solver_);

    // print_solver_state(iter);

    // check if solver is stuck
    if (status) {
        break;
    }
    status = gsl_multiroot_test_residual(solver_->f, tolerance_);
  } while (status == GSL_CONTINUE && iter < 1000);

  if (status != GSL_SUCCESS) {
    throw std::runtime_error(gsl_strerror(status));
  }
  return {gsl_vector_get(solver_->x, 0),
          gsl_vector_get(solver_->x, 1),
          gsl_vector_get(solver_->x, 2)};
}

void HadgasEos::print_solver_state(size_t iter) {
  std::cout <<
          "iter = " << iter << "," <<
          " x = "   << gsl_vector_get(solver_->x, 0) << " " <<
                       gsl_vector_get(solver_->x, 1) << " " <<
                       gsl_vector_get(solver_->x, 2) << ", " <<
          "f(x) = " << gsl_vector_get(solver_->f, 0) << " " <<
                       gsl_vector_get(solver_->f, 1) << " " <<
                       gsl_vector_get(solver_->f, 2) << std::endl;
}


}  // namespace Smash
