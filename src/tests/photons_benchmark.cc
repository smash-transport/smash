#include "unittest.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include "../include/kinematics.h"
#include "../include/particletype.h"
#include "../include/photoncrosssections.h"
#include "../include/tabulationnd.h"
#include "setup.h"
#include <vector>

#include <chrono>

using namespace Smash;

// idea: draw N random numbers in some range, call analytic and lookup, compare timing


std::vector<std::string> Process {
  "pi0_rho0_pi0",
  "pi0_rho_pi",
  "pi_pi0_rho",
  "pi_pi_rho0",
  "pi_rho0_pi",
  "pi_rho_pi0"
}; 

double random_double(double min, double max) {
  double random = ((double)rand()) / (double)RAND_MAX;
  double range = max - min;
  return (random * range) + min;
}

TEST(time_pi0_rho0) {

  int N = 10e7;
  std::vector<size_t> N_vec {1000, 10000, 100000, 1000000, 10000000, 100000000};

  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::fstream f_;
  f_.open("results.dat", std::fstream::out);

  for (size_t N = 128; N < 1e8; N *= 2) {
  // init arrays of randomly drawn s,t. m_rho is fixed
  std::vector<double> s_vec(N);
  std::vector<double> t1_vec(N);
  std::vector<double> t2_vec(N);
  double m_rho = 0.776, m_pi = 0.139;
  double s0 = std::sqrt(m_rho);
  double s1 = 20.0;
  for (int i = 0; i < N; i++) {
    double s = random_double(s0, s1);
    double sqrts = std::sqrt(s);
    std::array<double, 2> mandelstam_t = get_t_range(sqrts, m_pi, m_rho, m_pi, 0.0);
    const double t1 = mandelstam_t[1];
    const double t2 = mandelstam_t[0];
    
    s_vec[i] = s;
    t1_vec[i] = t1;
    t2_vec[i] = t2;

  }

  std::vector<double> xs_vec_tab(N), xs_vec_an(N);
  double xs = xs_tab.xs_diff_pi0_rho0_pi0(s_vec[0], t1_vec[0], m_rho);

  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < N; i++) {
    xs_vec_tab[i] = xs_tab.xs_diff_pi0_rho0_pi0(s_vec[i], t1_vec[i], m_rho);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto t_tab = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();

  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < N; i++) {
    xs_vec_an[i] = xs_an.xs_diff_pi0_rho0_pi0(s_vec[i], t1_vec[i], m_rho);
  }  
  t2 = std::chrono::high_resolution_clock::now();
  auto t_an = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();

  f_ << N << " " << t_tab << " " << t_an << std::endl;
}
  f_.close();
}