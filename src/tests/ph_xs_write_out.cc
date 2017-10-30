#include "unittest.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include "../include/kinematics.h"
#include "../include/photoncrosssections.h"
#include "../include/tabulationnd.h"
#include "setup.h"

using namespace Smash;

const std::string base_path =
    "/home/jonas/Master/cross_sections_tests/diff_xs/";

double random_double(double min, double max) {
  double random = ((double)rand()) / (double)RAND_MAX;
  double range = max - min;
  return (random * range) + min;
}

const double ds_val[] = {0.01, 0.05, 0.01};
const double dt_val[] = {0.005, 0.05, 0.01};

TEST(uniform_random_pi0_rho0_pi0) {
  // sample N points and interpolate value there
  const int N = 40000;

  double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02, dt = 0.01,
         ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  for (int i = 0; i < 3; i++) {
    dt = dt_val[i];
    ds_sq = ds_val[i];
    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                         xs_object.xs_diff_pi0_rho0_pi0);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0/"
       << "unirand_dt_" << dt << "_ds_" << ds_sq << ".dat";

    fs.precision(10);
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
    fs << "#s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference <\n";
    fs << "# " << tab2.get_size() << " points\n";
    // sample at n random points between tabulated values
    for (int i = 0; i < N; i++) {
      s_sqrt = random_double(s0_sqrt, s1_sqrt);
      s = s_sqrt * s_sqrt;
      t = random_double(t0, t1);

      xs = xs_object.xs_diff_pi0_rho0_pi0(s, t);
      xs_interp = tab2.get_linear(s, t);
      double divisor = (xs != 0) ? xs : 1;
      double rel_diff = std::abs((xs - xs_interp)) / divisor;
      fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
         << rel_diff << "\n";
    }

    fs.close();
  }
}
TEST(uniform_random_pi_pi_rho0) {
  // sample N points and interpolate value there
  const int N = 40000;
  double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02, dt = 0.01,
         ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  for (int i = 0; i < 3; i++) {
    dt = dt_val[i];
    ds_sq = ds_val[i];

    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                         xs_object.xs_diff_pi_pi_rho0);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi_pi_rho0/"
       << "unirand_dt_" << dt << "_ds_" << ds_sq << ".dat";

    fs.precision(10);
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
    fs << "#s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference \n";
    fs << "# " << tab2.get_size() << " points\n";
    // sample at n random points between tabulated values
    for (int i = 0; i < N; i++) {
      s_sqrt = random_double(s0_sqrt, s1_sqrt);
      s = s_sqrt * s_sqrt;
      t = random_double(t0, t1);

      xs = xs_object.xs_diff_pi_pi_rho0(s, t);
      xs_interp = tab2.get_linear(s, t);
      double divisor = (xs != 0) ? xs : 1;
      double rel_diff = std::abs((xs - xs_interp)) / divisor;
      fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
         << rel_diff << "\n";
    }

    fs.close();
  }
}
TEST(uniform_random_pi_pi0_rho) {
  // sample N points and interpolate value there
  const int N = 40000;
  double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02, dt = 0.01,
         ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  for (int i = 0; i < 3; i++) {
    dt = dt_val[i];
    ds_sq = ds_val[i];

    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                         xs_object.xs_diff_pi_pi0_rho);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi_pi0_rho/"
       << "unirand_dt_" << dt << "_ds_" << ds_sq << ".dat";

    fs.precision(10);
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
    fs << "# s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
    fs << "# " << tab2.get_size() << "\n";

    // sample at n random points between tabulated values
    for (int i = 0; i < N; i++) {
      s_sqrt = random_double(s0_sqrt, s1_sqrt);
      s = s_sqrt * s_sqrt;
      t = random_double(t0, t1);

      xs = xs_object.xs_diff_pi_pi0_rho(s, t);
      xs_interp = tab2.get_linear(s, t);
      double divisor = (xs != 0) ? xs : 1;
      double rel_diff = std::abs((xs - xs_interp)) / divisor;
      fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
         << rel_diff << "\n";
    }

    fs.close();
  }
}
TEST(uniform_random_pi0_rho_pi) {
  // sample N points and interpolate value there
  const int N = 40000;
  double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.20, t1 = -0.03, dt = 0.01,
         ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  for (int i = 0; i < 3; i++) {
    dt = dt_val[i];
    ds_sq = ds_val[i];

    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                         xs_object.xs_diff_pi0_rho_pi);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho_pi/"
       << "unirand_dt_" << dt << "_ds_" << ds_sq << ".dat";

    fs.precision(10);
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
    fs << "# s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
    fs << "# " << tab2.get_size() << "\n";
    // sample at n random points between tabulated values
    for (int i = 0; i < N; i++) {
      s_sqrt = random_double(s0_sqrt, s1_sqrt);
      s = s_sqrt * s_sqrt;
      t = random_double(t0, t1);

      xs = xs_object.xs_diff_pi0_rho_pi(s, t);
      xs_interp = tab2.get_linear(s, t);
      double divisor = (xs != 0) ? xs : 1;
      double rel_diff = std::abs((xs - xs_interp)) / divisor;
      fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
         << rel_diff << "\n";
    }

    fs.close();
  }
}


TEST(uniform_random_pi_rho_pi0) {
  // sample N points and interpolate value there
  const int N = 40000;
  double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02, dt = 0.01,
         ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  for (int i = 0; i < 3; i++) {
    dt = dt_val[i];
    ds_sq = ds_val[i];

    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                         xs_object.xs_diff_pi_rho_pi0);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi_rho_pi0/"
       << "unirand_dt_" << dt << "_ds_" << ds_sq << ".dat";

    fs.precision(10);
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
    fs << "# s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
    fs << "# " << tab2.get_size() << "\n";
    // sample at n random points between tabulated values
    for (int i = 0; i < N; i++) {
      s_sqrt = random_double(s0_sqrt, s1_sqrt);
      s = s_sqrt * s_sqrt;
      t = random_double(t0, t1);

      xs = xs_object.xs_diff_pi_rho_pi0(s, t);
      xs_interp = tab2.get_linear(s, t);
      double divisor = (xs != 0) ? xs : 1;
      double rel_diff = std::abs((xs - xs_interp)) / divisor;
      fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
         << rel_diff << "\n";
    }

    fs.close();
  }
}
TEST(pi_rho_pi0) {
    const double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02,
                 dt = 0.005, ds = 0.01;
    const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
    PhotonCrossSection<ComputationMethod::Analytic> xs_object;
    TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi_rho_pi0);
    std::fstream fs;
    std::stringstream ss;
    ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi_rho_pi0/"
       << "dt_" << dt << "_ds_" << ds << ".dat";
    fs.open(ss.str(), std::fstream::out);
    std::cout << ss.str();
    double t, xs, s, s_sqrt, xs_interp;
    fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
    for (s_sqrt = s0_sqrt; s_sqrt <= s1_sqrt; s_sqrt += 0.2 * ds) {
      for (t = t0; t < t1; t += 0.2 * dt) {
        s = s_sqrt * s_sqrt;
        xs = xs_object.xs_diff_pi_rho_pi0(s, t);
        xs_interp = tab2.get_linear(s, t);
        double divisor = (xs != 0) ? xs : 1;
        double rel_diff = std::abs((xs - xs_interp)) / divisor;
        fs << s_sqrt << " " << t << " " << xs << " " << xs_interp << " "
           << rel_diff << "\n";
      }
      fs << std::endl;
    }
    fs.close();
  }

