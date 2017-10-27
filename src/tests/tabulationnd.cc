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

double f_sin(double x, double y) { return x * x + y * y; }

/*


TEST(max_error) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.010;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  double ds = 0.1, dt = 0.5 * ds;

  double s, t, divisor, rel_diff, max_diff = 0.0, xs, xs_interp, s_sqrt;
  bool finished = false;
  auto max_coords = std::make_tuple(s, t, max_diff);

  while (!finished) {
    TabulationND<2> tab(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi0_rho0_pi0);
    for (s_sqrt = s0_sqrt; s_sqrt <= s1_sqrt; s_sqrt += ds * 0.2) {
      for (t = t0; t < t1; t += dt * 0.2) {
        s = s_sqrt * s_sqrt;
        xs = xs_object.xs_diff_pi0_rho0_pi0(s, t);
        xs_interp = tab.get_linear(s, t);
        divisor = (xs != 0) ? xs : 1;
        rel_diff = std::abs((xs - xs_interp)) / divisor;
        max_diff = (rel_diff > max_diff) ? rel_diff : max_diff;
        max_coords = std::make_tuple(s, t, max_diff);
      }
    }

    if (std::get<2>(max_coords) < 0.02) {
      finished = true;
      break;
    } else {
      std::cout
          << "-----------------------------------------------------------\n";
      std::cout << "ds= " << ds << ", dt= " << dt
                << " max.error= " << std::get<2>(max_coords) << " at "
                << std::get<0>(max_coords) << " " << std::get<1>(max_coords) <<
std::endl; dt *= 0.75; ds *= 0.75;
    }
  }
  std::cout << "ds, dt " << ds << " " << dt << "\n";
  std::cout << "Max. error " << max_diff << std::endl;
}
*/
TEST(rand_pi_rho_pi0) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.1, t0 = -0.22, t1 = -0.02, dt = 0.01,
               ds_sq = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  TabulationND<2> tab2(s0, s1, t0, t1, ds_sq * ds_sq, dt,
                       xs_object.xs_diff_pi_rho_pi0);
  std::fstream fs;
  std::stringstream ss;
  ss << "/home/jonas/Master/cross_sections_tests/diff_xs/pi_rho_pi0/"
     << "rand_dt_" << dt << "_ds_" << ds_sq << ".dat";

  fs.precision(10);
  fs.open(ss.str(), std::fstream::out);
  std::cout << ss.str();
  double t, xs, s, s_sqrt, xs_interp, s_s, t_s;
  fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
  // sample at n random points between tabulated values
  const int n_sample = 10;
  for (s_sqrt = s0_sqrt; s_sqrt <= s1_sqrt - 2 * ds_sq; s_sqrt += ds_sq) {
    for (t = t0; t < t1; t += dt) {
      for (int i = 0; i < n_sample; i++) {
        s = s_sqrt * s_sqrt;

        double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        double r2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        assert(r < 1 && r2 < 1);
        double s_sqrt_s = s_sqrt + r * ds_sq;
        s_s = s_sqrt_s * s_sqrt_s;
        t_s = t + r2 * dt;
        xs = xs_object.xs_diff_pi_rho_pi0(s_s, t_s);
        xs_interp = tab2.get_linear(s_s, t_s);
        double divisor = (xs != 0) ? xs : 1;
        double rel_diff = std::abs((xs - xs_interp)) / divisor;
        fs << s_sqrt_s << " " << t_s << " " << xs << " " << xs_interp << " "
           << rel_diff << "\n";
      }
    }
    fs << std::endl;
  }
  fs.close();
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
  // sample at ten points between tabulated values
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

/*
TEST(pi0_rho_pi) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.5, t0 = -0.22, t1 = -0.02,
               dt = 0.005, ds = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi0_rho_pi);
  std::fstream fs;

  fs.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho_pi/"
      "dt_005_ds_01.dat",
      std::fstream::out);

  double t, xs, s, s_sqrt, xs_interp;
  fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
  for (s_sqrt = s0_sqrt + 0.5 * ds; s_sqrt <= s1_sqrt; s_sqrt += ds) {
    for (t = t0 + 0.5 * dt; t < t1; t += dt) {
      s = s_sqrt * s_sqrt;
      xs = xs_object.xs_diff_pi0_rho_pi(s, t);
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

*/
/*
TEST(sin) {

  const double x0 = -1., x1 = +1., y0 = -1., y1 = 1., dx = 0.01, dy = 0.01;
  TabulationND<2> tab(x0, x1, y0, y1, dx, dy, f_sin);
  std::fstream fs;
  fs.open("/home/jonas/Master/cross_sections_tests/diff_xs/easy_functions/xxyy.dat",
std::fstream::out); double x, y; double x_poll = 0.1 * dx, y_poll = 0.1 * dy;
  for (x = x0; x < x1; x += x_poll)
  {
    for (y = y0; y < y1; y += y_poll)
    {
      double z_tab = tab.get_linear(x,y);
      double z_an = f_sin(x,y);
      double divisor = (z_an != 0) ? z_an : 1;
      double rel_diff = (std::abs(z_tab - z_an) / divisor);
      fs << x << " " << y << " " << z_an << " " << z_tab << " " << rel_diff <<
"\n";

    }
    fs << std::endl;
  }
  fs.close();


}

TEST(pi0_rho0_pi0) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.5, t0 = -0.22, t1 = -0.04,
               dt = 0.01, ds = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi0_rho0_pi0);
  std::fstream fs;

  fs.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0/"
      "dt_01_ds_01.dat",
      std::fstream::out);

  double t, xs, s, s_sqrt, xs_interp;
  fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
  for (s_sqrt = s0_sqrt; s_sqrt <= s1_sqrt; s_sqrt += 0.1*ds) {
    for (t = t0; t < t1; t += 0.1*dt) {
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


TEST(pi_pi0_rho) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.5, t0 = -0.22, t1 = -0.02,
               dt = 0.005, ds = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi_pi0_rho);
  std::fstream fs;

  fs.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi_pi0_rho/"
      "dt_005_ds_01.dat",
      std::fstream::out);

  double t, xs, s, s_sqrt, xs_interp;
  fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
  for (s_sqrt = s0_sqrt + 0.5 * ds; s_sqrt <= s1_sqrt; s_sqrt += ds) {
    for (t = t0 + 0.5 * dt; t < t1; t += dt) {
      s = s_sqrt * s_sqrt;
      xs = xs_object.xs_diff_pi_pi0_rho(s, t);
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

TEST(pi_pi_rho0) {
  const double s0_sqrt = 0.95, s1_sqrt = 1.5, t0 = -0.22, t1 = -0.02,
               dt = 0.005, ds = 0.01;
  const double s0 = s0_sqrt * s0_sqrt, s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, xs_object.xs_diff_pi_pi_rho0);
  std::fstream fs;

  fs.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi_pi_rho0/"
      "dt_005_ds_01.dat",
      std::fstream::out);

  double t, xs, s, s_sqrt, xs_interp;
  fs << "s_sqrt \t t \t xs_analytic \t xs_interp \t rel_difference\n";
  for (s_sqrt = s0_sqrt + 0.5 * ds; s_sqrt <= s1_sqrt; s_sqrt += ds) {
    for (t = t0 + 0.5 * dt; t < t1; t += dt) {
      s = s_sqrt * s_sqrt;
      xs = xs_object.xs_diff_pi_pi_rho0(s, t);
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
*/
/*
TEST(xs_write_1d) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  // later in loop

  // start at mass of rho
  const double s0_sqrt = 0.95;
  const double s1_sqrt = 3.2;
  const double s0 = s0_sqrt * s0_sqrt;
  const double s1 = s1_sqrt * s1_sqrt;

  const double ds = 0.04;
  double xs;

  TabulationND<1> tab(s0, s1, ds, xs_object.xs_pi0_rho0_pi0);
  std::cout << tab.n_points() << "\n";
  std::vector<double> xs_an;
  std::vector<double> xs_tab, xvalues_poll;
  double ds_poll = 0.0001;
  for (double s_sqrt = s0_sqrt; s_sqrt <= s1_sqrt; s_sqrt += ds_poll) {
    xs = xs_object.xs_pi0_rho0_pi0(s_sqrt * s_sqrt);
    // xs = 0;
    xs_an.push_back(xs);
    xs = tab.get_linear(s_sqrt * s_sqrt);

    xs_tab.push_back(xs);
    xvalues_poll.push_back(s_sqrt);
  }

  std::fstream fs;
  fs.open("/home/jonas/Master/cross_sections_tests/pi0_rho0_pi0.dat",
          std::fstream::out);
  for (int i = 0; i < xs_tab.size(); i++) {
    double diff = (xs_an[i] - xs_tab[i]) / xs_an[i];
    fs << xvalues_poll[i] << " " << xs_an[i] << " " << xs_tab[i] << " " <<
diff
       << "\n";
  }
  fs.close();
}

TEST(xs_write_2d) {
  // double s0 = 0., s1 = 10.0, t0 = 0., t1 = 10., ds = 0.1, dt = 0.1;
  // COMPARE(tab2.get_linear(0.0, 0.0), 2.);
  double s0_sqrt = 0.99, s1_sqrt = 1.1, t0 = -1.8, t1 = -0.8, ds = 0.01,
         dt = 0.01;
  double s0 = s0_sqrt * s0_sqrt;
  double s1 = s1_sqrt * s1_sqrt;
  PhotonCrossSection<ComputationMethod::Analytic> xs_diff;
  // TabulationND<2> tab2(s0, s1, t0, t1, ds, dt, [](double x, double
y){return
  // x+y;});
  TabulationND<2> tab2(s0, s1, t0, t1, ds, dt,
xs_diff.xs_diff_pi0_rho0_pi0);

  double s, t, xs;
  std::fstream fs, fs2, fs3, fs4;
  fs.open("/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0_an.dat",
          std::fstream::out);
  fs2.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0_int.dat",
      std::fstream::out);
  fs3.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0_tab.dat",
      std::fstream::out);
  fs4.open(
      "/home/jonas/Master/cross_sections_tests/diff_xs/pi0_rho0_pi0_diff.dat",
      std::fstream::out);

  double s_sqrt;
  for (s_sqrt = s0_sqrt + 0.5 * ds; s_sqrt < s1_sqrt; s_sqrt += ds) {
    for (t = t0 + 0.5 * dt; t < t1; t += dt) {
      s = s_sqrt * s_sqrt;
      double xs_an = xs_diff.xs_diff_pi0_rho0_pi0(s, t);
      fs << s_sqrt << " " << t << " " << xs_an << "\n";
      double xs_tab = tab2.get_linear(s, t);
      fs2 << s_sqrt << " " << t << " " << xs_tab << "\n";
      // fs3 << s << " " << t << " " << tab2.get_closest(s,t) << "\n";
      double divisor = (xs_an != 0) ? xs_an : 1.;
      fs4 << s_sqrt << " " << t << " " << (xs_an - xs_tab) / divisor <<
"\n";
    }
    fs << std::endl;
    fs2 << std::endl;
    fs3 << std::endl;
    fs4 << std::endl;
  }
  fs.close();
  fs2.close();
  fs3.close();
  fs4.close();
}

/*
TEST(constant) {
  const TabulationND<1> tab(0.0, 10.0, 0.1, [](double x) { return 1.; });
  COMPARE(tab.n_points(), 101);
  FUZZY_COMPARE(tab.get_linear(0.0), 1.);
  FUZZY_COMPARE(tab.get_linear(0.4), 1.);
  FUZZY_COMPARE(tab.get_linear(0.15), 1.);
  FUZZY_COMPARE(tab.get_linear(10.0), 1.);
  FUZZY_COMPARE(tab.get_linear(5.0), 1.);
}

const double& mpion = m_pion_, &mrho = m_rho_;
auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho_, m_pion_, 0.);
const double &t1 = t_mandelstam[1];
const double &t2 = t_mandelstam[0];TEST(linear) {
  const TabulationND<1> tab(0.0, 10.0, 0.1, [](double x) { return 2 * x; });
  COMPARE(tab.n_points(), 101);

  for (int i = 0; i < tab.n_points(); i++)
  {
    std::cout << tab.get_from_index(i) << "\n";
  }


  COMPARE_RELATIVE_ERROR(tab.get_linear(0.0), 0.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(1.0), 2.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(10.0), 20.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(5.5), 11.0, 0.05);
}
*/
