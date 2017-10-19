#include "unittest.h"

#include "../include/tabulationnd.h"
#include <fstream>
#include <iostream>
#include <string>
#include "../include/kinematics.h"
#include "../include/photoncrosssections.h"
#include "setup.h"


using namespace Smash;

TEST(xs_write) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_object;
  // later in loop

  // start at mass of rho
  const double s0_sqrt = 0.95;
  const double s1_sqrt = 3.2;
  const double s0 = s0_sqrt * s0_sqrt;
  const double s1 = s1_sqrt * s1_sqrt;

  const double ds = 0.01;
  double xs;

  TabulationND<1> tab(s0, s1, ds, xs_object.xs_pi0_rho0_pi0);

  std::vector<double> xs_an;
  std::vector<double> xs_tab, xvalues_poll;
  double ds_poll = 0.0001;
  for (double s_sqrt = s0_sqrt; s_sqrt <=s1_sqrt; s_sqrt += ds_poll) {
    xs = xs_object.xs_pi0_rho0_pi0(s_sqrt * s_sqrt);
    //xs = 0;
    xs_an.push_back(xs);
    xs = tab.get_linear(s_sqrt * s_sqrt);

    xs_tab.push_back(xs);
    xvalues_poll.push_back(s_sqrt);
  }

  std::fstream fs;
  fs.open("/home/jonas/Master/cross_sections_tests/pi0_rho0_pi0.dat",
          std::fstream::out);
  for (int i = 0; i < xs_tab.size(); i++)
  {
    double diff = (xs_an[i] - xs_tab[i]) / xs_an[i];
    fs << xvalues_poll[i] << " " << xs_an[i] << " " << xs_tab[i] << " " << diff << "\n";
  }
  fs.close();
}


TEST(constant) {
  const TabulationND<1> tab(0.0, 10.0, 0.1, [](double x) { return 1.; });
  COMPARE(tab.n_points(), 101);
  FUZZY_COMPARE(tab.get_linear(0.0), 1.);
  FUZZY_COMPARE(tab.get_linear(0.4), 1.);
  FUZZY_COMPARE(tab.get_linear(0.15), 1.);
  FUZZY_COMPARE(tab.get_linear(10.0), 1.);
  FUZZY_COMPARE(tab.get_linear(5.0), 1.);
}

TEST(linear) {
  const TabulationND<1> tab(0.0, 10.0, 0.1, [](double x) { return 2 * x; });
  COMPARE(tab.n_points(), 101);
  /*
  for (int i = 0; i < tab.n_points(); i++)
  {
    std::cout << tab.get_from_index(i) << "\n";
  }
  */

  COMPARE_RELATIVE_ERROR(tab.get_linear(0.0), 0.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(1.0), 2.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(10.0), 20.0, 0.05);
  COMPARE_RELATIVE_ERROR(tab.get_linear(5.5), 11.0, 0.05);
}

TEST(write_out) {
  const std::string path{"/home/jonas/Master/cross_sections_tests/file.dat"};
  std::fstream fs;
  fs.open(path, std::fstream::out);
  fs << "hi"
     << "\n";
  fs.close();
}
