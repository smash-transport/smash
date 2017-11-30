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

using namespace Smash;

const std::string basepath = "/home/jonas/Master/cross_sections_tests/stable/";

std::vector<std::string> Process {
  "pi0_rho0_pi0",
  "pi0_rho_pi",
  "pi_pi0_rho",
  "pi_pi_rho0",
  "pi_rho0_pi",
  "pi_rho_pi0"
}; 


void diff_stable(std::string proc)
{
  const double s0 = 0.1, s1 = 5.0, ds = 0.01;
  const double t0 = -5.0, t1 = 5.0, dt = 0.01;
  const double mrho = 0.776;
  std::cout << Process[0];
  std::stringstream ss;
  ss << basepath << "diff/ " << proc;
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;
  
  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  {
    // can not swtich on string...
  if (proc == "pi0_rho0_pi0")
  {
    xsAn = xs_an.xs_diff_pi0_rho0_pi0(s, t, mrho);
    xsTab = xs_tab.xs_diff_pi0_rho0_pi0(s, t, mrho);
  }
  else if (proc == "pi0_rho_pi")
  {
    xsAn = xs_an.xs_diff_pi0_rho_pi(s, t, mrho);
    xsTab = xs_tab.xs_diff_pi0_rho_pi(s, t, mrho);
    
  }
  else if (proc == "pi_pi0_rho")
  {
    xsAn = xs_an.xs_diff_pi_pi0_rho(s,t,mrho);
    xsTab = xs_tab.xs_diff_pi_pi0_rho(s,t,mrho);
    
  }
  else if (proc == "pi_pi_rho0")
  {
    xsAn = xs_an.xs_diff_pi_pi_rho0(s, t, mrho);
    xsTab = xs_tab.xs_diff_pi_pi_rho0(s, t, mrho);
  }
  else if (proc == "pi_rho0_pi")
  {
    xsAn = xs_an.xs_diff_pi_rho0_pi(s, t, mrho);
    xsTab = xs_tab.xs_diff_pi_rho0_pi(s, t, mrho);
    
  }
  else if (proc == "pi_rho_pi0")
  {
    xsAn = xs_an.xs_diff_pi_rho_pi0(s, t, mrho);
    xsTab = xs_tab.xs_diff_pi_rho_pi0(s, t, mrho);
  }

  fs << s << " " << t << " " << xsAn << " " << xsTab << "\n";


  }

  fs.close();
}


double random_double(double min, double max) {
  double random = ((double)rand()) / (double)RAND_MAX;
  double range = max - min;
  return (random * range) + min;
}

const double ds_val[] = {0.01, 0.05, 0.01};
const double dt_val[] = {0.005, 0.05, 0.01};

TEST(diff_stable__)
{
  for (const auto &s : Process) {
    diff_stable(s);
  }
}

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "π               0.138   7.7e-9      111     211\n"
      "σ               0.800   0.400   9000221\n"
      "ρ               0.776   0.149       113     213\n"
      "ω               0.783   8.49e-3     223\n"
      "N       0.938 0         2112    2212\n"
      "Δ       1.232 0.117    1114    2114    2214    2224\n"
      "Λ        1.116 0         3122\n"
      "Λ(1520)  1.520 0.0156    3124\n"
      "Λ(1690)  1.690 0.0600   13124\n"
      "Σ       1.189 0        3112    3212    3222\n"
      "e⁻ 0.000511 0 11\n");
}

TEST(stable_pi0_rho0_pi0) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi0_rho0_pi0_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi0_rho0_pi0(s, m_rho);
    const double xsT = xs_tab.xs_pi0_rho0_pi0(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}

TEST(stable_pi0_rho_pi) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi0_rho_pi_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi0_rho_pi(s, m_rho);
    const double xsT = xs_tab.xs_pi0_rho_pi(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}
TEST(stable_pi_pi0_rho) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi_pi0_rho_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi_pi0_rho(s, m_rho);
    const double xsT = xs_tab.xs_pi_pi0_rho(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}

TEST(stable_pi_pi_rho0) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi_pi_rho0_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi_pi_rho0(s, m_rho);
    const double xsT = xs_tab.xs_pi_pi_rho0(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}

TEST(stable_pi_rho0_pi) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi_rho0_pi_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi_rho0_pi(s, m_rho);
    const double xsT = xs_tab.xs_pi_rho0_pi(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}

TEST(stable_pi_rho_pi0) {
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  std::stringstream ss;
  ss << basepath << "pi_rho_pi0_tot.dat";
  std::fstream fs;
  fs.open(ss.str(), std::fstream::out);

  double m_rho = 0.776;
  double s0 = 0.1, s1 = 5.0, ds = 0.01;
  double t0 = -5.0, t1 = 5.0, dt = 0.01;

  /***************
   * Total XS
   **************
   */
  for (double s = s0; s < s1; s += ds) {
    const double xsA = xs_an.xs_pi_rho_pi0(s, m_rho);
    const double xsT = xs_tab.xs_pi_rho_pi0(s, m_rho);
    fs << s << " " << xsA << " " << xsT << "\n";
  }
  fs.close();
}


/**********************************
 * BROAD RHO DIFFERENTIAL
 * ********************************
 */

TEST(pi0_rho0_pi0_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi0_rho0_pi0.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi0_rho0_pi0(s, t, m);
xsTab = xs_tab.xs_diff_pi0_rho0_pi0(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
TEST(pi0_rho_pi_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi0_rho_pi.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi0_rho_pi(s, t, m);
xsTab = xs_tab.xs_diff_pi0_rho_pi(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
TEST(pi_pi0_rho_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi_pi0_rho.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi_pi0_rho(s, t, m);
xsTab = xs_tab.xs_diff_pi_pi0_rho(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
TEST(pi_pi_rho0_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi_pi_rho0.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi_pi_rho0(s, t, m);
xsTab = xs_tab.xs_diff_pi_pi_rho0(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
TEST(pi_rho0_pi_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi_rho0_pi.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi_rho0_pi(s, t, m);
xsTab = xs_tab.xs_diff_pi_rho0_pi(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
TEST(pi_rho_pi0_broad_rho_diff) { 
 const double s0 = 0.1, s1 = 5.0, ds = 0.1;
  const double t0 = -5.0, t1 = 5.0, dt = 0.1;
  const double m0 = 0.1, m1 = 1.0, dm = 0.1;
  std::cout << Process[0];
  std::stringstream ss;
  std::fstream fs; ss << "/home/jonas/Master/cross_sections_tests/broad/diff/" << "pi_rho_pi0.dat";
  fs.open(ss.str(), std::fstream::out);
  double xsAn, xsTab;
  PhotonCrossSection<ComputationMethod::Analytic> xs_an;
  PhotonCrossSection<ComputationMethod::Lookup> xs_tab;

  for (double s = s0; s < s1; s += ds)
  for (double t = t0; t < t1; t += dt)
  for (double m = m0; m < m1; m += dm)
  { 
xsAn = xs_an.xs_diff_pi_rho_pi0(s, t, m);
xsTab = xs_tab.xs_diff_pi_rho_pi0(s,t,m);

    fs << s << " " << t << " " << m << " " << xsAn << " " << xsTab << "\n";
    }
  fs.close();
  }
    
