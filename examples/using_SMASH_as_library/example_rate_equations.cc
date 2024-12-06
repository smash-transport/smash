#include <filesystem>

#include "gsl/gsl_sf_bessel.h"

#include "smash/decaymodes.h"
#include "smash/experiment.h"
#include "smash/forwarddeclarations.h"
#include "smash/hadgas_eos.h"
#include "smash/integrate.h"
#include "smash/isoparticletype.h"
#include "smash/kinematics.h"
#include "smash/particledata.h"
#include "smash/scatteraction.h"
#include "smash/sha256.h"

using namespace smash;

static double thermal_average_sigmavrel(const ParticleTypePtr A_type,
                                        const ParticleTypePtr B_type,
                                        double T) {
  static smash::Integrator integrate;

  ParticleData A(*A_type), B(*B_type);
  const double ma = A.pole_mass(), mb = B.pole_mass();
  double time = 0.0, string_formation_time = 0.0;
  bool isotropic = true;
  const double norm =
      1.0 / (4.0 * ma * ma * mb * mb * T * gsl_sf_bessel_Kn(2, ma / T) *
             gsl_sf_bessel_Kn(2, mb / T));
  /* The keys under General are not used in this example but some of them are
   * required to create the experiment parameters. For demonstrative purposes we
   * provide here a pretty complete configuration, which is not strictly
   * speaking needed to create the ScatterActionsFinderParameters object. */
  Configuration exp_config{R"(
  General:
    Modus: Box
    End_Time: 100
    Nevents: 1
    Randomseed: -1
  Collision_Term:
    Collision_Criterion: Stochastic
    Strings: false
    Elastic_NN_Cutoff_Sqrts: 0
    Included_2to2: []
    Multi_Particle_Reactions: [Deuteron_3to2, A3_Nuclei_4to2]
    NNbar_Treatment: "no annihilation"
  Modi:
    Box:
      Length: 10.0
      Temperature: 0.2
      Start_Time: 0.0
      Use_Thermal_Multiplicities: True
      Initial_Condition: "thermal momenta"
  )",
                           Configuration::InitializeFromYAMLString};
  Configuration finder_config{R"(
  Collision_Term:
    Only_Warn_For_High_Probability: true
  )",
                              Configuration::InitializeFromYAMLString};
  const ScatterActionsFinderParameters finder_params(
      finder_config, smash::create_experiment_parameters(exp_config));
  // Clear configurations since we did not consume them all (the Configuration
  // destructor will give an error if called on a non-empty instance).
  exp_config.clear();
  finder_config.clear();
  const auto integral = integrate(0.0, 1.0, [&](double x) {
    const double m = (ma + mb) / x;
    const double jac = m / x;
    const double pcm = pCM(m, ma, mb);
    A.set_4momentum(ma, pcm, 0.0, 0.0);
    B.set_4momentum(mb, -pcm, 0.0, 0.0);
    ScatterActionPtr act = std::make_unique<ScatterAction>(
        A, B, time, isotropic, string_formation_time);
    act->add_all_scatterings(finder_params);
    const double sigma = act->cross_section();
    if (std::isnan(sigma) || sigma <= 0.0)
      return 0.0;
    const double g =
        (m * m - (ma + mb) * (ma + mb)) * (m * m - (ma - mb) * (ma - mb));
    return sigma * g * gsl_sf_bessel_Kn(1, m / T) * jac;
  });

  return integral * norm;
}

int main() {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π    0.138   7.7e-9  -  111         211\n"
      "N    0.938   0       +  2112        2212\n"
      "Λ    1.116   0       +  3122\n"
      "d    1.8756  0       +  1000010020\n"
      "t    2.8089  0       +  1000010030\n"
      "he3  2.8084  0       +  1000020030\n"
      "H3L  2.9934  0       +  1010010030\n");
  DecayModes::load_decaymodes("");
  sha256::Context hash_context;
  const auto hash = hash_context.finalize();
  IsoParticleType::tabulate_integrals(hash, std::filesystem::path("."));

  // 1) Prepare particle types
  const ParticleTypePtr pip = ParticleType::try_find(pdg::pi_p);
  const ParticleTypePtr pim = ParticleType::try_find(pdg::pi_m);
  const ParticleTypePtr pi0 = ParticleType::try_find(pdg::pi_z);
  const ParticleTypePtr p = ParticleType::try_find(pdg::p);
  const ParticleTypePtr n = ParticleType::try_find(pdg::n);
  const ParticleTypePtr La = ParticleType::try_find(pdg::Lambda);
  const ParticleTypePtr d = ParticleType::try_find(pdg::deuteron);
  const ParticleTypePtr t = ParticleType::try_find(pdg::triton);
  const ParticleTypePtr he3 = ParticleType::try_find(pdg::he3);
  const ParticleTypePtr H3L = ParticleType::try_find(pdg::hypertriton);
  // Make sure that necessary particles are in the table
  if (!(pip && pim && pi0)) {
    std::cout << "Pion missing in the particle table" << std::endl;
  }
  if (!(p && n && La)) {
    std::cout << "Nucleon or Lambda missing in the particle table" << std::endl;
  }
  if (!d) {
    std::cout << "Deuteron missing in the particle table" << std::endl;
  }
  if (!(t && he3 && H3L)) {
    std::cout << "Triton, He-3, or hypertriton H3L missing"
              << " in the particle table" << std::endl;
  }
  assert(pip && pim && pi0 && p && n && La && d && t && he3 && H3L);

  // 2) Prepare initial condition
  // Assuming that T(time) = const.
  // There is in fact no gaurantee that this must be true in a box simulation.
  // But explicit checks show that it tends to be true.
  const double T = 0.155;       // GeV, temperature
  const double V = 1000.0;      // fm^3, volume
  const double Npi_init = 180;  // This includes pi+, pi-, pi0
  const double Np_init = 60;
  const double Nn_init = 60;
  const double NLa_init = 0;
  const double Nd_init = 0;
  const double Nt_init = 0;
  const double Nhe3_init = 0;
  const double NH3L_init = 0;

  const double dt = 0.01;    // fm, integration timestep
  const double tend = 20.0;  // fm, end time

  std::cout << "\nRate Equation Example\n---------------------" << std::endl;

  // 3) Prepare averages cross sections
  // Here and further assuming that cross sections do not depend on pion charge
  const double sig_pid = thermal_average_sigmavrel(pip, d, T);
  const double sig_nd = thermal_average_sigmavrel(n, d, T);
  const double sig_pd = thermal_average_sigmavrel(p, d, T);
  assert(std::abs(sig_nd - sig_pd) < 1e-9);
  std::cout << "<sigma vrel> [mb] for pi+d, n+d, p+d: " << sig_pid << " "
            << sig_nd << " " << sig_pd << std::endl;

  const double sig_pit = thermal_average_sigmavrel(pip, t, T);
  const double sig_nt = thermal_average_sigmavrel(n, t, T);
  const double sig_pt = thermal_average_sigmavrel(p, t, T);
  std::cout << "<sigma vrel> [mb] for pi+t, n+t, p+t: " << sig_pit << " "
            << sig_nt << " " << sig_pt << std::endl;

  const double sig_pihe3 = thermal_average_sigmavrel(pip, he3, T);
  const double sig_nhe3 = thermal_average_sigmavrel(n, he3, T);
  const double sig_phe3 = thermal_average_sigmavrel(p, he3, T);
  std::cout << "<sigma vrel> [mb] for pi+he3, n+he3, p+he3: " << sig_pihe3
            << " " << sig_nhe3 << " " << sig_phe3 << std::endl;

  const double sig_piH3L = thermal_average_sigmavrel(pip, H3L, T);
  const double sig_nH3L = thermal_average_sigmavrel(n, H3L, T);
  const double sig_pH3L = thermal_average_sigmavrel(p, H3L, T);
  std::cout << "<sigma vrel> [mb] for pi+H3L, n+H3L, p+H3L: " << sig_piH3L
            << " " << sig_nH3L << " " << sig_pH3L << std::endl;

  // 4) Prepare thermal densities at mu = 0
  const double nth_pi = HadronGasEos::partial_density(*pip, T, 0, 0, 0, true) +
                        HadronGasEos::partial_density(*pi0, T, 0, 0, 0, true) +
                        HadronGasEos::partial_density(*pim, T, 0, 0, 0, true),
               nth_p = HadronGasEos::partial_density(*p, T, 0, 0, 0, true),
               nth_n = HadronGasEos::partial_density(*n, T, 0, 0, 0, true),
               nth_La = HadronGasEos::partial_density(*La, T, 0, 0, 0, true),
               nth_d = HadronGasEos::partial_density(*d, T, 0, 0, 0, true),
               nth_t = HadronGasEos::partial_density(*t, T, 0, 0, 0, true),
               nth_he3 = HadronGasEos::partial_density(*he3, T, 0, 0, 0, true),
               nth_H3L = HadronGasEos::partial_density(*H3L, T, 0, 0, 0, true);

  double la_pi = Npi_init / (nth_pi * V), la_p = Np_init / (nth_p * V),
         la_n = Nn_init / (nth_n * V), la_La = NLa_init / (nth_La * V),
         la_d = Nd_init / (nth_d * V), la_t = Nt_init / (nth_t * V),
         la_he3 = Nhe3_init / (nth_he3 * V), la_H3L = NH3L_init / (nth_H3L * V);

  int step_max = std::ceil(tend / dt);
  std::cout << "# t[fm], p, n, La, d, t, he3, H3L, NB, NS" << std::endl;
  for (int step = 0; step < step_max; step++) {
    const double nB = la_p * nth_p + la_n * nth_n + 2 * la_d * nth_d +
                      la_La * nth_La +
                      3 * (la_t * nth_t + la_he3 * nth_he3 + la_H3L * nth_H3L);
    const double nS = (la_H3L * nth_H3L + la_La * nth_La);
    if (step % 10 == 0) {
      std::cout << step * dt << " " << la_p * nth_p * V << " "
                << la_n * nth_n * V << " " << la_La * nth_La * V << " "
                << la_d * nth_d * V << " " << la_t * nth_t * V << " "
                << la_he3 * nth_he3 * V << " " << la_H3L * nth_H3L * V << " "
                << nB * V << " " << nS * V << std::endl;
    }
    // Rates
    double Rd = nth_d *
                (sig_pid * nth_pi * la_pi + sig_pd * nth_p * la_p +
                 sig_nd * nth_n * la_n) *
                fm2_mb,
           Rt = nth_t *
                (sig_pit * nth_pi * la_pi + sig_pt * nth_p * la_p +
                 sig_nt * nth_n * la_n) *
                fm2_mb,
           Rhe3 = nth_he3 *
                  (sig_pihe3 * nth_pi * la_pi + sig_phe3 * nth_p * la_p +
                   sig_nhe3 * nth_n * la_n) *
                  fm2_mb,
           RH3L = nth_H3L *
                  (sig_piH3L * nth_pi * la_pi + sig_pH3L * nth_p * la_p +
                   sig_nH3L * nth_n * la_n) *
                  fm2_mb;
    // d(la)/dt, where la are fugacities
    double dla_d = Rd * (la_p * la_n - la_d);
    double dla_t = Rt * (la_p * la_n * la_n - la_t);
    double dla_he3 = Rhe3 * (la_p * la_p * la_n - la_he3);
    double dla_H3L = RH3L * (la_p * la_n * la_La - la_H3L);

    double dla_p = -dla_d - dla_t - 2 * dla_he3 - dla_H3L;
    double dla_n = -dla_d - 2 * dla_t - dla_he3 - dla_H3L;
    double dla_La = -dla_H3L;

    dla_p *= (dt / nth_p);
    dla_n *= (dt / nth_n);
    dla_La *= (dt / nth_La);
    dla_d *= (dt / nth_d);
    dla_t *= (dt / nth_t);
    dla_he3 *= (dt / nth_he3);
    dla_H3L *= (dt / nth_H3L);

    la_p += dla_p;
    la_n += dla_n;
    la_La += dla_La;
    la_d += dla_d;
    la_t += dla_t;
    la_he3 += dla_he3;
    la_H3L += dla_H3L;
  }
}
