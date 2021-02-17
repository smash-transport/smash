/*
 *
 *    Copyright (c) 2012-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <getopt.h>

#include <set>
#include <sstream>
#include <vector>

#include <boost/filesystem/fstream.hpp>

#include "smash/cxx14compat.h"
#include "smash/decaymodes.h"
#include "smash/experiment.h"
#include "smash/filelock.h"
#include "smash/random.h"
#include "smash/scatteractionsfinder.h"
#include "smash/setup_particles_decaymodes.h"
#include "smash/sha256.h"
#include "smash/stringfunctions.h"
/* build dependent variables */
#include "smash/config.h"

namespace smash {

namespace {
/**
 * Prints usage information and exits the program
 *
 * \param[out] rc Exit status to return
 * \param[in] progname Name of the program
 *
 * usage() is called when either the `--help` or `-h` command line
 * options are given to the program; in this case, the exit status is
 * EXIT_SUCCESS, or when an unknown option is given; in this case,
 * the exit status is EXIT_FAIL.
 */
void usage(const int rc, const std::string &progname) {
  /*!\Userguide
   * \page page_smash_invocation SMASH Invocation
   *
   * SMASH can be run simply by executing the binary without any options (i.e.
   * there are no required arguments). It does require an input file, though
   * (see \ref inputconfig).
   * Per default, the input file is expected in the current working directory
   * with the name '`config.yaml`'.
   *
   * The following options are supported:
   *
   * <table>
   * <tr><th>Short&nbsp;Variant <th>Long&nbsp;Variant <th>Documentation
   * <tr><td>`-h` <td>`--help`
   * <td>Prints usage information and quits the program.
   * <tr><td>`-v` <td>`--version`
   * <td>Prints the version of SMASH and quits the program.
   * <tr><td>`-i <file>` <td>`--inputfile <file>`
   * <td>Overrides the location of the default '`config.yaml`' input file. The
   *     input settings will be read from the specified file instead.
   * <tr><td>`-d <file>` <td>`--decaymodes <file>`
   * <td>The default decay modes are compiled in. With this argument you can
   *     override the decay modes to the exact set defined in the file. Multiple
   *     `-d` arguments are not supported.
   * <tr><td>`-p <file>` <td>`--particles <file>`
   * <td>The default particle data is compiled in. With this argument you can
   *     override the particles to the exact set defined in the file. Multiple
   *     `-p` arguments are not supported.
   * <tr><td>`-c <YAML string>` <td>`--config <YAML string>`
   * <td>The string argument to `-c` containts YAML markup to override input
   *     options from the input file (`-i`). Multiple `-c` arguments are
   *     supported. (Later settings may override preceding settings.) This can
   *     be a handy way to test different scenarios from a script.
   * <tr><td>`-m \<modus\>` <td>`--modus \<modus\>`
   * <td>This is a shortcut for `-c 'General: { Modus: \<modus\> }'`. Note that
   *     `-m` always overrides `-c`.
   * <tr><td>`-e \<time\>` <td>`--endtime \<time\>`
   * <td>This is a shortcut for `-c 'General: { End_Time: \<time\> }'`. Note
   * that
   *     `-e` always overrides `-c`.
   * <tr><td>`-o \<dir\>` <td>`--output \<dir\>`
   * <td>Sets the output directory. The default output directory is
   *     `./data/<runid>`, where `<rundid>` is an automatically incrementing
   *     integer. Note that this might cause races if several instances of SMASH
   *     run in parallel. In that case, make sure to specify a different output
   *     directory for every instance of SMASH.
   * <tr><td>`-l \<dir\>` <td>`--list-2-to-n \<dir\>`
   * <td>Dumps the list of all possible 2->n reactions (n > 1). Note that
   *     resonance decays and formations are NOT dumped. Every particle
   *     available in SMASH is collided against every and reactions with
   *     non-zero cross-section are dumped. Both colliding particles are
   *     assigned momenta from 0.1 to 10 GeV in the opposite directions to
   *     scan the possible sqrt(S).
   * <tr><td>`-r <pdg>` <td>`--resonance <pdg>`
   * <td> Dumps the width(m) and m * spectral function(m^2) versus resonance
   *     mass m.
   * <tr><td>`-s <pdg1>,<pdg2>[,mass1,mass2]`
   * <td>`--cross-sections <pdg1>,<pdg2>[,mass1,mass2[,plab1,...]]`
   * <td> Dumps all the partial cross-sections of pdg1 + pdg2 with
   *     masses mass1 and mass2. Masses are optional, default values are pole
   *     masses. Optionally, the lab frame momenta (fixed target) in GeV can be
   *     specified. (The value of plab depends on the order of the particles.
   *     The first particle is considered to be the projectile, the second one
   *     the target.)
   * <tr><td>`-f` <td>`--force`
   * <tr><td>`-S <pdg1>,<pdg2>[,mass1,mass2]`
   * <td>`--cross-sections-fs <pdg1>,<pdg2>[,mass1,mass2[,plab1,...]]`
   * <td> Dumps an approximation of the final-state cross-sections of pdg1 +
   *     pdg2 with masses mass1 and mass2. Masses are optional, default values
   *     are pole masses. Optionally, the lab frame momenta (fixed target) in
   *     GeV can be specified. (The value of plab depends on the order of the
   *     particles. The first is considered to be the projectile, the second
   *     one the target.) After the initial collision, only decays are
   *     considered and all resonances are assumed to have their pole mass. This
   *     may yield different results than a full simulation with SMASH, where
   *     the resonances masses are sampled from the spectral function.
   *     Typically, this results in errors of less than 1 mb in the worst case.
   *     Also, contributions from strings are not considered.
   * <tr><td>`-f` <td>`--force`
   * <td>Forces overwriting files in the output directory. Normally, if you
   *     specifiy an output directory with `-o`, the directory must be empty.
   *     With `-f` this check is skipped.
   * <tr><td>`-q` <td>`--quiet`
   * <td> Quiets the disclaimer for scenarios where no printout is wanted. To
   * get no printout, you also need to disable logging from the config.
   * </table>
   */
  std::printf("\nUsage: %s [option]\n\n", progname.c_str());
  std::printf(
      "  -h, --help              usage information\n"
      "\n"
      "  -i, --inputfile <file>  path to input configuration file\n"
      "                          (default: ./config.yaml)\n"
      "  -d, --decaymodes <file> override default decay modes from file\n"
      "  -p, --particles <file>  override default particles from file\n"
      "\n"
      "  -c, --config <YAML>     specify config value overrides\n"
      "                          (multiple -c arguments are supported)\n"
      "  -m, --modus <modus>     shortcut for -c 'General: { Modus: <modus> "
      "}'\n"
      "  -e, --endtime <time>    shortcut for -c 'General: { End_Time: <time> "
      "}'"
      "\n"
      "\n"
      "  -o, --output <dir>      output directory (default: ./data/<runid>)\n"
      "  -l, --list-2-to-n       list all possible 2->n reactions (with n>1)\n"
      "  -r, --resonance <pdg>   dump width(m) and m*spectral function(m^2)"
      " for resonance pdg\n"
      "  -s, --cross-sections    <pdg1>,<pdg2>[,mass1,mass2[,plab1,...]] \n"
      "                          dump all partial cross-sections of "
      "pdg1 + pdg2 reactions versus sqrt(s).\n"
      "  -S, --cross-sections-fs <pdg1>,<pdg2>[,mass1,mass2[,plab1,...]] \n"
      "                          dump an approximation of the final-state"
      " cross-sections\n"
      "                          of pdg1 + pdg2 reactions versus sqrt(s).\n"
      "                          Contributions from strings are not considered"
      " for the final state.\n"
      "                          Masses are optional, by default pole masses"
      " are used.\n"
      "                          Note the required comma and no spaces.\n"
      "  -f, --force             force overwriting files in the output "
      "directory"
      "\n"
      "  -x, --dump_iSS          Dump particle table in iSS format\n"
      "                          This format is used in MUSIC and CLVisc\n"
      "                          relativistic hydro codes\n"
      "  -q, --quiet             Supress disclaimer print-out\n"
      "  -n, --no-cache          Don't cache integrals on disk\n"
      "  -v, --version\n\n");
  std::exit(rc);
}

/// Print the disclaimer.
void print_disclaimer() {
  std::cout
      << "###################################################################"
      << "############"
      << "\n"
      << "\n"
      << "            :::s.\n"
      << "       ....  ''ss::                                                  "
         "ah:\n"
      << "     a::''.   ..:sss                                                "
         ".HH.\n"
      << "   mss'..     ...s'm.   sSA##As  mAhh##hsh##s.   .hA##ASa  sS###As  "
         "hMAh##h.\n"
      << "  :a':'.      .'':as   sM#'     .HHm''HMS''AMs  mMA' .AMs aMA.     "
         "'HHa..HM:\n"
      << "  .a:s'.    ..''ss     'h#H#S.  mMm  'M#' .HH. 'MH.  :MA  'h#H#S.  "
         "hMm  :M#\n"
      << "   .::ss'.  ....          'SMm .HH.  SMa  hMm  sM#..mHMs     'AMa "
         "'MH.  SMa\n"
      << "      .s::'           #SMASHh  aMS  .MH:  HM:   #MMMmMM. #SMASHh  "
         "aMS  .MM.\n"
      << "\n"
      << "###################################################################"
      << "############"
      << "\n"
      << " This is SMASH version: " << VERSION_MAJOR << "\n"
      << " Simulating Many Accelerated Strongly-interacting Hadrons"
      << "\n"
      << "\n"
      << " Distributed under the GNU General Public License 3.0"
      << " (GPLv3 or later)."
      << "\n"
      << " See LICENSE for details."
      << "\n"
      << " For the full list of contributors see AUTHORS."
      << "\n"
      << "\n"
      << " When using SMASH, please cite"
      << "\n"
      << "      J. Weil et al., Phys.Rev. C94 (2016) no.5, 054905"
      << "\n"
      << " together with the software DOI for the specific code version "
      << "employed:"
      << "\n"
      << "      https://doi.org/10.5281/zenodo.3484711."
      << "\n"
      << " In addition, if Pythia is used please cite"
      << "\n"
      << "      T. Sjöstrand, S. Mrenna and P. Skands, JHEP05 (2006) 026,"
      << "\n"
      << "              Comput. Phys. Comm. 178 (2008) 852."
      << "\n"
      << "\n"
      << " Webpage: https://smash-transport.github.io"
      << "\n"
      << "\n"
      << " Report issues at https://github.com/smash-transport/smash"
      << "\n"
      << " or contact us by email at elfner@itp.uni-frankfurt.de"
      << "\n"
      << "\n"
      << "###################################################################"
      << "############"
      << "\n"
      << "\n";
}

/**
 * \ingroup exception
 * Exception class that is thrown, if the requested output directory
 * already exists and `-f` was not specified on the command line.
 */
struct OutputDirectoryExists : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
/**
 * \ingroup exception
 * Exception class that is thrown, if no new output path can be
 * generated (there is a directory name for each positive integer
 * value)
 */
struct OutputDirectoryOutOfIds : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

/// \return the default path for output.
bf::path default_output_path() {
  const bf::path p = bf::absolute("data");
  if (!bf::exists(p)) {
    return p / "0";
  }
  bf::path p2;
  for (int id = 0; id < std::numeric_limits<int>::max(); ++id) {
    p2 = p / std::to_string(id);
    if (!bf::exists(p2)) {
      break;
    }
  }
  if (p == p2) {
    throw OutputDirectoryOutOfIds("no unique data subdir ID left");
  }
  return p2;
}

/**
 * Ensures the output path is valid.
 *
 * \throw OutputDirectoryExists if the Output directory already exists.
 * \param[in] path The output path to be written to
 */
void ensure_path_is_valid(const bf::path &path) {
  if (bf::exists(path)) {
    if (!bf::is_directory(path)) {
      throw OutputDirectoryExists("The given path (" + path.native() +
                                  ") exists, but it is not a directory.");
    }
  } else {
    if (!bf::create_directories(path)) {
      throw OutputDirectoryExists(
          "Race condition detected: The directory " + path.native() +
          " did not exist a few cycles ago, but was created in the meantime by "
          "a different process.");
    }
  }
}

/**
 * Prepares ActionsFinder for cross-section and reaction dumps.
 *
 * \param[in] configuration Necessary parameters to switch reactions on/off
 * \return The constructed Scatteractionsfinder.
 */
ScatterActionsFinder actions_finder_for_dump(Configuration configuration) {
  // The following parameters are only relevant for nucleus-nucleus collisions.
  // Setting them to the valuing values makes sure they don't have any effect.
  const std::vector<bool> nucleon_has_interacted = {};
  const int N_tot = 0;
  const int N_proj = 0;

  ExperimentParameters params = create_experiment_parameters(configuration);
  return ScatterActionsFinder(configuration, params, nucleon_has_interacted,
                              N_tot, N_proj);
}

/** Checks if the SMASH version is compatible with the version of the
 * configuration file
 *
 * \param[in] configuration The configuration object
 *
 * \throws Runtime error if versions do not match or if config version is
 * invalid
 */
void check_config_version_is_compatible(Configuration configuration) {
  const std::string smash_version = "1.8";
  const std::set<std::string> compatible_config_versions = {"1.8"};

  const std::string config_version = configuration.read({"Version"});

  if (compatible_config_versions.find(config_version) ==
      compatible_config_versions.end()) {
    std::stringstream err;
    err << "The version of the configuration file (" << config_version
        << ") is not compatible with the SMASH version (" << smash_version
        << ").\nThe following config versions are supported:\n";
    for (auto it : compatible_config_versions) {
      err << it << " ";
    }
    err << "\nPlease consider updating your config or using a compatible SMASH"
           " version.";
    throw std::runtime_error(err.str());
  }
}

/**
 * Checks if there are unused config values.
 */
void check_for_unused_config_values(const Configuration &configuration) {
  const std::string report = configuration.unused_values_report();

  if (report != "{}") {
    throw std::runtime_error(
        "The following configuration values were not used:\n" + report);
  }
}

/**
 * Remove all config values that are only needed for simulations.
 *
 * This is useful when checking for unused config value when SMASH only
 * outputs cross sections, resonance properties or possible reactions.
 */
void ignore_simulation_config_values(Configuration &configuration) {
  for (const std::string s :
       {"Version", "particles", "decaymodes", "Modi", "General", "Output",
        "Lattice", "Potentials", "Forced_Thermalization"}) {
    if (configuration.has_value({s.c_str()})) {
      configuration.take({s.c_str()});
    }
  }
}

/// Initialize the particles and decays from the configuration.
void initialize_particles_and_decays(Configuration &configuration) {
  ParticleType::create_type_list(configuration.take({"particles"}));
  DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
  ParticleType::check_consistency();
}

/** Initialize the particles and decays from the configuration,
 *  the hash and the path to the cashed resonance integrals.
 */
void initialize_particles_and_decays(Configuration &configuration,
                                     sha256::Hash hash,
                                     bf::path tabulations_path) {
  initialize_particles_and_decays(configuration);
  logg[LMain].info("Tabulating cross section integrals...");
  IsoParticleType::tabulate_integrals(hash, tabulations_path);
}

}  // unnamed namespace

}  // namespace smash

/**
 * Main program
 * Smashes Many Accelerated Strongly-Interacting Hadrons :-)
 *
 * Do command line parsing and hence decides modus
 *
 * \param[in] argc Number of arguments on command-line
 * \param[in] argv List of arguments on command-line
 * \return Either 0 or EXIT_FAILURE.
 */
int main(int argc, char *argv[]) {
  using namespace smash;  // NOLINT(build/namespaces)

  constexpr option longopts[] = {
      {"config", required_argument, 0, 'c'},
      {"decaymodes", required_argument, 0, 'd'},
      {"endtime", required_argument, 0, 'e'},
      {"force", no_argument, 0, 'f'},
      {"help", no_argument, 0, 'h'},
      {"inputfile", required_argument, 0, 'i'},
      {"modus", required_argument, 0, 'm'},
      {"particles", required_argument, 0, 'p'},
      {"output", required_argument, 0, 'o'},
      {"list-2-to-n", no_argument, 0, 'l'},
      {"resonance", required_argument, 0, 'r'},
      {"cross-sections", required_argument, 0, 's'},
      {"cross-sections-fs", required_argument, 0, 'S'},
      {"dump-iSS", no_argument, 0, 'x'},
      {"version", no_argument, 0, 'v'},
      {"no-cache", no_argument, 0, 'n'},
      {"quiet", no_argument, 0, 'q'},
      {nullptr, 0, 0, 0}};

  // strip any path to progname
  const std::string progname = bf::path(argv[0]).filename().native();

  try {
    bool force_overwrite = false;
    bf::path output_path = default_output_path(), input_path("./config.yaml");
    std::vector<std::string> extra_config;
    char *particles = nullptr, *decaymodes = nullptr, *modus = nullptr,
         *end_time = nullptr, *pdg_string = nullptr, *cs_string = nullptr;
    bool list2n_activated = false;
    bool resonance_dump_activated = false;
    bool cross_section_dump_activated = false;
    bool final_state_cross_sections = false;
    bool particles_dump_iSS_format = false;
    bool cache_integrals = true;

    // parse command-line arguments
    int opt;
    bool suppress_disclaimer = false;
    while ((opt = getopt_long(argc, argv, "c:d:e:fhi:m:p:o:lr:s:S:xvnq",
                              longopts, nullptr)) != -1) {
      switch (opt) {
        case 'c':
          extra_config.emplace_back(optarg);
          break;
        case 'd':
          decaymodes = optarg;
          break;
        case 'f':
          force_overwrite = true;
          break;
        case 'i':
          input_path = optarg;
          break;
        case 'h':
          usage(EXIT_SUCCESS, progname);
          break;
        case 'm':
          modus = optarg;
          break;
        case 'p':
          particles = optarg;
          break;
        case 'e':
          end_time = optarg;
          break;
        case 'o':
          output_path = optarg;
          break;
        case 'l':
          list2n_activated = true;
          suppress_disclaimer = true;
          break;
        case 'r':
          resonance_dump_activated = true;
          pdg_string = optarg;
          suppress_disclaimer = true;
          break;
        case 'S':
          final_state_cross_sections = true;
          // fallthrough
        case 's':
          cross_section_dump_activated = true;
          cs_string = optarg;
          suppress_disclaimer = true;
          break;
        case 'x':
          particles_dump_iSS_format = true;
          suppress_disclaimer = true;
          break;
        case 'v':
          std::printf(
              "%s\n"
              "Branch   : %s\nSystem   : %s\nCompiler : %s %s\n"
              "Build    : %s\nDate     : %s\n",
              VERSION_MAJOR, GIT_BRANCH, CMAKE_SYSTEM, CMAKE_CXX_COMPILER_ID,
              CMAKE_CXX_COMPILER_VERSION, CMAKE_BUILD_TYPE, BUILD_DATE);
          std::exit(EXIT_SUCCESS);
        case 'n':
          cache_integrals = false;
          break;
        case 'q':
          suppress_disclaimer = true;
          break;
        default:
          usage(EXIT_FAILURE, progname);
      }
    }

    // Abort if there are unhandled arguments left.
    if (optind < argc) {
      std::cout << argv[0] << ": invalid argument -- '" << argv[optind]
                << "'\n";
      usage(EXIT_FAILURE, progname);
    }

    if (!suppress_disclaimer) {
      print_disclaimer();
    }

    // Read in config file
    Configuration configuration(input_path.parent_path(),
                                input_path.filename());
    // Merge config passed via command line
    for (const auto &config : extra_config) {
      configuration.merge_yaml(config);
    }

    // Set up logging
    set_default_loglevel(
        configuration.take({"Logging", "default"}, einhard::ALL));
    create_all_loggers(configuration["Logging"]);

    setup_default_float_traps();

    // check if version matches before doing anything else
    check_config_version_is_compatible(configuration);

    logg[LMain].trace(SMASH_SOURCE_LOCATION, " create ParticleType and DecayModes");

    auto particles_and_decays =
        load_particles_and_decaymodes(particles, decaymodes);
    /* For particles and decaymodes: external file is superior to config.
     * Hovever, warn in case of conflict.
     */
    if (configuration.has_value({"particles"}) && particles) {
      logg[LMain].warn(
          "Ambiguity: particles from external file ", particles,
          " requested, but there is also particle list in the config."
          " Using particles from ",
          particles);
    }
    if (!configuration.has_value({"particles"}) || particles) {
      configuration["particles"] = particles_and_decays.first;
    }

    if (configuration.has_value({"decaymodes"}) && decaymodes) {
      logg[LMain].warn(
          "Ambiguity: decaymodes from external file ", decaymodes,
          " requested, but there is also decaymodes list in the config."
          " Using decaymodes from",
          decaymodes);
    }
    if (!configuration.has_value({"decaymodes"}) || decaymodes) {
      configuration["decaymodes"] = particles_and_decays.second;
    }

    // Calculate a hash of the SMASH version, the particles and decaymodes.
    const std::string version(VERSION_MAJOR);
    const std::string particle_string = configuration["particles"].to_string();
    const std::string decay_string = configuration["decaymodes"].to_string();
    sha256::Context hash_context;
    hash_context.update(version);
    hash_context.update(particle_string);
    hash_context.update(decay_string);
    const auto hash = hash_context.finalize();
    logg[LMain].info() << "Config hash: " << sha256::hash_to_string(hash);

    bf::path tabulations_path;
    if (cache_integrals) {
      tabulations_path = output_path.parent_path() / "tabulations";
      bf::create_directories(tabulations_path);
      logg[LMain].info() << "Tabulations path: " << tabulations_path;
    } else {
      tabulations_path = "";
    }
    if (list2n_activated) {
      /* Print only 2->n, n > 1. Do not dump decays, which can be found in
       * decaymodes.txt anyway */
      configuration.merge_yaml("{Collision_Term: {Two_to_One: False}}");
      logg[LMain].info() << "Tabulations path: " << tabulations_path;
      initialize_particles_and_decays(configuration, hash, tabulations_path);
      auto scat_finder = actions_finder_for_dump(configuration);

      ignore_simulation_config_values(configuration);
      check_for_unused_config_values(configuration);

      scat_finder.dump_reactions();
      std::exit(EXIT_SUCCESS);
    }
    if (particles_dump_iSS_format) {
      initialize_particles_and_decays(configuration);
      ParticleTypePtrList list;
      list.clear();
      for (const auto &ptype : ParticleType::list_all()) {
        list.push_back(&ptype);
      }
      std::sort(list.begin(), list.end(),
                [](ParticleTypePtr a, ParticleTypePtr b) {
                  return a->mass() < b->mass();
                });
      for (const ParticleTypePtr &ptype : list) {
        if (ptype->pdgcode().is_lepton() || ptype->baryon_number() < 0) {
          continue;
        }
        const auto &decay_modes = ptype->decay_modes();
        const auto &modelist = decay_modes.decay_mode_list();
        int ndecays = ptype->is_stable() ? 1 : modelist.size();
        std::printf("%13i %s %10.5f %10.5f %5i %5i %5i %5i %5i %5i %5i %5i\n",
                    ptype->pdgcode().get_decimal(),
                    smash::utf8::fill_left(ptype->name(), 12, ' ').c_str(),
                    ptype->mass(), ptype->width_at_pole(),
                    ptype->pdgcode().spin_degeneracy(), ptype->baryon_number(),
                    ptype->strangeness(), ptype->pdgcode().charmness(),
                    ptype->pdgcode().bottomness(), ptype->isospin() + 1,
                    ptype->charge(), ndecays);
        if (!ptype->is_stable()) {
          for (const auto &decay : modelist) {
            auto ptypes = decay->particle_types();
            std::printf("%13i %13i %20.5f %13i %13i %13i %13i %13i\n",
                        ptype->pdgcode().get_decimal(), 2, decay->weight(),
                        ptypes[0]->pdgcode().get_decimal(),
                        ptypes[1]->pdgcode().get_decimal(), 0, 0, 0);
          }
        } else {
          std::printf("%13i %13i %20.5f %13i %13i %13i %13i %13i\n",
                      ptype->pdgcode().get_decimal(), 1, 1.0,
                      ptype->pdgcode().get_decimal(), 0, 0, 0, 0);
        }
      }
      std::exit(EXIT_SUCCESS);
    }
    if (resonance_dump_activated) {
      // Ignore config values that don't make sense.
      initialize_particles_and_decays(configuration, hash, tabulations_path);
      const auto _dummy = ExperimentBase::create(configuration, "");
      ignore_simulation_config_values(configuration);
      check_for_unused_config_values(configuration);

      PdgCode pdg(pdg_string);
      const ParticleType &res = ParticleType::find(pdg);
      res.dump_width_and_spectral_function();
      std::exit(EXIT_SUCCESS);
    }
    if (cross_section_dump_activated) {
      initialize_particles_and_decays(configuration, hash, tabulations_path);
      std::string arg_string(cs_string);
      std::vector<std::string> args = split(arg_string, ',');
      const unsigned int n_arg = args.size();
      if (n_arg != 2 && n_arg != 4 && n_arg < 5) {
        throw std::invalid_argument("-s usage: pdg1,pdg2[,m1,m2[,sqrts1,...]]");
      }
      PdgCode pdg_a(args[0]), pdg_b(args[1]);
      const ParticleType &a = ParticleType::find(pdg_a);
      const ParticleType &b = ParticleType::find(pdg_b);
      if (n_arg < 4) {
        for (unsigned int i = 0; i < 4 - n_arg; i++) {
          args.push_back("");
        }
      }
      double ma = (args[2] == "") ? a.mass() : std::stod(args[2]);
      double mb = (args[3] == "") ? b.mass() : std::stod(args[3]);
      if (a.is_stable() && args[2] != "" && std::stod(args[2]) != a.mass()) {
        ma = a.mass();
        std::cerr << "Warning: pole mass is used for stable particle "
                  << a.name() << " instead of " << args[2] << std::endl;
      }
      if (b.is_stable() && args[3] != "" && std::stod(args[3]) != b.mass()) {
        mb = b.mass();
        std::cerr << "Warning: pole mass is used for stable particle "
                  << b.name() << " instead of " << args[3] << std::endl;
      }
      const size_t plab_size = n_arg <= 4 ? 0 : n_arg - 4;
      std::vector<double> plab;
      plab.reserve(plab_size);
      for (size_t i = 4; i < n_arg; i++) {
        plab.push_back(std::stod(args.at(i)));
      }
      auto scat_finder = actions_finder_for_dump(configuration);

      ignore_simulation_config_values(configuration);
      check_for_unused_config_values(configuration);

      scat_finder.dump_cross_sections(a, b, ma, mb, final_state_cross_sections,
                                      plab);
      std::exit(EXIT_SUCCESS);
    }
    if (modus) {
      configuration["General"]["Modus"] = std::string(modus);
    }
    if (end_time) {
      configuration["General"]["End_Time"] = std::abs(std::atof(end_time));
    }

    int64_t seed = configuration.read({"General", "Randomseed"});
    if (seed < 0) {
      configuration["General"]["Randomseed"] = random::generate_63bit_seed();
    }

    // Check output path
    ensure_path_is_valid(output_path);
    const bf::path lock_path = output_path / "smash.lock";
    FileLock lock(lock_path);
    if (!lock.acquire()) {
      throw std::runtime_error(
          "Another instance of SMASH is already writing to the specified "
          "output directory. If you are sure this is not the case, remove \"" +
          lock_path.native() + "\".");
    }
    logg[LMain].debug("output path: ", output_path);
    if (!force_overwrite && bf::exists(output_path / "config.yaml")) {
      throw std::runtime_error(
          "Output directory would get overwritten. Select a different output "
          "directory, clean up, or tell SMASH to ignore existing files.");
    }

    /* Keep a copy of the configuration that was used in the output directory
     * also save information about SMASH build as a comment */
    bf::ofstream(output_path / "config.yaml")
        << "# " << VERSION_MAJOR << '\n'
        << "# Branch   : " << GIT_BRANCH << '\n'
        << "# System   : " << CMAKE_SYSTEM << '\n'
        << "# Compiler : " << CMAKE_CXX_COMPILER_ID << ' '
        << CMAKE_CXX_COMPILER_VERSION << '\n'
        << "# Build    : " << CMAKE_BUILD_TYPE << '\n'
        << "# Date     : " << BUILD_DATE << '\n'
        << configuration.to_string() << '\n';
    logg[LMain].trace(SMASH_SOURCE_LOCATION, " create ParticleType and DecayModes");
    initialize_particles_and_decays(configuration, hash, tabulations_path);

    // Create an experiment
    logg[LMain].trace(SMASH_SOURCE_LOCATION, " create Experiment");
    auto experiment = ExperimentBase::create(configuration, output_path);

    // Version value is not used in experiment. Get rid of it to prevent
    // warning.
    configuration.take({"Version"});
    check_for_unused_config_values(configuration);

    // Run the experiment
    logg[LMain].trace(SMASH_SOURCE_LOCATION, " run the Experiment");
    experiment->run();
  } catch (std::exception &e) {
    logg[LMain].fatal() << "SMASH failed with the following error:\n"
                        << e.what();
    return EXIT_FAILURE;
  }

  logg[LMain].trace() << SMASH_SOURCE_LOCATION << " about to return from main";
  return 0;
}
