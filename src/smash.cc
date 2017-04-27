/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <getopt.h>

#include <sstream>
#include <vector>

#include <boost/filesystem/fstream.hpp>

#include "include/cxx14compat.h"
#include "include/decaymodes.h"
#include "include/experiment.h"
#include "include/filelock.h"
#include "include/inputfunctions.h"
#include "include/random.h"
#include "include/scatteractionsfinder.h"
/* build dependent variables */
#include "include/config.h"

namespace Smash {

namespace {
/** prints usage information and exits the program
 *
 * \param rc Exit status to return
 * \param progname Name of the program
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
   * (see \ref inputoptions).
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
   * <tr><td>`-m <modus>` <td>`--modus <modus>`
   * <td>This is a shortcut for `-c 'General: { Modus: <modus> }'`. Note that
   *     `-m` always overrides `-c`.
   * <tr><td>`-e <time>` <td>`--endtime <time>`
   * <td>This is a shortcut for `-c 'General: { End_Time: <time> }'`. Note that
   *     `-e` always overrides `-c`.
   * <tr><td>`-o <dir>` <td>`--output <dir>`
   * <td>Sets the output directory. The default output directory is
   *     `./data/<runid>`, where `<rundid>` is an automatically incrementing
   *     integer. Note that this might cause races if several instances of SMASH
   *     run in parallel. In that case, make sure to specify a different output
   *     directory for every instance of SMASH.
   * <tr><td>`-l <dir>` <td>`--list-2-to-n <dir>`
   * <td>Dumps the list of all possible 2->n reactions (n > 1). Note that
   *     resonance decays and formations are NOT dumped. Every particle
   *     available in SMASH is collided against every and reactions with
   *     non-zero cross-section are dumped. Both colliding particles are
   *     assigned momenta from 0.1 to 10 GeV in the opposite directions to
   *     scan the possible sqrt(S).
   * <tr><td>`-r <pdg>` <td>`--resonance <pdg>`
   * <td> Dumps the width(m) and m * spectral function(m^2) versus resonance
   *     mass m.
   * <tr><td>`-f` <td>`--force`
   * <td>Forces overwriting files in the output directory. Normally, if you
   *     specifiy an output directory with `-o`, the directory must be empty.
   *     With `-f` this check is skipped.
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
      "  -l, --list-2-to-n       list all possible 2->2 reactions\n"
      "  -r, --resonance <pdg>   dump width(m) and m*spectral function(m^2)"
      " for resonance pdg\n"
      "  -f, --force             force overwriting files in the output "
      "directory"
      "\n"
      "  -v, --version\n\n");
  std::exit(rc);
}

/** \ingroup exception
 *  Exception class that is thrown if the requested output directory
 * already exists and `-f` was not specified on the command line.
 */
struct OutputDirectoryExists : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
/** \ingroup exception
 *  Exception class that is thrown if no new output path can be
 * generated (there is a directory name for each positive integer
 * value)
 */
struct OutputDirectoryOutOfIds : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

/// returns the default path for output.
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

/// makes sure the output path is valid (throws if not)
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

}  // unnamed namespace

}  // namespace Smash

/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  using namespace Smash;  // NOLINT(build/namespaces)
  setup_default_float_traps();

  const auto &log = logger<LogArea::Main>();

  constexpr option longopts[] = {{"config", required_argument, 0, 'c'},
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
                                 {"version", no_argument, 0, 'v'},
                                 {nullptr, 0, 0, 0}};

  /* strip any path to progname */
  const std::string progname = bf::path(argv[0]).filename().native();

  try {
    bool force_overwrite = false;
    bf::path output_path = default_output_path(), input_path("./config.yaml");
    std::vector<std::string> extra_config;
    char *particles = nullptr, *decaymodes = nullptr, *modus = nullptr,
         *end_time = nullptr, *pdg_string = nullptr;
    // This variable remembers if --list-2-to-n option is activated
    bool list2n_activated = false;
    bool resonance_dump_activated = false;

    /* parse command-line arguments */
    int opt;
    while ((opt = getopt_long(argc, argv, "c:d:e:fhi:m:p:o:lr:v", longopts,
                              nullptr)) != -1) {
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
          break;
        case 'r':
          resonance_dump_activated = true;
          pdg_string = optarg;
          break;
        case 'v':
          std::printf(
              "%s\n"
              "Branch   : %s\nSystem   : %s\nCompiler : %s %s\n"
              "Build    : %s\nDate     : %s\n",
              VERSION_MAJOR, GIT_BRANCH, CMAKE_SYSTEM, CMAKE_CXX_COMPILER_ID,
              CMAKE_CXX_COMPILER_VERSION, CMAKE_BUILD_TYPE, BUILD_DATE);
          std::exit(EXIT_SUCCESS);
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

    /* read in config file */
    Configuration configuration(input_path.parent_path(),
                                input_path.filename());
    for (const auto &config : extra_config) {
      configuration.merge_yaml(config);
    }
    if (particles) {
      std::cout << "particles: " << particles << std::endl;
      if (!bf::exists(particles)) {
        std::stringstream err;
        err << "The particles file was expected at '" << particles
            << "', but the file does not exist.";
        throw std::runtime_error(err.str());
      }
      configuration["particles"] = read_all(bf::ifstream{particles});
    }
    if (decaymodes) {
      if (!bf::exists(decaymodes)) {
        std::stringstream err;
        err << "The decay modes file was expected at '" << decaymodes
            << "', but the file does not exist.";
        throw std::runtime_error(err.str());
      }
      configuration["decaymodes"] = read_all(bf::ifstream{decaymodes});
    }
    if (list2n_activated) {
      // Do not make all elastic cross-sections a fixed number
      constexpr float elastic_parameter = -1.0f;
      // Does not matter here, just dummy
      constexpr int ntest = 1;
      // Print only 2->n, n > 1. Do not dump decays, which can be found in
      // decaymodes.txt anyway
      constexpr bool two_to_one = false;
      ParticleType::create_type_list(configuration.take({"particles"}));
      DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
      std::vector<bool> nucleon_has_interacted = {};
      auto scat_finder = make_unique<ScatterActionsFinder>(elastic_parameter,
                                     ntest, nucleon_has_interacted, two_to_one);
      scat_finder->dump_reactions();
      std::exit(EXIT_SUCCESS);
    }
    if (resonance_dump_activated) {
      ParticleType::create_type_list(configuration.take({"particles"}));
      DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
      PdgCode pdg(pdg_string);
      const ParticleType &res = ParticleType::find(pdg);
      res.dump_width_and_spectral_function();
      std::exit(EXIT_SUCCESS);
    }
    if (modus) {
      configuration["General"]["Modus"] = std::string(modus);
    }
    if (end_time) {
      configuration["General"]["End_Time"] = std::abs(std::atof(end_time));
    }

    /* set up logging */
    set_default_loglevel(
        configuration.take({"Logging", "default"}, einhard::ALL));
    create_all_loggers(configuration["Logging"]);
    log.info(progname, " (", VERSION_MAJOR, ')');

    int64_t seed = configuration.read({"General", "Randomseed"});
    if (seed < 0) {
      // Seed with a real random value, if available
      std::random_device rd;
      seed = rd();
      configuration["General"]["Randomseed"] = seed;
    }

    /* check output path*/
    ensure_path_is_valid(output_path);
    const bf::path lock_path = output_path / "smash.lock";
    FileLock lock(lock_path);
    if (!lock.acquire()) {
      throw std::runtime_error(
          "Another instance of SMASH is already writing to the specified "
          "output directory. If you are sure this is not the case, remove \"" +
          lock_path.native() + "\".");
    }
    log.debug("output path: ", output_path);
    if (!force_overwrite && bf::exists(output_path / "config.yaml")) {
      throw std::runtime_error(
          "Output directory would get overwritten. Select a different output "
          "directory, clean up, or tell SMASH to ignore existing files.");
    }

    // keep a copy of the configuration that was used in the output directory
    bf::ofstream(output_path / "config.yaml") << configuration.to_string()
                                              << '\n';

    // take the seed setting only after the configuration was stored to file
    seed = configuration.take({"General", "Randomseed"});
    Random::set_seed(seed);
    log.info() << "Random number seed: " << seed;

    log.trace(source_location, " create ParticleType and DecayModes");
    ParticleType::create_type_list(configuration.take({"particles"}));
    DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
    ParticleType::check_consistency();

    // create an experiment
    log.trace(source_location, " create Experiment");
    auto experiment = ExperimentBase::create(configuration, output_path);
    const std::string report = configuration.unused_values_report();
    if (report != "{}") {
      log.warn() << "The following configuration values were not used:\n"
                 << report;
    }

    // run the experiment
    log.trace(source_location, " run the Experiment");
    experiment->run();
  } catch (std::exception &e) {
    log.fatal() << "SMASH failed with the following error:\n" << e.what();
    return EXIT_FAILURE;
  }

  log.trace() << source_location << " about to return from main";
  return 0;
}
