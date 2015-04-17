/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <fenv.h>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <list>
#include <stdexcept>
#include <string>

#include "include/configuration.h"
#include "include/decaymodes.h"
#include "include/experiment.h"
#include "include/forwarddeclarations.h"
#include "include/inputfunctions.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particletype.h"
/* Outputs */
#include "include/binaryoutputcollisions.h"
#include "include/binaryoutputparticles.h"
#include "include/densityoutput.h"
#include "include/oscaroutput.h"
#ifdef SMASH_USE_ROOT
#  include "include/rootoutput.h"
#endif
#include "include/vtkoutput.h"
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
   *     integer.
   * <tr><td>`-f` <td>`--force`
   * <td>Forces overwriting files in the output directory. Normally, if you
   *     specifiy an output directory with `-o`, the directory must be empty.
   *     With `-f` this check is skipped.
   * </table>
   */
  printf("\nUsage: %s [option]\n\n", progname.c_str());
  printf("Calculate transport box\n"
    "  -h, --help              usage information\n"
    "\n"
    "  -i, --inputfile <file>  path to input configuration file\n"
    "                          (default: ./config.yaml)\n"
    "  -d, --decaymodes <file> override default decay modes from file\n"
    "  -p, --particles <file>  override default particles from file\n"
    "\n"
    "  -c, --config <YAML>     specify config value overrides\n"
    "                          (multiple -c arguments are supported)\n"
    "  -m, --modus <modus>     shortcut for -c 'General: { Modus: <modus> }'\n"
    "  -e, --endtime <time>    shortcut for -c 'General: { End_Time: <time> }'"
    "\n"
    "\n"
    "  -o, --output <dir>      output directory (default: ./data/<runid>)\n"
    "  -f, --force             force overwriting files in the output directory"
    "\n"
    "  -v, --version\n\n");
  exit(rc);
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

void setup_floatingpoint_traps() {
  const auto &log = logger<LogArea::Main>();

  // standard C/C++ don't have a function to modify the trapping behavior. You
  // can only save and restore the setup. With glibc you can change it via
  // feenableexcept and fedisableexcept.
#ifdef _GNU_SOURCE
  // pole error occurred in a floating-point operation:
  if (-1 == feenableexcept(FE_DIVBYZERO)) {
    log.warn("Failed to setup trap on pole error.");
  }

  // domain error occurred in an earlier floating-point operation:
  if (-1 == feenableexcept(FE_INVALID)) {
    log.warn("Failed to setup trap on domain error.");
  }

  // the result of the earlier floating-point operation was too large to be
  // representable:
  if (-1 == feenableexcept(FE_OVERFLOW)) {
    log.warn("Failed to setup trap on overflow.");
  }

  // the result of the earlier floating-point operation was subnormal with a
  // loss of precision:
  if (-1 == feenableexcept(FE_UNDERFLOW)) {
    log.warn("Failed to setup trap on underflow.");
  }

  // there's also FE_INEXACT, but this traps if "rounding was necessary to store
  // the result of an earlier floating-point operation". This is common and not
  // really an error condition.
#else
  log.warn("Failed to setup floating-point traps.");
#endif
}

}  // unnamed namespace

}  // namespace Smash

/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  using namespace Smash;
  setup_floatingpoint_traps();

  const auto &log = logger<LogArea::Main>();

  constexpr option longopts[] = {
    { "config",     required_argument,      0, 'c' },
    { "decaymodes", required_argument,      0, 'd' },
    { "endtime",    required_argument,      0, 'e' },
    { "force",      no_argument,            0, 'f' },
    { "help",       no_argument,            0, 'h' },
    { "inputfile",  required_argument,      0, 'i' },
    { "modus",      required_argument,      0, 'm' },
    { "particles",  required_argument,      0, 'p' },
    { "output",     required_argument,      0, 'o' },
    { "version",    no_argument,            0, 'v' },
    { nullptr,      0,                      0,  0  }
  };

  /* strip any path to progname */
  const std::string progname = bf::path(argv[0]).filename().native();

  try {
    bool force_overwrite = false;
    bf::path output_path = default_output_path(),
             input_path("./config.yaml");
    std::vector<std::string> extra_config;
    char *particles = nullptr, *decaymodes = nullptr, *modus = nullptr,
         *end_time = nullptr;

    /* parse command-line arguments */
    int opt;
    while ((opt = getopt_long(argc, argv, "c:d:e:fhi:m:p:o:v", longopts,
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
        case 'v':
          printf("%s\nSystem   : %s\nCompiler : %s %s\n"
                     "Build    : %s\nDate     : %s\n",
                 VERSION_MAJOR, CMAKE_SYSTEM, CMAKE_CXX_COMPILER_ID,
                 CMAKE_CXX_COMPILER_VERSION, CMAKE_BUILD_TYPE, BUILD_DATE);
          exit(EXIT_SUCCESS);
        default:
          usage(EXIT_FAILURE, progname);
      }
    }

    /* read in config file */
    Configuration configuration(input_path.parent_path(),
                                input_path.filename());
    for (const auto &config : extra_config) {
      configuration.merge_yaml(config);
    }
    if (particles)
      configuration["particles"] = read_all(bf::ifstream{particles});
    if (decaymodes)
      configuration["decaymodes"] = read_all(bf::ifstream{decaymodes});
    if (modus)
      configuration["General"]["Modus"] = std::string(modus);
    if (end_time)
      configuration["General"]["End_Time"] = std::abs(atof(end_time));

    /* set up logging */
    set_default_loglevel(configuration.take({"Logging", "default"},
                                            einhard::ALL));
    create_all_loggers(configuration["Logging"]);
    log.info(progname, " (", VERSION_MAJOR, ')');

    int64_t seed = configuration.read({"General",  "Randomseed"});
    if (seed < 0) {
      // Seed with a real random value, if available
      std::random_device rd;
      seed = rd();
      configuration["General"]["Randomseed"] = seed;
    }

    /* check output path*/
    ensure_path_is_valid(output_path);
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
    seed = configuration.take({"General",  "Randomseed"});
    Random::set_seed(seed);
    log.info() << "Random number seed: " << seed;

    log.trace(source_location, " create ParticleType and DecayModes");
    ParticleType::create_type_list(configuration.take({"particles"}));
    DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
    ParticleType::check_consistency();

    // create outputs
    log.trace(source_location, " create OutputInterface objects");
    OutputsList output_list;
    /*!\Userguide
     * \page input_output_options_ Output
     *
     * \key Output: \n
     * Below this key the configuration for the different output formats is
     * defined. To enable a certain output, set the 'Enable' key below the
     * selected format identifier. The identifiers are described below.
     * The following outputs exist:
     * \li \subpage input_oscar_particlelist
     * \li \subpage input_oscar_collisions
     * \li \subpage input_vtk
     * \li \subpage input_binary_collisions
     * \li \subpage input_binary_particles
     * \li \subpage input_root
     */
    auto output_conf = configuration["Output"];

    /*!\Userguide
     * \page output_general_ Output files
     * There are different optional formats for SMASH output that are explained
     * below in more detail. Per default, the selected output files will be
     * saved in the directory ./data/\<run_id\>, where \<run_id\> is an integer
     * number starting from 0. At the beginning
     * of a run SMASH checks, if the ./data/0 directory exists. If it does not exist, it
     * is created and all output files are written there. If the directory
     * already exists, SMASH tries for ./data/1, ./data/2 and so on until it
     * finds a free number. The user can change output directory by a command
     * line option, if desired:
     * \code smash -o <user_output_dir> \endcode
     * SMASH supports several kinds of configurable output formats.
     * They are called OSCAR1999, OSCAR2013, binary OSCAR2013, VTK and ROOT
     * outputs. Every format can be switched on/off using option Enable in the
     * configuration file config.yaml. For more information on configuring the
     * output see corresponding pages: \ref input_oscar_particlelist,
     * \ref input_oscar_collisions, \ref input_binary_collisions,
     * \ref input_binary_particles, \ref input_root, \ref input_vtk.
     *
     * \key Details of output formats are explained here: \n
     * \li General block structure of OSCAR formats: \n
     *     \subpage oscar_general_
     * \li A family of OSCAR ASCII outputs.\n
     *     \subpage format_oscar_particlelist\n
     *     \subpage format_oscar_collisions
     * \li Binary outputs analoguous to OSCAR format\n
     *     \subpage format_binary_\n
     * \li Output in vtk format suitable for an easy
     *     visualization using paraview software:\n \subpage format_vtk
     * \li Formatted binary output that uses ROOT software
     *     (http://root.cern.ch).\n Fast to read and write, requires less
     *     disk space.\n \subpage format_root
     */

    // loop until all OSCAR outputs are created (create_oscar_output will return
    // nullptr then).
    while (std::unique_ptr<OutputInterface> oscar =
               create_oscar_output(output_path, output_conf)) {
      output_list.emplace_back(std::move(oscar));
    }
    if (static_cast<bool>(output_conf.take({"Vtk", "Enable"}))) {
      output_list.emplace_back(new VtkOutput(output_path, output_conf["Vtk"]));
    } else {
      output_conf.take({"Vtk"});
    }
    if (static_cast<bool>(output_conf.take({"Binary_Collisions", "Enable"}))) {
      output_list.emplace_back(new BinaryOutputCollisions(output_path,
                                       output_conf["Binary_Collisions"]));
    } else {
      output_conf.take({"Binary_Collisions"});
    }
    if (static_cast<bool>(output_conf.take({"Binary_Particles", "Enable"}))) {
      output_list.emplace_back(new BinaryOutputParticles(output_path,
                                       output_conf["Binary_Particles"]));
    } else {
      output_conf.take({"Binary_Particles"});
    }
    if (static_cast<bool>(output_conf.take({"Root", "Enable"}))) {
#ifdef SMASH_USE_ROOT
      output_list.emplace_back(new RootOutput(
                               output_path, output_conf["Root"]));
#else
      log.error() << "You requested Root output, but Root support has not been "
                     "compiled in. To enable Root support call: cmake -D "
                     "USE_ROOT=ON <path>.";
      output_conf.take({"Root"});
#endif
    } else {
      output_conf.take({"Root"});
    }
    if (static_cast<bool>(output_conf.take({"Density", "Enable"}))) {
      output_list.emplace_back(new DensityOutput(output_path,
                               output_conf["Density"]));
    } else {
      output_conf.take({"Density"});
    }

    // create an experiment
    log.trace(source_location, " create Experiment");
    auto experiment = ExperimentBase::create(configuration);
    const std::string report = configuration.unused_values_report();
    if (report != "{}") {
      log.warn() << "The following configuration values were not used:\n"
                 << report;
    }

    // set outputs to experiment
    experiment->set_outputs(std::move(output_list));

    // run the experiment
    log.trace(source_location, " run the Experiment");
    experiment->run();
  }
  catch(std::exception &e) {
    log.fatal() << "SMASH failed with the following error:\n" << e.what();
    return EXIT_FAILURE;
  }

  DecayModes::clear_decaymodes();

  log.trace() << source_location << " about to return from main";
  return 0;
}

