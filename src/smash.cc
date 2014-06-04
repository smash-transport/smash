/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <list>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "include/configuration.h"
#include "include/experiment.h"
#include "include/forwarddeclarations.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particletype.h"

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
  printf("\nUsage: %s [option]\n\n", progname.c_str());
  printf("Calculate transport box\n"
         "  -h, --help              usage information\n"
         "\n"
         "  -i, --inputfile <file>  path to input configuration file\n"
         "  -d, --decaymodes <file> override default decay modes from file\n"
         "  -p, --particles <file>  override default particles from file\n"
         "\n"
         "  -c, --config <YAML>     specify config value overrides\n"
         "  -m, --modus <modus>     shortcut for -c 'General: { MODUS: <modus> }'\n"
         "  -e, --endtime <time>    shortcut for -c 'General: { END_TIME: <time> }'\n"
         "\n"
         "  -o, --output <dir>      output directory (default: $PWD/data/<runid>)\n"
         "  -f, --force             force overwriting files in the output directory\n"
         "  -v, --version\n\n");
  exit(rc);
}

/** Exception class that is thrown if the requested output directory
 * already exists and `-f` was not specified on the command line.
 */
struct OutputDirectoryExists : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
/** Exception class that is thrown if no new output path can be
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

/**
 * Utility function to read a complete input stream (e.g. file) into one string.
 *
 * \param input The input stream. Since it reads until EOF und thus "uses up the
 * whole input stream" the function takes an rvalue reference to the stream
 * object (just pass a temporary).
 *
 * \note There's no slicing here: the actual istream object is a temporary that
 * is not destroyed until read_all returns.
 */
std::string read_all(std::istream &&input) {
  return {std::istreambuf_iterator<char>{input},
          std::istreambuf_iterator<char>{}};
}

}  // unnamed namespace

}  // namespace Smash

/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  using namespace Smash;
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
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  const std::string progname = bf::path(argv[0]).filename().native();
  printf("%s (%d)\n", progname.c_str(), VERSION_MAJOR);

  try {
    bool force_overwrite = false;
    bf::path output_path = default_output_path();

    /* read in config file */
    Configuration configuration(".");

    /* check for overriding command line arguments */
    int opt;
    while ((opt = getopt_long(argc, argv, "c:d:e:fhi:m:p:o:v", longopts,
                              nullptr)) != -1) {
      switch (opt) {
        case 'c':
          configuration.merge_yaml(optarg);
          break;
        case 'd': {
          configuration["decaymodes"] = read_all(bf::ifstream{optarg});
        } break;
        case 'f':
          force_overwrite = true;
          break;
        case 'i': {
          const bf::path file(optarg);
          configuration = Configuration(file.parent_path(), file.filename());
        } break;
        case 'h':
          usage(EXIT_SUCCESS, progname);
          break;
        case 'm':
          configuration["General"]["MODUS"] = std::string(optarg);
          break;
        case 'p': {
          configuration["particles"] = read_all(bf::ifstream{optarg});
        } break;
        case 'e':
          configuration["General"]["END_TIME"] = abs(atof(optarg));
          break;
        case 'o':
          output_path = optarg;
          break;
        case 'v':
          exit(EXIT_SUCCESS);
        default:
          usage(EXIT_FAILURE, progname);
      }
    }

    ensure_path_is_valid(output_path);
    printd("output path: %s\n", output_path.native().c_str());

    if (!force_overwrite && bf::exists(output_path / "config.yaml")) {
      throw std::runtime_error(
          "Output directory would get overwritten. Select a different output "
          "directory, clean up, or tell SMASH to ignore existing files.");
    }

    // keep a copy of the configuration that was used in the output directory
    bf::ofstream(output_path / "config.yaml") << configuration.to_string()
                                              << '\n';

    ParticleType::create_type_list(configuration.take({"particles"}));
    auto experiment = ExperimentBase::create(configuration);
    const std::string report = configuration.unused_values_report();
    if (report != "{}") {
      printf("The following configuration values were not used:\n%s\n",
             report.c_str());
    }
    experiment->run(output_path);
  }
  catch(std::exception &e) {
    printf("SMASH failed with the following error:\n%s\n", e.what());
    return EXIT_FAILURE;
  }

  return 0;
}
