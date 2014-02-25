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
#include "include/parameters.h"
#include "include/macros.h"
#include "include/outputroutines.h"

/* build dependent variables */
#include "include/config.h"

namespace bf = boost::filesystem;

namespace {
void usage(int rc, const std::string &progname) {
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
         "  -s, --steps <steps>     shortcut for -c 'General: { STEPS: <steps> }'\n"
         "\n"
         "  -o, --output <dir>      output directory (default: $PWD/data/<runid>)\n"
         "  -v, --version\n\n");
  exit(rc);
}

struct OutputDirectoryExists : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
struct OutputDirectoryOutOfIds : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

bf::path default_output_path() {
  const bf::path p = bf::absolute("data");
  bf::create_directory(p);
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
  if (!bf::create_directory(p2)) {
    throw OutputDirectoryExists(
        "Race condition detected: The directory " + p2.native() +
        " did not exist a few cycles ago, but was created in the meantime by a "
        "different process.");
  }
  return p2;
}

void ensure_path_is_valid(const bf::path &path) {
  if (bf::exists(path)) {
    if (!bf::is_directory(path)) {
      throw OutputDirectoryExists("The given path (" + path.native() +
                                  ") exists, but it is not a directory.");
    }
  } else {
    if (!bf::create_directory(path)) {
      throw OutputDirectoryExists(
          "Race condition detected: The directory " + path.native() +
          " did not exist a few cycles ago, but was created in the meantime by "
          "a different process.");
    }
  }
}

std::string read_all(std::istream &input) {
  return {std::istreambuf_iterator<char>{input},
          std::istreambuf_iterator<char>{}};
}

}  // unnamed namespace

/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  constexpr option longopts[] = {
    { "config",     required_argument,      0, 'c' },
    { "decaymodes", required_argument,      0, 'd' },
    { "help",       no_argument,            0, 'h' },
    { "inputfile",  required_argument,      0, 'i' },
    { "modus",      required_argument,      0, 'm' },
    { "particles",  required_argument,      0, 'p' },
    { "steps",      required_argument,      0, 's' },
    { "output",     required_argument,      0, 'o' },
    { "version",    no_argument,            0, 'v' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  const std::string progname = bf::path(argv[0]).filename().native();
  printf("%s (%d)\n", progname.c_str(), VERSION_MAJOR);

  try {
    bf::path output_path = default_output_path();

    /* read in config file */
    Configuration configuration(".");

    /* check for overriding command line arguments */
    int opt;
    while ((opt = getopt_long(argc, argv, "c:d:hi:m:p:s:o:v", longopts,
                              nullptr)) != -1) {
      switch (opt) {
        case 'c':
          configuration.merge_yaml(optarg);
          break;
        case 'd': {
          bf::ifstream file{optarg};
          configuration["decaymodes"] = read_all(file);
        } break;
        case 'i': {
          const boost::filesystem::path file(optarg);
          configuration = Configuration(file.parent_path(), file.filename());
        } break;
        case 'h':
          usage(EXIT_SUCCESS, progname);
          break;
        case 'm':
          configuration["General"]["MODUS"] = std::string(optarg);
          break;
        case 'p': {
          bf::ifstream file{optarg};
          configuration["particles"] = read_all(file);
        } break;
        case 's':
          configuration["General"]["STEPS"] = abs(atoi(optarg));
          break;
        case 'o':
          output_path = optarg;
          ensure_path_is_valid(output_path);
          break;
        case 'v':
          exit(EXIT_SUCCESS);
        default:
          usage(EXIT_FAILURE, progname);
      }
    }

    printd("output path: %s\n", output_path.native().c_str());

    // keep a copy of the configuration that was used in the output directory
    bf::ofstream(output_path / "config.yaml") << configuration.to_string()
                                              << '\n';

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
