/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/outputformatter.h"

#include "setup.h"

using namespace smash;

TEST(ASCII_converter) {
  ToASCII converter;
  // int goes to an string containing int
  VERIFY(converter.as_integer(42) == std::string{"42"});
  // float precision behaves correctly
  const double pi = 3.1415926535897932384;
  VERIFY(converter.as_double(pi) == std::string{"3.14159"});
  VERIFY(converter.as_precise_double(pi) == std::string{"3.14159265"});
  // strings (literal and not) behave correctly
  const std::string smash_str{"smash"};
  VERIFY(converter.as_string(smash_str) == smash_str);
  VERIFY(converter.as_string("smash") == smash_str);
}

TEST_CATCH(empty_quantities, std::invalid_argument) {
  std::vector<std::string> empty{};
  OutputFormatterASCII formatter(empty);
}

TEST_CATCH(invalid_quantity, std::invalid_argument) {
  std::vector<std::string> invalid_quantities = {"gibberish"};
  OutputFormatterASCII formatter(invalid_quantities);
}

TEST_CATCH(repeated_quantity, std::invalid_argument) {
  std::vector<std::string> repeated_quantities = {"t", "t"};
  OutputFormatterASCII formatter(repeated_quantities);
}

TEST_CATCH(incompatible_quantity, std::invalid_argument) {
  std::vector<std::string> repeated_quantities = {"id", "ID"};
  OutputFormatterASCII formatter(repeated_quantities);
}

TEST(valid_line_maker) {
  Test::create_smashon_particletypes();
  ParticleData p = Test::smashon_random();

  std::vector<std::string> valid_quantities = {"t",
                                               "x",
                                               "y",
                                               "z",
                                               "mass",
                                               "p0",
                                               "px",
                                               "py",
                                               "pz",
                                               "pdg",
                                               "ID",
                                               "charge",
                                               "ncoll",
                                               "form_time",
                                               "xsecfac",
                                               "proc_id_origin",
                                               "proc_type_origin",
                                               "time_last_coll",
                                               "pdg_mother1",
                                               "pdg_mother2",
                                               "baryon_number",
                                               "strangeness"};

  OutputFormatterASCII formatter(valid_quantities);

  const auto quantities_line = std::accumulate(
      std::begin(valid_quantities), std::end(valid_quantities), std::string{},
      [](const std::string& ss, const std::string& s) {
        return ss.empty() ? s : ss + " " + s;
      });
  VERIFY(quantities_line == formatter.quantities_line());

  std::string units_line{
      "fm fm fm fm GeV GeV GeV GeV GeV none none e none fm none none none fm "
      "none none none none"};
  VERIFY(units_line == formatter.unit_line());

  std::stringstream correct_line{};
  for (int i = 0; i < 4; ++i) {
    correct_line << p.position()[i] << " ";
  }
  correct_line << p.effective_mass() << " ";
  // For momentum only, the default precision is 9
  correct_line.precision(9);
  for (int i = 0; i < 4; ++i) {
    correct_line << p.momentum()[i] << " ";
  }
  correct_line.precision(6);
  correct_line << p.pdgcode().string() << " ";
  correct_line << p.id() << " ";
  correct_line << p.type().charge() << " ";
  correct_line << p.get_history().collisions_per_particle << " ";
  correct_line << p.formation_time() << " ";
  correct_line << p.xsec_scaling_factor() << " ";
  correct_line << p.get_history().id_process << " ";
  correct_line << static_cast<int>(p.get_history().process_type) << " ";
  correct_line << p.get_history().time_last_collision << " ";
  correct_line << p.get_history().p1.string() << " ";
  correct_line << p.get_history().p2.string() << " ";
  correct_line << p.pdgcode().baryon_number() << " ";
  correct_line << p.pdgcode().strangeness() << "\n";
  VERIFY(correct_line.str() == formatter.data_line(p));
}

TEST(binary_data_line) {
  ParticleData p = Test::smashon_random();
  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatterBinary formatter(quantities);

  std::vector<char> chunk = formatter.data_line(p);

  double t = *reinterpret_cast<double*>(chunk.data());
  double x = *reinterpret_cast<double*>(chunk.data() + sizeof(double));
  double y = *reinterpret_cast<double*>(chunk.data() + 2 * sizeof(double));
  double z = *reinterpret_cast<double*>(chunk.data() + 3 * sizeof(double));
  int32_t id = *reinterpret_cast<int32_t*>(chunk.data() + 4 * sizeof(double));

  VERIFY(t == p.position()[0]);
  VERIFY(x == p.position()[1]);
  VERIFY(y == p.position()[2]);
  VERIFY(z == p.position()[3]);
  VERIFY(id == p.id());
}

TEST(per_particle_binary_data_line_same_as_particles_chunk) {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);

  for (std::size_t i = 0; i < N; ++i) {
    particles.push_back(Test::smashon_random());  // no default-ctor needed
  }

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatterBinary formatter(quantities);

  // Chunk for whole container
  std::vector<char> particles_chunk = formatter.particles_chunk(particles);

  // Chunk by concatenating per-particle
  std::vector<char> per_particle_chunk;
  per_particle_chunk.reserve(particles_chunk.size());
  for (const auto& p : particles) {
    auto data = formatter.data_line(p);
    per_particle_chunk.insert(per_particle_chunk.end(), data.begin(),
                              data.end());
  }

  VERIFY(particles_chunk == per_particle_chunk);
}
