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
  OutputFormatter<ToASCII> formatter(empty);
}

TEST_CATCH(invalid_quantity, std::invalid_argument) {
  std::vector<std::string> invalid_quantities = {"gibberish"};
  OutputFormatter<ToASCII> formatter(invalid_quantities);
}

TEST_CATCH(repeated_quantity, std::invalid_argument) {
  std::vector<std::string> repeated_quantities = {"t", "t"};
  OutputFormatter<ToASCII> formatter(repeated_quantities);
}

TEST_CATCH(incompatible_quantity, std::invalid_argument) {
  std::vector<std::string> repeated_quantities = {"id", "ID"};
  OutputFormatter<ToASCII> formatter(repeated_quantities);
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

  OutputFormatter<ToASCII> formatter(valid_quantities);

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

  VERIFY(correct_line.str() == formatter.particle_line(p));
}

TEST(binary_particle_line) {
  ParticleData p = Test::smashon_random();
  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<ToBinary> formatter(quantities);

  std::vector<char> chunk = formatter.particle_line(p);

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

TEST(per_particle_particle_line_same_as_particles_chunk) {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);

  for (std::size_t i = 0; i < N; ++i) {
    particles.push_back(Test::smashon_random());
  }

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<ToBinary> formatter_binary(quantities);
  OutputFormatter<ToASCII> formatter_ascii(quantities);

  // ---- Binary ----
  const std::vector<char> particles_chunk_binary =
      formatter_binary.particles_chunk(particles);

  std::vector<char> per_particle_chunk_binary;
  per_particle_chunk_binary.reserve(particles_chunk_binary.size());
  for (const auto& p : particles) {
    const auto data = formatter_binary.particle_line(p);
    per_particle_chunk_binary.insert(per_particle_chunk_binary.end(),
                                     data.begin(), data.end());
  }

  VERIFY(particles_chunk_binary == per_particle_chunk_binary);

  // ---- ASCII ----
  const std::string particles_chunk_ascii =
      formatter_ascii.particles_chunk(particles);

  std::string per_particle_chunk_ascii;
  per_particle_chunk_ascii.reserve(particles_chunk_ascii.size());
  for (const auto& p : particles) {
    const auto line = formatter_ascii.particle_line(p);
    per_particle_chunk_ascii.append(line);
  }

  VERIFY(particles_chunk_ascii == per_particle_chunk_ascii);
}

TEST(compute_single_size_equals_particle_line_size) {
  ParticleData p = Test::smashon_random();
  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};

  OutputFormatter<ToBinary> formatter_binary(quantities);
  const auto line_binary = formatter_binary.particle_line(p);
  const std::size_t computed_binary = formatter_binary.compute_single_size(p);
  VERIFY(computed_binary == line_binary.size());

  OutputFormatter<ToASCII> formatter_ascii(quantities);
  const auto line_ascii = formatter_ascii.particle_line(p);
  const std::size_t computed_ascii = formatter_ascii.compute_single_size(p);
  VERIFY(computed_ascii == line_ascii.size());
}

TEST(sum_of_single_sizes_equals_particles_chunk_size) {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    particles.push_back(Test::smashon_random());
  }

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};

  OutputFormatter<ToBinary> formatter_binary(quantities);
  const auto chunk_binary = formatter_binary.particles_chunk(particles);
  std::size_t sum_binary = 0;
  for (const auto& p : particles)
    sum_binary += formatter_binary.compute_single_size(p);
  VERIFY(chunk_binary.size() == sum_binary);

  OutputFormatter<ToASCII> formatter_ascii(quantities);
  const auto chunk_ascii = formatter_ascii.particles_chunk(particles);
  std::size_t sum_ascii = 0;
  for (const auto& p : particles)
    sum_ascii += formatter_ascii.compute_single_size(p);
  VERIFY(chunk_ascii.size() == sum_ascii);

  const ParticleData sample = Test::smashon_random();
  VERIFY(formatter_ascii.compute_single_size(sample) ==
         formatter_ascii.particle_line(sample).size());
}
