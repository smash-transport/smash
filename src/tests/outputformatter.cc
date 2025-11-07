/*
 *
 *    Copyright (c) 2024-2025
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

TEST_CATCH(invalid_quantity, OutputFormatter<ToASCII>::UnknownQuantity) {
  std::vector<std::string> invalid_quantities = {"gibberish"};
  OutputFormatter<ToASCII> formatter(invalid_quantities);
}

TEST_CATCH(repeated_quantity, OutputFormatter<ToASCII>::RepeatedQuantity) {
  std::vector<std::string> repeated_quantities = {"t", "t"};
  OutputFormatter<ToASCII> formatter(repeated_quantities);
}

TEST_CATCH(incompatible_quantity_1, OutputFormatter<ToASCII>::AliasesQuantity) {
  std::vector<std::string> repeated_quantities = {"id", "ID"};
  OutputFormatter<ToASCII> formatter(repeated_quantities);
}

TEST_CATCH(incompatible_quantity_2, OutputFormatter<ToASCII>::AliasesQuantity) {
  std::vector<std::string> repeated_quantities = {"eta", "eta_s"};
  OutputFormatter<ToASCII> formatter(repeated_quantities);
}

TEST_CATCH(incompatible_quantity_3, OutputFormatter<ToASCII>::AliasesQuantity) {
  std::vector<std::string> repeated_quantities = {"Rap", "y_rap"};
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

  auto quantities_line = std::accumulate(
      std::begin(valid_quantities), std::end(valid_quantities), std::string{},
      [](const std::string& ss, const std::string& s) {
        return ss.empty() ? s : ss + " " + s;
      });
  quantities_line.append("\n");
  VERIFY(quantities_line == formatter.quantities_line());

  std::string units_line{
      "fm fm fm fm GeV GeV GeV GeV GeV none none e none fm none none none fm "
      "none none none none\n"};
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

  VERIFY(correct_line.str() == formatter.single_particle_data(p));
}

TEST(binary_single_particle_data) {
  ParticleData p = Test::smashon_random();
  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<ToBinary> formatter(quantities);

  std::vector<char> chunk = formatter.single_particle_data(p);

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

template <typename Converter>
static bool chunks_are_equal() {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    particles.push_back(Test::smashon_random());
  }

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<Converter> formatter(quantities);

  const typename Converter::type full_chunk =
      formatter.particles_data_chunk(particles);

  typename Converter::type concatenated_chunk;
  concatenated_chunk.reserve(full_chunk.size());
  for (const auto& p : particles) {
    const auto line = formatter.single_particle_data(p);
    concatenated_chunk.insert(concatenated_chunk.end(), line.begin(),
                              line.end());
  }

  return full_chunk == concatenated_chunk;
}

TEST(per_particle_single_particle_data_same_as_particles_data_chunk) {
  VERIFY(chunks_are_equal<ToBinary>());
  VERIFY(chunks_are_equal<ToASCII>());
}

template <typename Converter>
static bool single_size_matches_line_size() {
  const ParticleData p = Test::smashon_random();

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<Converter> formatter(quantities);

  const auto line = formatter.single_particle_data(p);
  const std::size_t computed_size = formatter.compute_single_size(p);

  return computed_size == line.size();
}

TEST(compute_single_size_equals_single_particle_data_size) {
  VERIFY(single_size_matches_line_size<ToBinary>());
  VERIFY(single_size_matches_line_size<ToASCII>());
}

template <typename Converter>
static bool sum_of_single_sizes_equals_chunk_size() {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    particles.push_back(Test::smashon_random());
  }

  std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  OutputFormatter<Converter> formatter(quantities);

  const auto chunk = formatter.particles_data_chunk(particles);

  std::size_t total_size = 0;
  for (const auto& p : particles)
    total_size += formatter.compute_single_size(p);

  return chunk.size() == total_size;
}

TEST(sum_of_single_sizes_equals_particles_data_chunk_size) {
  VERIFY(sum_of_single_sizes_equals_chunk_size<ToBinary>());
  VERIFY(sum_of_single_sizes_equals_chunk_size<ToASCII>());
}

template <typename Converter>
static bool write_in_chunks_produces_same_output(
    std::size_t buffer_size, const std::vector<std::string>& quantities) {
  ParticleList particles;
  const std::size_t N = 33;
  particles.reserve(N);
  for (std::size_t i = 0; i < N; ++i)
    particles.push_back(Test::smashon_random());

  OutputFormatter<Converter> formatter(quantities);
  const typename Converter::type one_chunk =
      formatter.particles_data_chunk(particles);

  typename Converter::type multi_chunk;
  details::write_in_chunk_impl<Converter>(
      particles, formatter,
      [&](const typename Converter::type& buf) {
        multi_chunk.insert(multi_chunk.end(), buf.begin(), buf.end());
      },
      buffer_size);

  return one_chunk == multi_chunk;
}

TEST(write_in_chunks_same_as_particles_data_chunk) {
  const std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  constexpr std::size_t bytes_per_particle_est =
      sizeof(double) * 4 + sizeof(int);
  constexpr std::size_t buffer_size = 3 * bytes_per_particle_est;
  VERIFY(
      write_in_chunks_produces_same_output<ToBinary>(buffer_size, quantities));
  VERIFY(
      write_in_chunks_produces_same_output<ToASCII>(buffer_size, quantities));
}

TEST_CATCH(write_in_chunks_throws_if_buffer_too_small, std::runtime_error) {
  const std::vector<std::string> quantities = {"t", "x", "y", "z", "ID"};
  const std::size_t tiny_buffer_size = sizeof(double) * 4 + sizeof(int);
  write_in_chunks_produces_same_output<ToBinary>(tiny_buffer_size, quantities);
  write_in_chunks_produces_same_output<ToASCII>(tiny_buffer_size, quantities);
}
