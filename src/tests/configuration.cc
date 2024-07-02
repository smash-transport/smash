/*
 *
 *    Copyright (c) 2014-2015,2017-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/configuration.h"

#include <filesystem>

#include "setup.h"
#include "smash/forwarddeclarations.h"
#include "smash/input_keys.h"
#include "smash/macros.h"

using namespace smash;

static Configuration make_test_configuration() {
  return Configuration{
      std::filesystem::path{TEST_CONFIG_PATH} / "src" / "tests",
      "test_config.yaml"};
}

template <typename T>
static Key<T> get_key(KeyLabels labels) {
  return Key<T>{labels, {"1.0"}};
}

TEST(create_object) {
  Configuration conf = make_test_configuration();
  conf.clear();
}

TEST(merge_does_override) {
  Configuration conf = make_test_configuration();
  const auto key_1 = get_key<int>({"fireballs", "arena"});
  const auto key_2 = get_key<int>({"fireballs", "classify"});
  COMPARE(conf.read(key_1), 1000);
  COMPARE(conf.read(key_2), 1);
  conf.merge_yaml("fireballs: { classify: 2 }");
  COMPARE(conf.read(key_1), 1000);
  COMPARE(conf.read(key_2), 2);
  conf.clear();
}

TEST_CATCH(merge_with_incorrect_indent, Configuration::ParseError) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("fireballs:\n foo: 1\n  test: 1\n");
}

TEST(take) {
  Configuration conf = make_test_configuration();
  const double d = conf.take(get_key<double>({"tamer", "pipit", "bushelling"}));
  COMPARE(d, 5.);
  conf.clear();
}

TEST(take_multiple) {
  Configuration conf = make_test_configuration();
  double d = conf.take(get_key<double>({"tamer", "Altaic", "Meccas"}));
  COMPARE(d, 10.);
  d = conf.take(get_key<double>({"tamer", "Altaic", "Kathleen"}));
  COMPARE(d, 0.2);
  const int i = conf.take(get_key<int>({"tamer", "Altaic", "Brahmins"}));
  COMPARE(i, 1);
  conf.clear();
}

TEST_CATCH(take_twice_optional_key, Configuration::TakeSameKeyTwice) {
  Configuration conf = make_test_configuration();
  const Key<double> optional_key{
      {"tamer", "pipit", "bushelling"}, 3.14, {"1.0"}};
  [[maybe_unused]] double d = conf.take(optional_key);
  d = conf.take(optional_key);
}

TEST_CATCH(take_incorrect_type, Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  const Key<int> key = get_key<int>({"tamer", "pipit", "bushelling"});
  [[maybe_unused]] int i = conf.take(key);
}

TEST(take_cannot_always_convert_to_string) {
  Configuration conf = make_test_configuration();
  auto key = get_key<double>({"tamer", "pipit", "bushelling"});
  const bool b = std::is_convertible_v<decltype(conf.take(key)), std::string>;
  VERIFY(!b);
  conf.clear();
}

TEST_CATCH(take_not_existing_key, std::invalid_argument) {
  Configuration conf = make_test_configuration();
  conf.take(get_key<int>({"not existing key"}));
}

TEST_CATCH(take_not_existing_key_in_existing_section, std::invalid_argument) {
  Configuration conf = make_test_configuration();
  conf.take(get_key<int>({"tamer", "not existing key"}));
}

TEST(take_array) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123, 456, 789]}");
  std::array<int, 3> x = conf.take(get_key<std::array<int, 3>>({"test"}));
  VERIFY(x[0] == 123 && x[1] == 456 && x[2] == 789);
  conf.clear();
}

TEST(take_enum_key) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: on}");
  FermiMotion x = conf.take(get_key<FermiMotion>({"test"}));
  VERIFY(x == FermiMotion::On);
  conf.clear();
}

TEST_CATCH(take_array_wrong_n, Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123, 456, 789]}");
  auto key = get_key<std::array<int, 4>>({"test"});
  conf.take(key);
}

TEST(take_reactions_bitset) {
  // Make sure that only the right bits are set
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [NN_to_NR, KN_to_KN]}");
  ReactionsBitSet bs = conf.take(get_key<ReactionsBitSet>({"test"}));
  for (std::size_t i = 0; i < bs.size(); i++) {
    if (i == IncludedReactions::NN_to_NR || i == IncludedReactions::KN_to_KN) {
      VERIFY(bs.test(i));
    } else {
      VERIFY(!bs.test(i));
    }
  }
  // Make sure that all bits are set
  conf.merge_yaml("{test2: [All]}");
  ReactionsBitSet bs2 = conf.take(get_key<ReactionsBitSet>({"test2"}));
  for (std::size_t i = 0; i < bs2.size(); i++) {
    VERIFY(bs2.test(i));
  }
  // All means really ALL reactions are on
  conf.merge_yaml("{test3: [NN_to_NR, All]}");
  ReactionsBitSet bs3 = conf.take(get_key<ReactionsBitSet>({"test3"}));
  for (std::size_t i = 0; i < bs3.size(); i++) {
    VERIFY(bs3.test(i));
  }
  conf.clear();
}

TEST(take_removes_entry) {
  Configuration conf = make_test_configuration();
  const auto key = get_key<double>({"tamer", "pipit", "bushelling"});
  VERIFY(conf.has_value(key));
  conf.take(key);
  VERIFY(!conf.has_value(key));
  conf.clear();
}

TEST(take_removes_empty_section) {
  Configuration conf{R"(
    Section:
      Sub-section:
        Key: "Value"
  )"};
  const auto section = get_key<std::string>({"Section"}),
             subsection = get_key<std::string>({"Section", "Sub-section"}),
             key = get_key<std::string>({"Section", "Sub-section", "Key"});
  VERIFY(conf.has_value(section));
  VERIFY(conf.has_value(subsection));
  conf.take(key);
  VERIFY(!conf.has_value(section));
  VERIFY(!conf.has_value(subsection));
}

TEST(take_removes_empty_section_but_not_empty_lists) {
  Configuration conf{R"(
    Section:
      Sub-section:
        Key: "Value"
        Empty_list: []
  )"};
  conf.take(get_key<std::string>({"Section", "Sub-section", "Key"}));
  VERIFY(conf.has_value(get_key<std::string>({"Section", "Sub-section"})));
  conf.clear();
}

TEST_CATCH(read_failed_sequence_conversion,
           Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123 456]}");
  [[maybe_unused]] auto x = conf.read(get_key<std::vector<int>>({"test"}));
}

TEST(read_check_config_general_contents) {
  Configuration conf = make_test_configuration();

  auto modus = conf.read(get_key<std::string>({"fireballs", "extorting"}));
  COMPARE(modus, "feathered");
  COMPARE(conf.read(get_key<double>({"fireballs", "infection"})), 0.01);
  COMPARE(conf.read(get_key<int>({"fireballs", "arena"})), 1000);
  COMPARE(conf.read(get_key<int>({"fireballs", "pendulous"})), 10);
  COMPARE(conf.read(get_key<int>({"fireballs", "scudded"})), 1);
  COMPARE(conf.read(get_key<double>({"fireballs", "firebrands"})), 10.0);
  COMPARE(conf.read(get_key<int>({"fireballs", "joker"})), 1);
  COMPARE(conf.read(get_key<int>({"fireballs", "classify"})), 1);
  conf.clear();
}

TEST(read_check_config_collider_contents) {
  Configuration conf = make_test_configuration();
  COMPARE(conf.read(get_key<int>({"tamer", "schmoozed", "warbler"})), 211);
  COMPARE(conf.read(get_key<int>({"tamer", "schmoozed", "neglects"})), -211);
  COMPARE(conf.read(get_key<double>({"tamer", "schmoozed", "reedier"})), 1.0);
  conf.clear();
}

TEST(read_does_not_take_and_reading_multiple_times_is_fine) {
  Configuration conf = make_test_configuration();
  const auto key = get_key<int>({"fireballs", "classify"});
  auto nevents = conf.read(key);
  COMPARE(nevents, 1);
  nevents = conf.read(key);
  COMPARE(nevents, 1);
  nevents = conf.read(key);
  COMPARE(nevents, 1);
  conf.clear();
}

TEST(set_existing_value) {
  Configuration conf = make_test_configuration();
  const double new_value = 3.1415;
  const auto key = get_key<double>({"tamer", "Altaic", "Meccas"});
  conf.set_value(key, new_value);
  COMPARE(conf.read(key), new_value);
  conf.clear();
}

TEST(set_new_value_on_non_empty_conf) {
  Configuration conf = make_test_configuration();
  const auto key = get_key<double>({"Test"});
  VERIFY(!conf.has_value(key));
  conf.set_value(key, 1.);
  VERIFY(conf.has_value(key));
  COMPARE(conf.read(key), 1.);
  conf.clear();
}

TEST(set_value_on_empty_conf) {
  auto conf = Configuration("");
  const auto key = get_key<int>({"New section", "New key"});
  conf.set_value(key, 42);
  VERIFY(conf.has_value(get_key<std::string>({"New section"})));
  VERIFY(conf.has_value(key));
  conf.clear();
}

TEST(set_value_on_conf_created_with_empty_file) {
  auto tmp_dir = std::filesystem::temp_directory_path();
  auto tmp_file = "empty_config.yaml";
  std::ofstream ofs(tmp_dir / tmp_file);
  if (!ofs) {
    FAIL() << "Unable to create empty temporary file!";
  }
  auto conf = Configuration{tmp_dir, tmp_file};
  const auto key = get_key<int>({"New section", "New key"});
  conf.set_value(key, 42);
  VERIFY(conf.has_value(get_key<std::string>({"New section"})));
  VERIFY(conf.has_value(key));
  ofs.close();
  std::filesystem::remove(tmp_dir / tmp_file);
  conf.clear();
}

TEST(remove_all_entries_in_section_but_one) {
  Configuration conf = make_test_configuration();
  conf.remove_all_entries_in_section_but_one("pipit", {"tamer"});
  conf.remove_all_entries_in_section_but_one("tamer", {});
  COMPARE(conf.to_string(), "tamer:\n  pipit:\n    bushelling: 5.0");
  conf.clear();
}

TEST(extract_sub_configuration) {
  Configuration conf = make_test_configuration();
  Configuration sub_conf = conf.extract_sub_configuration({"tamer", "pipit"});
  VERIFY(!conf.has_value(get_key<std::string>({"tamer", "pipit"})));
  COMPARE(sub_conf.to_string(), "bushelling: 5.0");
  sub_conf = conf.extract_sub_configuration({"fireballs"});
  const auto list_of_keys = sub_conf.list_upmost_nodes();
  const std::vector<std::string> reference = {
      "extorting", "infection",  "arena", "pendulous",
      "scudded",   "firebrands", "joker", "classify"};
  COMPARE(list_of_keys.size(), reference.size());
  for (std::size_t i = 0; i < list_of_keys.size(); ++i) {
    COMPARE(list_of_keys[i], reference[i]);
  }
  conf.clear();
  sub_conf.clear();
}

TEST_CATCH(extract_scalar_key_as_section, std::runtime_error) {
  Configuration conf = make_test_configuration();
  auto sub_conf = conf.extract_sub_configuration({"fireballs", "joker"});
}

TEST_CATCH(extract_sequence_key_as_section, std::runtime_error) {
  Configuration conf = make_test_configuration();
  auto sub_conf =
      conf.extract_sub_configuration({"tamer", "feathered", "stopcock"});
}

TEST_CATCH(extract_empty_map_as_section, std::runtime_error) {
  Configuration conf{"section: {}"};
  auto sub_conf = conf.extract_sub_configuration({"section"});
}

TEST_CATCH(extract_key_without_value_section, std::runtime_error) {
  Configuration conf{"section:"};
  auto sub_conf = conf.extract_sub_configuration({"section"});
}

TEST_CATCH(extract_not_existing_section, std::runtime_error) {
  Configuration conf = make_test_configuration();
  auto sub_conf = conf.extract_sub_configuration({"Not existing section"});
}

TEST(extract_not_existing_section_as_empty_conf) {
  Configuration conf = make_test_configuration();
  auto sub_conf = conf.extract_sub_configuration({"Not existing section"},
                                                 Configuration::GetEmpty::Yes);
  COMPARE(sub_conf.to_string(), "");
  conf.clear();
}

TEST(extract_complete_sub_configuration) {
  Configuration conf = make_test_configuration();
  KeyLabels section{"tamer", "feathered", "dove", "floozy"};
  Configuration sub_conf = conf.extract_complete_sub_configuration(section);
  conf.clear();
  const std::string result_as_string(R"(tamer:
  feathered:
    dove:
      floozy: {2212: 1, 2112: 1})");
  COMPARE(sub_conf.to_string(), result_as_string);
  for ([[maybe_unused]] const auto &label : section) {
    auto key = get_key<int>(
        {section.begin(), section.end()});  // Type of key here irrelevant
    VERIFY(sub_conf.has_value(key));
    section.pop_back();
  }
  sub_conf.clear();
}

TEST(extract_complete_not_existing_section_as_empty_conf) {
  Configuration conf = make_test_configuration();
  KeyLabels section{"Not", "existing", "section"};
  auto sub_conf = conf.extract_complete_sub_configuration(
      section, Configuration::GetEmpty::Yes);
  conf.clear();
  const std::string result_as_string(R"(Not:
  existing:
    section:
      {})");
  COMPARE(sub_conf.to_string(), result_as_string);
  for ([[maybe_unused]] const auto &label : section) {
    auto key = get_key<int>(
        {section.begin(), section.end()});  // Type of key here irrelevant
    VERIFY(sub_conf.has_value(key));
    section.pop_back();
  }
  sub_conf.clear();
}

TEST(has_value_including_empty) {
  Configuration conf = Configuration{"Empty:"};
  const auto key = get_key<std::string>({"Empty"});
  VERIFY(!conf.has_value(key));
  VERIFY(conf.has_key(key));
  conf.clear();
}

TEST(has_value) {
  Configuration conf = make_test_configuration();
  VERIFY(conf.has_value(get_key<double>({"tamer", "pipit", "bushelling"})));
  conf.clear();
}

TEST(is_empty) {
  Configuration conf{""};
  VERIFY(conf.is_empty());
  conf = Configuration{"Key: Value"};
  VERIFY(!conf.is_empty());
  conf.clear();
}

TEST(to_string) {
  Configuration conf = Configuration{""};
  COMPARE(conf.to_string(), "{}");
}

TEST(validate) {
  // Disable logger output -> reenable if needed to e.g. debug
  logg[LogArea::Configuration::id].setVerbosity(einhard::OFF);
  Configuration invalid_conf = make_test_configuration();
  VERIFY(invalid_conf.validate(false) == Configuration::Is::Invalid);
  VERIFY(invalid_conf.validate(true) == Configuration::Is::Invalid);
  Configuration deprecated_conf = Configuration{"Version: 1.8"};
  VERIFY(deprecated_conf.validate() == Configuration::Is::Deprecated);
  invalid_conf.clear();
  deprecated_conf.clear();
  // Reenable logger output (it is global)
  logg[LogArea::Configuration::id].setVerbosity(einhard::ALL);
}

TEST(validate_shipped_input_files) {
  const std::filesystem::path codebase_path{TEST_CONFIG_PATH};
  const std::string input_folder_name{"input"};
  const std::string extension(".yaml");
  for (auto &input_file : std::filesystem::recursive_directory_iterator(
           codebase_path / input_folder_name)) {
    if (input_file.path().extension() == extension) {
      /* Use 0 just because the logging area is irrelevant here and
         in any case not set, since create_all_loggers is not called! */
      logg[0].debug() << " Validating " << input_file.path() << '\n';
      Configuration config{input_file.path().parent_path(),
                           input_file.path().filename()};
      VERIFY(config.validate(false) == Configuration::Is::Valid);
      VERIFY(config.validate(true) == Configuration::Is::Valid);
      config.clear();
    }
  }
}

// Sorry, but I have to put this in the std namespace, otherwise it doesn't
// compile. That's because the << operator is called from inside the vir::test
// namespace and all involved types are in the std namespace.
namespace std {
static ostream &operator<<(ostream &s, const vector<string> &v) {
  s << '{';
  for (const auto &x : v) {
    s << x << ", ";  // I'm too lazy to get the commas right
  }
  return s << '}';
}
}  // namespace std

static void expect_lines(std::vector<std::string> expected,
                         std::istream &stream) {
  std::string line;
  while (!expected.empty()) {
    getline(stream, line);
    const auto pos = find(expected.begin(), expected.end(), line);
    VERIFY(pos != expected.end()) << line << " was not in " << expected;
    expected.erase(pos);
  }
}

// Test not on Configuration functionality but more as integration test
TEST(check_unused_report) {
  std::string reference;
  Configuration conf = make_test_configuration();
  conf.take({"fireballs", "extorting"});
  conf.take({"fireballs", "infection"});
  conf.take({"fireballs", "arena"});
  conf.take({"fireballs", "pendulous"});
  conf.take({"fireballs", "scudded"});
  conf.take({"fireballs", "firebrands"});
  conf.take({"fireballs", "joker"});
  conf.take({"fireballs", "classify"});
  conf.take({"tamer", "Altaic", "Meccas"});
  conf.take({"tamer", "Altaic", "Kathleen"});
  conf.take({"tamer", "Altaic", "Brahmins"});
  conf.take({"tamer", "feathered"});
  {
    std::istringstream unused(conf.to_string());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    if (line == "  schmoozed:") {
      expect_lines(
          {"    neglects: -211", "    warbler: 211", "    reedier: 1.0"},
          unused);
      getline(unused, line);
      COMPARE(line, "  pipit:");
      getline(unused, line);
      COMPARE(line, "    bushelling: 5.0");
    } else {
      COMPARE(line, "  pipit:");
      getline(unused, line);
      COMPARE(line, "    bushelling: 5.0");
      getline(unused, line);
      COMPARE(line, "  schmoozed:");
      expect_lines(
          {"    neglects: -211", "    warbler: 211", "    reedier: 1.0"},
          unused);
    }
    VERIFY(unused.eof());
  }

  conf.take({"tamer", "pipit", "bushelling"});
  {
    std::istringstream unused(conf.to_string());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    expect_lines({"    neglects: -211", "    warbler: 211", "    reedier: 1.0"},
                 unused);
    VERIFY(unused.eof());
  }

  conf.take({"tamer", "schmoozed", "warbler"});
  {
    std::istringstream unused(conf.to_string());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    expect_lines({"    neglects: -211", "    reedier: 1.0"}, unused);
    VERIFY(unused.eof());
  }

  conf.take({"tamer", "schmoozed", "reedier"});
  {
    std::istringstream unused(conf.to_string());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    getline(unused, line);
    COMPARE(line, "    neglects: -211");
    VERIFY(unused.eof());
  }

  conf.take({"tamer", "schmoozed", "neglects"});
  reference = "{}";
  COMPARE(conf.to_string(), reference);
}
