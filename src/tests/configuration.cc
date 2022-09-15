/*
 *
 *    Copyright (c) 2014-2015,2017-2022
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
#include "smash/macros.h"

using namespace smash;

static Configuration make_test_configuration() {
  return Configuration{
      std::filesystem::path{TEST_CONFIG_PATH} / "src" / "tests",
      "test_config.yaml"};
}

TEST(create_object) { Configuration conf = Test::configuration(); }

TEST(check_config_general_contents) {
  Configuration conf = make_test_configuration();

  std::string modus = conf.read({"fireballs", "extorting"});
  COMPARE(modus, "feathered");
  COMPARE(double(conf.read({"fireballs", "infection"})), 0.01);
  COMPARE(int(conf.read({"fireballs", "arena"})), 1000);
  COMPARE(int(conf.read({"fireballs", "pendulous"})), 10);
  COMPARE(int(conf.read({"fireballs", "scudded"})), 1);
  COMPARE(double(conf.read({"fireballs", "firebrands"})), 10.0);
  COMPARE(int(conf.read({"fireballs", "joker"})), 1);
  COMPARE(int(conf.read({"fireballs", "classify"})), 1);
}

TEST(check_config_collider_contents) {
  Configuration conf = make_test_configuration();
  COMPARE(int(conf.read({"tamer", "schmoozed", "warbler"})), 211);
  COMPARE(int(conf.read({"tamer", "schmoozed", "neglects"})), -211);
  COMPARE(double(conf.read({"tamer", "schmoozed", "reedier"})), 1.0);
}

TEST(test_take) {
  Configuration conf = make_test_configuration();
  double d = conf.take({"tamer", "pipit", "bushelling"});
  COMPARE(d, 5.);
}

TEST(test_take_multiple) {
  Configuration conf = make_test_configuration();
  double d = conf.take({"tamer", "Altaic", "Meccas"});
  COMPARE(d, 10.);
  d = conf.take({"tamer", "Altaic", "Kathleen"});
  COMPARE(d, 0.2);
  int i = conf.take({"tamer", "Altaic", "Brahmins"});
  COMPARE(i, 1);
}

TEST_CATCH(take_incorrect_type, Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  int i = conf.take({"tamer", "pipit", "bushelling"});
  COMPARE(i, 5);
}

TEST(take_always_converts_to_string) {
  Configuration conf = make_test_configuration();
  std::string s = conf.take({"tamer", "pipit", "bushelling"});
  COMPARE(s, "5.0");
}

TEST_CATCH(take_not_existing_key, std::runtime_error) {
  Configuration conf = make_test_configuration();
  conf.take({"not existing key"});
}

TEST_CATCH(take_not_existing_key_in_existing_section, std::runtime_error) {
  Configuration conf = make_test_configuration();
  conf.take({"tamer", "not existing key"});
}

TEST(has_value) {
  Configuration conf = make_test_configuration();
  VERIFY(conf.has_value({"tamer", "pipit", "bushelling"}));
}

TEST(take_removes_entry) {
  Configuration conf = make_test_configuration();
  VERIFY(conf.has_value({"tamer", "pipit", "bushelling"}));
  conf.take({"tamer", "pipit", "bushelling"});
  VERIFY(!conf.has_value({"tamer", "pipit", "bushelling"}));
}

TEST(take_removes_empty_section) {
  Configuration conf{R"(
    Section:
      Sub-section:
        Key: "Value"
  )"};
  VERIFY(conf.has_value({"Section"}));
  VERIFY(conf.has_value({"Section", "Sub-section"}));
  conf.take({"Section", "Sub-section", "Key"});
  VERIFY(!conf.has_value({"Section", "Sub-section"}));
  VERIFY(!conf.has_value({"Section"})) << "\n" << conf.to_string();
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

TEST(test_config_read) {
  Configuration conf = make_test_configuration();
  int nevents = conf.read({"fireballs", "classify"});
  COMPARE(nevents, 1);
  nevents = conf.read({"fireballs", "classify"});
  COMPARE(nevents, 1);
  nevents = conf.take({"fireballs", "classify"});
  COMPARE(nevents, 1);
}

TEST(check_setting_new_value) {
  Configuration conf = make_test_configuration();
  VERIFY(!conf.has_value({"Test"}));
  conf.set_value({"Test"}, 1.);
  VERIFY(conf.has_value({"Test"}));
  COMPARE(double(conf.read({"Test"})), 1.);
}

TEST(check_changing_existing_value) {
  Configuration conf = make_test_configuration();
  const double new_value = 3.1415;
  conf.set_value({"tamer", "Altaic", "Meccas"}, new_value);
  COMPARE(double(conf.read({"tamer", "Altaic", "Meccas"})), new_value);
}

TEST(merge_override) {
  Configuration conf = make_test_configuration();
  COMPARE(int(conf.read({"fireballs", "arena"})), 1000);
  COMPARE(int(conf.read({"fireballs", "classify"})), 1);
  conf.merge_yaml("fireballs: { classify: 2 }");
  COMPARE(int(conf.read({"fireballs", "arena"})), 1000);
  COMPARE(int(conf.read({"fireballs", "classify"})), 2);
}

TEST(remove_all_entries_in_section_but_one) {
  Configuration conf = make_test_configuration();
  conf.remove_all_entries_in_section_but_one("pipit", {"tamer"});
  conf.remove_all_entries_in_section_but_one("tamer", {});
  COMPARE(conf.to_string(), "tamer:\n  pipit:\n    bushelling: 5.0");
}

TEST(extract_sub_configuration) {
  Configuration conf = make_test_configuration();
  Configuration sub_conf = conf.extract_sub_configuration({"tamer", "pipit"});
  VERIFY(!conf.has_value({"tamer", "pipit"}));
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
}

TEST(test_sub_config_objects) {
  Configuration conf = make_test_configuration();
  Configuration general = conf.extract_sub_configuration({"fireballs"});
  const Configuration box = conf.extract_sub_configuration({"tamer", "Altaic"});
  VERIFY(general.has_value({"classify"}));
  int nevents = general.read({"classify"});
  VERIFY(general.has_value({"classify"}));
  COMPARE(nevents, 1);
  nevents = general.take({"classify"});
  VERIFY(!general.has_value({"classify"}));
  COMPARE(nevents, 1);
  COMPARE(double(box.read({"Meccas"})), 10.);
}

TEST_CATCH(failed_sequence_conversion,
           Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123 456]}");
  std::vector<int> x = conf.read({"test"});
}

TEST_CATCH(incorrect_indent, Configuration::ParseError) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("fireballs:\n foo: 1\n  test: 1\n");
  int x = conf.read({"fireballs", "test"});
  COMPARE(x, 1);
}

TEST(take_array) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123, 456, 789]}");
  std::array<int, 3> x = conf.take({"test"});
  VERIFY(x[0] == 123 && x[1] == 456 && x[2] == 789);
}

TEST_CATCH(take_array_wrong_n, Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [123, 456, 789]}");
  std::array<int, 4> x = conf.take({"test"});
  SMASH_UNUSED(x);
}

TEST(reactions_bitset) {
  // Make sure that only the right bits are set
  Configuration conf = make_test_configuration();
  conf.merge_yaml("{test: [NN_to_NR, KN_to_KN]}");
  ReactionsBitSet bs = conf.take({"test"});
  for (std::size_t i = 0; i < bs.size(); i++) {
    if (i == IncludedReactions::NN_to_NR || i == IncludedReactions::KN_to_KN) {
      VERIFY(bs.test(i));
    } else {
      VERIFY(!bs.test(i));
    }
  }
  // Make sure that all bits are set
  conf.merge_yaml("{test2: [All]}");
  ReactionsBitSet bs2 = conf.take({"test2"});
  for (std::size_t i = 0; i < bs2.size(); i++) {
    VERIFY(bs2.test(i));
  }
  // All means really ALL reactions are on
  conf.merge_yaml("{test3: [NN_to_NR, All]}");
  ReactionsBitSet bs3 = conf.take({"test3"});
  for (std::size_t i = 0; i < bs3.size(); i++) {
    VERIFY(bs3.test(i));
  }
}

TEST(configuration_validation) {
  const char *valid_conf_literal = R"(
      General:
          Modus: Box
      )";
  // The valid_conf_literal must be replaced with SMASH config.yaml file
  Configuration valid_conf{valid_conf_literal};
  Configuration invalid_conf = make_test_configuration();
  VERIFY(validate(valid_conf, true));
  VERIFY(validate(valid_conf, false));
  VERIFY(!validate(invalid_conf, true));
  VERIFY(!validate(invalid_conf, false));
}
