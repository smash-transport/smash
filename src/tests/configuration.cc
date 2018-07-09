/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/configuration.h"
#include "../include/smash/forwarddeclarations.h"
#include "../include/smash/macros.h"

#include <boost/filesystem.hpp>

using namespace smash;

static Configuration make_test_configuration() {
  return Configuration{bf::path{TEST_CONFIG_PATH} / "src" / "tests",
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
  Configuration modi = conf["tamer"];
  double d = modi.take({"Altaic", "Meccas"});
  COMPARE(d, 10.);
  d = modi.take({"Altaic", "Kathleen"});
  COMPARE(d, 0.2);
  int i = modi.take({"Altaic", "Brahmins"});
  COMPARE(i, 1);
}

TEST_CATCH(take_incorrect_type, Configuration::IncorrectTypeInAssignment) {
  Configuration conf = make_test_configuration();
  Configuration modi = conf["tamer"];
  int i = modi.take({"pipit", "bushelling"});
  COMPARE(i, 5);
}

TEST(take_always_converts_to_string) {
  Configuration conf = make_test_configuration();
  Configuration modi = conf["tamer"];
  std::string s = modi.take({"pipit", "bushelling"});
  COMPARE(s, "5.0");
}

TEST(has_value) {
  Configuration conf = make_test_configuration();
  Configuration modi = conf["tamer"];
  VERIFY(modi.has_value({"pipit", "bushelling"}));
  VERIFY(modi.has_value({"pipit", "bushelling"}));
}

TEST(take_removes_entry) {
  Configuration conf = make_test_configuration();
  Configuration modi = conf["tamer"];
  VERIFY(modi.has_value({"pipit", "bushelling"}));
  modi.take({"pipit", "bushelling"});
  VERIFY(!modi.has_value({"pipit", "bushelling"}));
}

// Sorry, but I have to put this in the std namespace, otherwise it doesn't
// compile. That's because the << operator is called from inside the UnitTest
// namespace and all involved types are in the std namespace.
namespace std {
static ostream &operator<<(ostream &s, const vector<string> &v) {
  s << '{';
  for (const auto x : v) {
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
  Configuration modi = conf["tamer"];
  conf.take({"particles"});
  conf.take({"decaymodes"});
  conf.take({"fireballs", "extorting"});
  conf.take({"fireballs", "infection"});
  conf.take({"fireballs", "arena"});
  conf.take({"fireballs", "pendulous"});
  conf.take({"fireballs", "scudded"});
  conf.take({"fireballs", "firebrands"});
  conf.take({"fireballs", "joker"});
  conf.take({"fireballs", "classify"});
  modi.take({"Altaic", "Meccas"});
  modi.take({"Altaic", "Kathleen"});
  modi.take({"Altaic", "Brahmins"});
  modi.take({"feathered"});
  {
    std::istringstream unused(conf.unused_values_report());
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

  modi.take({"pipit", "bushelling"});
  {
    std::istringstream unused(conf.unused_values_report());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    expect_lines({"    neglects: -211", "    warbler: 211", "    reedier: 1.0"},
                 unused);
    VERIFY(unused.eof());
  }

  modi.take({"schmoozed", "warbler"});
  {
    std::istringstream unused(conf.unused_values_report());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    expect_lines({"    neglects: -211", "    reedier: 1.0"}, unused);
    VERIFY(unused.eof());
  }

  modi.take({"schmoozed", "reedier"});
  {
    std::istringstream unused(conf.unused_values_report());
    std::string line;
    getline(unused, line);
    COMPARE(line, "tamer:");
    getline(unused, line);
    COMPARE(line, "  schmoozed:");
    getline(unused, line);
    COMPARE(line, "    neglects: -211");
    VERIFY(unused.eof());
  }

  modi.take({"schmoozed", "neglects"});
  reference = "{}";
  COMPARE(conf.unused_values_report(), reference);
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

TEST(test_sub_config_objects) {
  Configuration conf = make_test_configuration();
  Configuration general = conf["fireballs"];
  const Configuration box = conf["tamer"]["Altaic"];
  VERIFY(general.has_value({"classify"}));
  int nevents = general.read({"classify"});
  VERIFY(general.has_value({"classify"}));
  COMPARE(nevents, 1);
  nevents = general.take({"classify"});
  VERIFY(!general.has_value({"classify"}));
  COMPARE(nevents, 1);
  COMPARE(double(box.read({"Meccas"})), 10.);
}

TEST(check_setting_new_value) {
  Configuration conf = make_test_configuration();
  VERIFY(!conf.has_value({"Test"}));
  conf["Test"] = 1.;
  VERIFY(conf.has_value({"Test"}));
  COMPARE(double(conf.read({"Test"})), 1.);
}

TEST(merge_override) {
  Configuration conf = make_test_configuration();
  COMPARE(int(conf.read({"fireballs", "arena"})), 1000);
  COMPARE(int(conf.read({"fireballs", "classify"})), 1);
  conf.merge_yaml("fireballs: { classify: 2 }");
  COMPARE(int(conf.read({"fireballs", "arena"})), 1000);
  COMPARE(int(conf.read({"fireballs", "classify"})), 2);
}

TEST(remove_all_but) {
  Configuration conf = make_test_configuration();
  Configuration modi = conf["tamer"];
  modi.remove_all_but("pipit");
  conf.remove_all_but("tamer");
  COMPARE(conf.to_string(), "tamer:\n  pipit:\n    bushelling: 5.0");
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
