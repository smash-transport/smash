/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/input_keys.h"

using namespace smash;

template <typename T>
static Key<T> get_test_key(std::optional<T> def_val = std::nullopt) {
  const Version v_intro = "1.0.0";
  return Key<T>{{"Test", "Key"}, def_val, {v_intro}};
}

TEST(create_keys) {
  get_test_key<int>();
  get_test_key<double>(3.14);
}

TEST_CATCH(ask_not_available_default, std::bad_optional_access) {
  const auto key = get_test_key<int>();
  key.default_value();
}

TEST_CATCH(key_without_version, smash::Key<int>::WrongNumberOfVersions) {
  Key<int> key({"Test"}, {});
}

TEST_CATCH(key_with_too_many_versions, smash::Key<int>::WrongNumberOfVersions) {
  Key<int> key({"Test"}, {"v1", "v2", "v3", "v4"});
}

TEST(check_introduction_version) {
  const auto key = get_test_key<int>();
  COMPARE(key.introduced_in(), "1.0.0");
}

TEST_CATCH(not_deprecated_key, std::bad_optional_access) {
  const auto key = get_test_key<int>();
  VERIFY(key.is_allowed());
  VERIFY(!key.is_deprecated());
  key.deprecated_in();
}

TEST(deprecated_key) {
  const Version v_deprecation = "1.2.0";
  const Key<std::string> key{{"Test"}, "Hello world", {"1.0.0", v_deprecation}};
  VERIFY(key.is_allowed());
  VERIFY(key.is_deprecated());
  COMPARE(key.deprecated_in(), v_deprecation);
}

TEST_CATCH(not_removed_key, std::bad_optional_access) {
  const auto key = get_test_key<int>();
  VERIFY(key.is_allowed());
  key.removed_in();
}

TEST(removed_key) {
  const Version v_deprecation = "1.2.0";
  const Version v_removal = "1.2.0";
  const Key<std::string> key{
      {"Test"}, "Hello world", {"1.0.0", v_deprecation, v_removal}};
  VERIFY(!key.is_allowed());
  COMPARE(key.removed_in(), v_removal);
}

TEST(has_same_labels) {
  const auto key = get_test_key<int>();
  VERIFY(key.has_same_labels({"Test", "Key"}));
  VERIFY(!key.has_same_labels({}));
  VERIFY(!key.has_same_labels({"Test"}));
  VERIFY(!key.has_same_labels({"Test", "Key", "Extra"}));
  VERIFY(!key.has_same_labels({"Test", "Other_Key"}));
}

TEST(to_string) {
  const auto key = get_test_key<int>();
  const std::string result = "\"Test: Key\"";
  COMPARE(static_cast<std::string>(key), result);
}
