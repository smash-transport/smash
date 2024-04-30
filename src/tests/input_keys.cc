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
using namespace std::string_literals;

template <typename T>
static Key<T> get_test_key(std::optional<T> def_val = std::nullopt) {
  const Version v_intro = "1.0.0";
  if (def_val)
    return Key<T>{{"Test", "Key"}, *def_val, {v_intro}};
  else
    return Key<T>{{"Test", "Key"}, {v_intro}};
}

TEST(create_keys) {
  get_test_key<int>();
  get_test_key<double>(3.14);
}

TEST_CATCH(key_with_invalid_default_type_I, std::logic_error) {
  Key<int> key({"Test"}, DefaultType::Null, {"1.0"});
}

TEST_CATCH(key_with_invalid_default_type_II, std::logic_error) {
  Key<int> key({"Test"}, DefaultType::Value, {"1.0"});
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
  const Key<std::string> key{
      {"Test"}, "Hello world"s, {"1.0.0", v_deprecation}};
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
      {"Test"}, "Hello world"s, {"1.0.0", v_deprecation, v_removal}};
  VERIFY(!key.is_allowed());
  COMPARE(key.removed_in(), v_removal);
}

TEST(dependent_default) {
  const auto key_1 = get_test_key<int>();
  VERIFY(!key_1.has_dependent_default());
  const auto key_2 = get_test_key<int>(42);
  VERIFY(!key_2.has_dependent_default());
  Key<double> dependent_key{{"Test"}, DefaultType::Dependent, {"1.0.0"}};
  VERIFY(dependent_key.has_dependent_default());
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

/*

// The following code is useful to print all keys in the database for debug
// purposes and it is intentionally left as part of the codebase commented out
// for future needs.

template <typename T, typename = void>
auto constexpr ostreamable_v = false;

template <typename T>
auto constexpr ostreamable_v<
    T, std::void_t<decltype(std::cout << std::declval<T>())>> = true;

template <typename T>
std::enable_if_t<ostreamable_v<T>> print(T const& in) {
  std::cout << in << '\n';
}

template <typename T>
std::enable_if_t<!ostreamable_v<T>> print(T const& in) {
  std::cout << "NOT PRINTABLE (" << typeid(in).name() << ")\n";
}

TEST(check_list) {
  for (const auto& key : InputKeys::list) {
    try {
      std::visit(
          [](auto&& arg) {
            std::cout << std::setw(70) << std::string(arg.get()) << "  ";
            if (arg.get().has_dependent_default())
              std::cout << "--> dependent default!\n";
            else {
              const auto def = arg.get().default_value();
              print(def);
            }
          },
          key);
    } catch (const std::bad_optional_access&) {
      std::visit(
          [](auto&&) {
            std::cout << "--> Mandatory\n";
          },
          key);
    }
  }
}

*/
