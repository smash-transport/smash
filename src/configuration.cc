/*
 *
 *    Copyright (c) 2014-2019,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/configuration.h"

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "smash/forwarddeclarations.h"
#include "smash/input_keys.h"
#include "smash/inputfunctions.h"
#include "smash/logging.h"
#include "smash/stringfunctions.h"

namespace smash {
static constexpr int LConf = LogArea::YAML_Configuration::id;

// internal helper functions
namespace {
/**
 * Remove all empty maps of a YAML::Node.
 *
 * \param[in] root YAML::Node that contains empty maps.
 * \return YAML::Node from above without empty maps.
 */
YAML::Node remove_empty_maps(YAML::Node root) {
  if (root.IsMap()) {
    std::vector<std::string> to_remove(root.size());
    for (auto n : root) {
      remove_empty_maps(n.second);
      if ((n.second.IsMap() || n.second.IsSequence()) && n.second.size() == 0) {
        to_remove.emplace_back(n.first.Scalar());
      }
    }
    for (const auto &key : to_remove) {
      root.remove(key);
    }
  }
  return root;
}

/**
 * Merge two YAML::Nodes
 *
 * \param[in] a YAML::Node into which b is merged.
 * \param[in] b YAML::Node that is merged into a.
 * \return YAML::Node which is the merge of a and b.
 */
YAML::Node operator|=(YAML::Node a, const YAML::Node &b) {
  if (b.IsMap()) {
    for (auto n0 : b) {
      a[n0.first.Scalar()] |= n0.second;
    }
  } else {
    a = b;
  }
  return a;
}

/**
 * Build a string with a list of keys as specified in the code.
 *
 * @param keys The list of keys.
 * @return A \c std::string with the desired result.
 */
std::string join_quoted(std::initializer_list<const char *> keys) {
  return std::accumulate(keys.begin(), keys.end(), std::string{"{"},
                         [](const std::string &ss, const std::string &s) {
                           return ss + ((ss.size() == 1) ? "\"" : ", \"") + s +
                                  "\"";
                         }) +
         "}";
}

}  // unnamed namespace

// Default constructor
Configuration::Configuration(const std::filesystem::path &path)
    : Configuration(path, "config.yaml") {}

// Constructor checking for validity of input
Configuration::Configuration(const std::filesystem::path &path,
                             const std::filesystem::path &filename) {
  const auto file_path = path / filename;
  if (!std::filesystem::exists(file_path)) {
    throw FileDoesNotExist("The configuration file was expected at '" +
                           file_path.native() +
                           "', but the file does not exist.");
  }
  if (has_crlf_line_ending(read_all(std::ifstream((file_path))))) {
    throw std::runtime_error(
        "The configuration file has CR LF line endings. Please use LF "
        "line endings.");
  }
  try {
    root_node_ = YAML::LoadFile(file_path.native());
  } catch (YAML::ParserException &e) {
    if (e.msg == "illegal map value" || e.msg == "end of map not found") {
      const auto line = std::to_string(e.mark.line + 1);
      throw ParseError("YAML parse error at\n" + file_path.native() + ':' +
                       line + ": " + e.msg +
                       " (check that the indentation of map keys matches)");
    }
    throw;
  }
}

void Configuration::merge_yaml(const std::string &yaml) {
  try {
    root_node_ |= YAML::Load(yaml);
  } catch (YAML::ParserException &e) {
    if (e.msg == "illegal map value" || e.msg == "end of map not found") {
      const auto line = std::to_string(e.mark.line + 1);
      throw ParseError("YAML parse error in:\n" + yaml + "\nat line " + line +
                       ": " + e.msg +
                       " (check that the indentation of map keys matches)");
    }
    throw;
  }
}

std::vector<std::string> Configuration::list_upmost_nodes() {
  std::vector<std::string> r;
  for (auto i : root_node_) {
    r.emplace_back(i.first.Scalar());
  }
  return r;
}

Configuration::Value Configuration::take(
    std::initializer_list<const char *> keys) {
  assert(keys.size() > 0);
  auto last_key_it = keys.end() - 1;
  auto node = find_node_at(root_node_, {keys.begin(), last_key_it});
  const auto r = node[*last_key_it];
  node.remove(*last_key_it);
  root_node_ = remove_empty_maps(root_node_);
  return {r, *last_key_it};
}

Configuration::Value Configuration::read(
    std::initializer_list<const char *> keys) const {
  return {find_node_at(root_node_, keys), keys.begin()[keys.size() - 1]};
}

void Configuration::remove_all_entries_in_section_but_one(
    const std::string &key, std::initializer_list<const char *> section) {
  auto node = find_node_at(root_node_, section);
  std::vector<std::string> to_remove{};
  for (auto i : node) {
    if (i.first.Scalar() != key) {
      to_remove.push_back(i.first.Scalar());
    }
  }
  for (auto i : to_remove) {
    node.remove(i);
  }
}

Configuration Configuration::extract_sub_configuration(
    std::initializer_list<const char *> keys,
    Configuration::GetEmpty empty_if_not_existing) {
  assert(keys.size() > 0);
  auto last_key_it = keys.end() - 1;
  auto node = find_node_at(root_node_, {keys.begin(), last_key_it});
  auto sub_conf_root_node = node[*last_key_it];
  if (!sub_conf_root_node.IsDefined()) {
    if (empty_if_not_existing == Configuration::GetEmpty::Yes)
      return Configuration(YAML::Node{});
    else
      throw std::runtime_error("Attempt to extract not existing section " +
                               join_quoted(keys));
  } else if (sub_conf_root_node.IsNull() ||
             (sub_conf_root_node.IsMap() && sub_conf_root_node.size() == 0)) {
    // Here we put together the cases of a key without value or with
    // an empty map {} as value (no need at the moment to distinguish)
    throw std::runtime_error("Attempt to extract empty section " +
                             join_quoted(keys));
  } else if (sub_conf_root_node.IsMap() && sub_conf_root_node.size() != 0) {
    Configuration sub_config{sub_conf_root_node};
    node.remove(*last_key_it);
    return sub_config;
  } else {  // sequence or scalar or any future new YAML type
    throw std::runtime_error("Tried to extract configuration section at " +
                             join_quoted(keys) +
                             " to get a key value. Use take instead!");
  }
}

bool Configuration::has_value_including_empty(
    std::initializer_list<const char *> keys) const {
  const auto n = find_node_at(root_node_, keys);
  return n.IsDefined();
}

bool Configuration::has_value(std::initializer_list<const char *> keys) const {
  const auto n = find_node_at(root_node_, keys);
  return n.IsDefined() && (!n.IsNull());
}

std::string Configuration::to_string() const {
  std::stringstream s;
  s << root_node_;
  return s.str();
}

YAML::Node Configuration::find_node_at(YAML::Node node,
                                       std::vector<const char *> keys) const {
  /* Here we do not assert(keys.size()>0) and allow to pass in an empty vector,
     in which case the passed in YAML:Node is simply returned. This might happen
     e.g. in the take or extract_sub_configuration methods if called with a
     label of a key at top level of the configuration file. */
  for (auto key : keys) {
    /* Node::reset does what you might expect Node::operator= to do. But
       operator= assigns a value to the node. So
         node = node[*keyIt]
       leads to modification of the data structure, not simple traversal. */
    node.reset(node[key]);
  }
  return node;
}

// internal helper functions
namespace {
/**
 * @brief Implementation of the algorithm to translate a configuration into
 * lists of labels, each identifying a key from the root YAML node.
 *
 * Since the level of nesting sections in a YAML input file is arbitrary, this
 * is a typical task to be solved using recursion. The main idea here is to
 * take advantage of Configuration functionality and in particular of its
 * Configuration::list_upmost_nodes method as well as the
 * Configuration::operator[]. In particular, from the configuration top-level
 * all upmost nodes are extracted and for each of them, recursively, the
 * same procedure is done over and over again. Before the recursive call, the
 * present section label is added to the list of labels in process to be
 * created and after the recursive call, the list of labels is added to the
 * main list if a value has been found in the previous recursive call.
 * The only aspect to pay attention to is to store for each recursive call
 * whether a value has been found and accordingly push a new element in the
 * main list or not. Such a \c bool flag is implemented as return value.
 *
 * @param[in] configuration The input configuration to extract from.
 * @param[inout] list The list of lists of labels to be filled.
 * @param[inout] new_list_entry New list of labels in process to be filled
 *                              during recursion.
 * @param[inout] push_to_main_list Whether to push or not to the main list.
 *                                 This is an array and not a single \c bool
 *                                 value, in order to keep track of all levels
 *                                 of recursion.
 * @return \c true if the passed configuration has no sections (just a value).
 * @return \c false if the passed configuration has sections.
 */
bool fill_list_of_labels_per_key_from_configuration_implementation(
    Configuration &configuration, std::vector<std::vector<std::string>> &list,
    std::vector<std::string> &new_list_entry,
    std::vector<bool> &push_to_main_list) {
  std::vector<std::string> list_of_sections{};
  try {
    list_of_sections = configuration.list_upmost_nodes();
  } catch (YAML::InvalidNode &) {
    /*
     * This happens for keys that are not sections and whose value
     * cannot be interpreted as a YAML Scalar, e.g. lists [2112, 2212].
     */
    return true;
  };

  if (list_of_sections.empty()) {
    // This happens when a key is not a section and has a Scalar value.
    return true;
  }

  for (const auto &section : list_of_sections) {
    new_list_entry.push_back(section);
    Configuration sub_configuration =
        configuration.extract_sub_configuration({section.c_str()});
    push_to_main_list.push_back(
        fill_list_of_labels_per_key_from_configuration_implementation(
            sub_configuration, list, new_list_entry, push_to_main_list));
    if (push_to_main_list.back()) {
      list.push_back(new_list_entry);
    }
    new_list_entry.pop_back();
    push_to_main_list.pop_back();
  }
  /* If the function is executed till here, it means that only sections and
   * not keys with values have been found. These are not to be pushed into the
   * main list. Hence, return false.
   */
  return false;
}

/**
 * @brief Clear and fill the given list object using the passed configuration.
 *
 * Given a \c Configuration object, for each YAML key having a value, all
 * labels to reach the given key from the root node are collected and a \c
 * std::vector containing them is built and inserted into the given list. This
 * function is calling the actual implementation preparing auxiliary needed
 * variables.
 *
 * @param[in] configuration The input configuration to extract from.
 * @param[inout] list The list of lists of labels to be filled.
 */
void fill_list_of_labels_per_key_from_configuration(
    Configuration &configuration, std::vector<std::vector<std::string>> &list) {
  std::vector<std::string> aux_v1{};
  std::vector<bool> aux_v2{};
  list.clear();
  SMASH_UNUSED(fill_list_of_labels_per_key_from_configuration_implementation(
      configuration, list, aux_v1, aux_v2));
}

/**
 * @brief Given some YAML labels (assumed to be in order from the top
 * section), it is checked whether any valid SMASH key with the same key
 * exists.
 *
 * All possible checks are done in a way such that the user is informed about
 *  - if the key has never been valid;
 *  - if the key was valid in the past but it has been removed;
 *  - if the key is valid but deprecated.
 *
 * @param[in] labels The series of labels identifying the key.
 *
 * @return \c true if the key is valid and
 * @return \c false otherwise.
 */
bool validate_key_based_on_labels(const std::vector<std::string> &labels) {
  using key_ref_var = smash::InputKeys::key_references_variant;
  auto key_ref_var_it = std::find_if(
      smash::InputKeys::list.begin(), smash::InputKeys::list.end(),
      [&labels](auto key) {
        return std::visit(
            [&labels](auto &&arg) { return arg.get().has_same_labels(labels); },
            key);
      });
  if (key_ref_var_it == smash::InputKeys::list.end()) {
    logg[LConf].error("Key ", smash::quote(smash::join(labels, " -> ")),
                      " is not a valid SMASH input key.");
    return false;
  }

  key_ref_var found_variant = *key_ref_var_it;
  const auto key_labels =
      std::visit([](auto &&arg) { return static_cast<std::string>(arg.get()); },
                 found_variant);

  if (std::visit([](auto &&arg) { return !arg.get().is_allowed(); },
                 found_variant)) {
    const auto v_removal = std::visit(
        [](auto &&arg) { return arg.get().removed_in(); }, found_variant);
    logg[LConf].error("Key ", key_labels, " has been removed in version ",
                      v_removal, " and it is not valid anymore.");
    return false;
  }
  if (std::visit([](auto &&arg) { return arg.get().is_deprecated(); },
                 found_variant)) {
    const auto v_deprecation = std::visit(
        [](auto &&arg) { return arg.get().deprecated_in(); }, found_variant);
    logg[LConf].warn("Key ", key_labels, " has been deprecated in version ",
                     v_deprecation);
  } else
    logg[LConf].debug("Key ", key_labels, " valid!");
  return true;
}

}  // namespace

bool validate(Configuration &configuration, bool full_validation) {
  std::vector<std::vector<std::string>> list{};
  fill_list_of_labels_per_key_from_configuration(configuration, list);
  if (full_validation)
    return std::transform_reduce(
        list.begin(), list.end(), bool{true},
        [](bool result, bool value) { return result && value; },
        validate_key_based_on_labels);
  else {
    return std::find_if_not(list.begin(), list.end(),
                            validate_key_based_on_labels) == list.end();
  }
}

}  // namespace smash
