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
 * Implementation of the algorithm to translate a YAML tree into
 * lists of labels, each identifying a key from the YAML root node.
 *
 * Since the level of nesting sections in a YAML input file is arbitrary, this
 * is a typical task to be solved using recursion. The main idea here is to
 * take advantage of YAML functionality and in particular of the possibility
 * to iterate over trees and test for nature of a node (is it a Map or not?).
 * Roughly speaking, from the tree top-level all upmost nodes are extracted
 * and for each of them, recursively, the same procedure is done over and
 * over again if they are maps. If a non-map node is found, i.e. a key value
 * is found, then recursion ends a new entry to \c list is added.
 *
 * @param[in] root_node The root YAML node to extract from.
 * @param[inout] list The list of lists of labels to be filled.
 * @param[inout] new_list_entry New list of labels in process to be filled
 *                              during recursion.
 */
void fill_list_of_labels_per_key_in_yaml_tree(const YAML::Node &root_node,
                                              std::vector<KeyLabels> &list,
                                              KeyLabels &new_list_entry) {
  // Here sub_node is an iterator value, i.e. a key/value pair of nodes,
  // not a single YAML node (that's how YAML library works)
  for (const auto &sub_node : root_node) {
    new_list_entry.push_back(sub_node.first.as<std::string>());
    if (sub_node.second.IsMap())
      fill_list_of_labels_per_key_in_yaml_tree(sub_node.second, list,
                                               new_list_entry);
    else
      list.push_back(new_list_entry);
    new_list_entry.pop_back();
  }
}

/**
 * Create a list of lists of key labels present in the passed YAML node
 * considered to be the root one of a YAML tree.
 *
 * Given a \c YAML::Node, for each key having a value, all labels to reach
 * the given key from the passed node are collected and a \c std::vector
 * containing them is built and inserted into the given list. This function
 * is calling the actual implementation preparing auxiliary needed variables.
 *
 * @param[in] root_node The root node of the YAML tree to be considered.
 *
 * @return A \c std::vector<KeyLabels> containing the desired information.
 */
auto get_list_of_labels_per_key_in_yaml_tree(const YAML::Node &root_node) {
  std::vector<KeyLabels> list{};
  KeyLabels aux{};
  fill_list_of_labels_per_key_in_yaml_tree(root_node, list, aux);
  return list;
}

#if 1
/**
 * \brief A utility type to be specialized to check if a type is a \c std::map .
 *
 * \tparam T A generic template parameter.
 */
template <typename T>
struct IsStdMap {
  /**
   * A boolean value to indicate whether \c T is a map or not. Here it is always
   * \c false, because there is another template specialization that will select
   * those cases where \c value is going to be \c true.
   */
  static constexpr bool value = false;
};

/**
 * \brief A specialization of \c IsStdMap<T> for cases where the
 * boolean value should be set to \c true.
 *
 * \tparam MapKey A type to indicate map keys.
 * \tparam MapValue A type to indicate map values.
 */
template <typename MapKey, typename MapValue>
struct IsStdMap<std::map<MapKey, MapValue>> {
  /**
   * A boolean value to indicate whether \c T is a map or not. Here it is always
   * \c true, because this is a template specialization for maps only.
   */
  static constexpr bool value = true;
};

/**
 * \brief Extract from the \c InputKeys database the labels of keys that
 * have a \c std::map as type.
 *
 * \return A list of key labels.
 */
auto collect_input_keys_taken_as_maps() {
  std::vector<KeyLabels> labels_of_keys_taken_as_map{};
  for (const auto &keys_variant : smash::InputKeys::list) {
    std::visit(
        [&labels_of_keys_taken_as_map](auto &&var) {
          /*
           * The following if checks if the SMASH input key has a map as value
           * and it deserves some explanation about the type extraction:
           *
           *   - arg -> object of type: std::cref(const Key<T>)
           *   - arg.get() -> object of type: const Key<T>&
           *   - decltype(arg.get()) ->  type: const Key<T>&
           *   - std::decay_t<decltype(arg.get())>::type -> type: Key<T>
           *   - std::decay_t<decltype(arg.get())>::type::value -> type: T
           */
          if constexpr (IsStdMap<typename std::decay_t<
                            decltype(var.get())>::type>::value)
            labels_of_keys_taken_as_map.push_back(var.get().labels());
        },
        keys_variant);
  }
  return labels_of_keys_taken_as_map;
}
#endif

/**
 * \brief Remove last labels of keys that are taken as maps in SMASH and remove
 * duplicates from the resulting list.
 *
 * The keys that are taken as maps in SMASH are here collected using the
 * database \c InputKeys and the list of keys contained in the configuration
 * must be adjusted by hand. This is a corner case, since YAML nodes that are
 * maps are **by definition** sections and cannot be distinguished from keys
 * "with a map value" in the recursive process to create the list of key labels.
 *
 * \param[in,out] list_of_input_key_labels The list of key labels to adjust.
 */
void adjust_list_of_labels_dealing_with_keys_taken_as_maps(
    std::vector<KeyLabels> &list_of_input_key_labels) {
  const std::vector<KeyLabels> labels_of_keys_taken_as_map =
      collect_input_keys_taken_as_maps();
  for (const auto &labels : labels_of_keys_taken_as_map) {
    std::for_each(list_of_input_key_labels.begin(),
                  list_of_input_key_labels.end(),
                  [&labels](KeyLabels &labels_of_input_key) {
                    if (std::equal(labels.begin(), labels.end(),
                                   labels_of_input_key.begin(),
                                   labels_of_input_key.begin() + labels.size()))
                      labels_of_input_key = labels;
                  });
  }
  // The identical keys in list are now next to each other and we do
  // not need/want to sort the list before calling std::unique.
  list_of_input_key_labels.erase(std::unique(list_of_input_key_labels.begin(),
                                             list_of_input_key_labels.end()),
                                 list_of_input_key_labels.end());
}

/**
 * Given some YAML labels (assumed to be in order from the top section),
 * it is checked whether any valid SMASH key with the same key exists.
 *
 * All possible checks are done in a way such that the user is informed about
 *  - if the key has never been valid;
 *  - if the key was valid in the past but it has been removed;
 *  - if the key is valid but deprecated.
 *
 * \param[in] labels The series of labels identifying the key.
 *
 * \return \c Configuration::Is::Valid if the key is valid;
 * \return \c Configuration::Is::Deprecated if the key is Deprecated and
 * \return \c Configuration::Is::Invalid if the key is invalid.
 */
Configuration::Is validate_key(const KeyLabels &labels) {
  auto key_ref_var_it = std::find_if(
      smash::InputKeys::list.begin(), smash::InputKeys::list.end(),
      [&labels](auto key) {
        return std::visit(
            [&labels](auto &&arg) { return arg.get().has_same_labels(labels); },
            key);
      });
  if (key_ref_var_it == smash::InputKeys::list.end()) {
    logg[LConf].error("Key ", smash::quote(smash::join(labels, ": ")),
                      " is not a valid SMASH input key.");
    return Configuration::Is::Invalid;
  }

  smash::InputKeys::key_references_variant found_variant = *key_ref_var_it;
  const auto key_labels =
      std::visit([](auto &&var) { return static_cast<std::string>(var.get()); },
                 found_variant);

  if (std::visit([](auto &&var) { return !var.get().is_allowed(); },
                 found_variant)) {
    const auto v_removal = std::visit(
        [](auto &&var) { return var.get().removed_in(); }, found_variant);
    logg[LConf].error("Key ", key_labels, " has been removed in version ",
                      v_removal, " and it is not valid anymore.");
    return Configuration::Is::Invalid;
  }
  if (std::visit([](auto &&var) { return var.get().is_deprecated(); },
                 found_variant)) {
    const auto v_deprecation = std::visit(
        [](auto &&var) { return var.get().deprecated_in(); }, found_variant);
    logg[LConf].warn("Key ", key_labels, " has been deprecated in version ",
                     v_deprecation);
    return Configuration::Is::Deprecated;
  } else {
    logg[LConf].debug("Key ", key_labels, " is valid!");
    return Configuration::Is::Valid;
  }
}

}  // namespace

Configuration::Is Configuration::validate(bool full_validation) const {
  auto list = get_list_of_labels_per_key_in_yaml_tree(root_node_);
  adjust_list_of_labels_dealing_with_keys_taken_as_maps(list);
  if (full_validation)
    return std::transform_reduce(
        list.begin(), list.end(), Is{Is::Valid},
        [](Is result_so_far, Is present_value) {
          switch (result_so_far) {
            case Is::Invalid:
              return result_so_far;
              break;
            case Is::Deprecated:
              return (present_value == Is::Valid) ? result_so_far
                                                  : present_value;
              break;
            case Is::Valid:
              return present_value;
              break;
          }
        },
        validate_key);
  else {
    for (const auto &key_labels : list) {
      Is key_state = validate_key(key_labels);
      if (key_state != Is::Valid)
        return key_state;
    }
    return Is::Valid;
  }
}

}  // namespace smash
