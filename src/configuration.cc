/*
 *
 *    Copyright (c) 2014-2019
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
#include "smash/inputfunctions.h"

namespace smash {

// internal helper functions
namespace {
/**
 * Removes all empty maps of a YAML::Node.
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
 * Builds a string with a list of keys as specified in the code.
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

std::string Configuration::unused_values_report() const {
  std::stringstream s;
  s << root_node_;
  return s.str();
}

std::string Configuration::to_string() const {
  std::stringstream s;
  s << root_node_;
  return s.str();
}

YAML::Node Configuration::find_node_at(YAML::Node node,
                                       std::vector<const char *> keys) const {
  assert(keys.size() > 0);
  for (auto key : keys) {
    // Node::reset does what you might expect Node::operator= to do. But
    // operator= assigns a value to the node. So
    //   node = node[*keyIt]
    // leads to modification of the data structure, not simple traversal.
    node.reset(node[key]);
  }
  return node;
}

}  // namespace smash
