/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/configuration.h"

#include <cstdio>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>  // NOLINT(build/include_order)
#include <boost/filesystem.hpp>

#include "smash/forwarddeclarations.h"

namespace smash {

// internal helper functions
namespace {
/**
 * Finds a node, copies its structure and replaces the previous keys by
 * the newly provided keys.
 *
 * \param[in] node YAML::node whose structure will be copied and whose key
 * values will be replaced.
 * \param[in] keys Keys that will replace the previous keys in the copied node.
 * \return Node with the same structure as the input node but with updated keys.
 */
YAML::Node find_node_at(YAML::Node node,
                        std::initializer_list<const char *> keys) {
  assert(keys.size() > 0);
  for (auto key : keys) {
    // see comment in take on Node::reset
    node.reset(node[key]);
  }
  return node;
}

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

#ifndef DOXYGEN
namespace particles_txt {
#include <particles.txt.h>
}  // namespace particles_txt
namespace decaymodes_txt {
#include <decaymodes.txt.h>
}  // namespace decaymodes_txt
#endif

}  // unnamed namespace

// Default constructor
Configuration::Configuration(const bf::path &path)
    : Configuration(path, "config.yaml") {}

// Constructor checking for validity of input
Configuration::Configuration(const bf::path &path, const bf::path &filename) {
  const auto file_path = path / filename;
  if (!bf::exists(file_path)) {
    throw FileDoesNotExist("The configuration file was expected at '" +
                           file_path.native() +
                           "', but the file does not exist.");
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

  if (!root_node_["decaymodes"].IsDefined()) {
    root_node_["decaymodes"] = decaymodes_txt::data;
  }
  if (!root_node_["particles"].IsDefined()) {
    root_node_["particles"] = particles_txt::data;
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
  auto node = root_node_;
  auto keyIt = begin(keys);
  std::size_t i = 0;
  for (; i < keys.size() - 1; ++i, ++keyIt) {
    // Node::reset does what you might expect Node::operator= to do. But
    // operator= assigns a value to the node. So
    //   node = node[*keyIt]
    // leads to modification of the data structure, not simple traversal.
    node.reset(node[*keyIt]);
  }
  const auto r = node[*keyIt];
  node.remove(*keyIt);
  return {r, keys.begin()[keys.size() - 1]};
}

Configuration::Value Configuration::read(
    std::initializer_list<const char *> keys) const {
  return {find_node_at(root_node_, keys), keys.begin()[keys.size() - 1]};
}

void Configuration::remove_all_but(const std::string &key) {
  std::vector<std::string> to_remove;
  for (auto i : root_node_) {
    if (i.first.Scalar() != key) {
      to_remove.push_back(i.first.Scalar());
    }
  }
  for (auto i : to_remove) {
    root_node_.remove(i);
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
  s << remove_empty_maps(root_node_);
  return s.str();
}

std::string Configuration::to_string() const {
  std::stringstream s;
  s << root_node_;
  return s.str();
}

}  // namespace smash
