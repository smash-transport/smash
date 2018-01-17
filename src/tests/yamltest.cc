/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include <yaml-cpp/yaml.h>
#include <iostream>

namespace YAML {
static Node operator|=(Node a, const Node &b) {
  if (b.IsMap()) {
    for (auto n0 : b) {
      a[n0.first.Scalar()] |= n0.second;
    }
  } else {
    a = b;
  }
  return a;
}
}  // namespace YAML

TEST(merge) {
  std::string yaml0 =
      "Foo:\n"
      "  Bar: Hallo\n"
      "  Baz: 1.2\n"
      "  Length: 1.2\n"
      "  Height: 21.8\n"
      "Test: 1\n";
  std::string yaml1 =
      "Foo:\n"
      "  Bar: Hallo Welt\n"
      "  Length: [1, 2, 3]\n"
      "  Extra:\n"
      "    radius: 5.2\n"
      "    height: 1.3\n"
      "Test: 2\n";
  YAML::Node a = YAML::Load(yaml0);
  YAML::Node b = YAML::Load(yaml1);
  std::cout << a << "\n\n";
  a |= b;
  std::cout << a << "\n\n";

  const auto foo = a["Foo"];
  const auto extra = foo["Extra"];
  COMPARE(foo["Bar"].Scalar(), "Hallo Welt");
  COMPARE(foo["Baz"].Scalar(), "1.2");
  COMPARE(foo["Length"].as<std::vector<int>>(), (std::vector<int>{1, 2, 3}));
  COMPARE(foo["Height"].Scalar(), "21.8");
  COMPARE(extra["radius"].Scalar(), "5.2");
  COMPARE(extra["height"].Scalar(), "1.3");
  COMPARE(a["Test"].Scalar(), "2");
}
