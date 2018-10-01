/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include <map>
#include <string>
#include "../include/smash/nucleus.h"
#include "../include/smash/pdgcode.h"

using namespace smash;

// we'll have equal particle tyoes for now so that center() ==
// center_of_mass (which is not yet implemented).
std::map<PdgCode, int> list = {{0x2212, 208}};

constexpr bool PRINT = false;

// we want to find the difference between initializing a vanilla
// woods-saxon distribution and a centered distribution, where the
// vanilla distribution is shifted so that its center-of-mass is at
// 0/0/0.
TEST(woods_saxon) {
  // this is where we store the distributions. Components (0,1,2,3) =
  // (x,y,z,r)
  std::map<int, int> hist_vanilla[4]{};
  std::map<int, int> hist_centerd[4]{};
  char name[4] = {'x', 'y', 'z', 'r'};
  // binning width for the distribution:
  constexpr double dx = 0.05;
  // this is the number of nuclei we create.
  constexpr int N_TEST = 10000;
  for (int i = 0; i < N_TEST; i++) {
    // initialize nucleus.
    Nucleus projectile(list, 1);
    projectile.arrange_nucleons();
    // com = center of mass.
    FourVector com = projectile.center();
    // for all particles in the nucleus, save position in histograms.
    for (auto p : projectile) {
      double r = p.position().abs3();
      ++hist_vanilla[0][floor(p.position().x1() / dx)];
      ++hist_vanilla[1][floor(p.position().x2() / dx)];
      ++hist_vanilla[2][floor(p.position().x3() / dx)];
      ++hist_vanilla[3][r / dx];
      // we'll "center" the nucleus "by hand", i.e., we subtract com by
      // hand instead of doing it in a dedicated function / loop.
      FourVector centered = p.position() - com;
      double R = centered.abs3();
      ++hist_centerd[0][floor(centered.x1() / dx)];
      ++hist_centerd[1][floor(centered.x2() / dx)];
      ++hist_centerd[2][floor(centered.x3() / dx)];
      ++hist_centerd[3][R / dx];
    }
  }
  constexpr int sigmabins = 4;
  // these are all a little smaller than erf((i+1)/sqrt(2)).
  constexpr double allowed[sigmabins] = {.682 * .99, .954 * .99, .997 * .99,
                                         1.0};
  // mnemonic: c = component
  for (int c = 0; c < 4; ++c) {
    int diffbad[sigmabins] = {0};
    int total = 0;
    for (auto b : hist_vanilla[c]) {
      // entries from vanilla distribution
      int vanilla = b.second;
      // corresponding entry in centered distribution
      int centerd = hist_centerd[c][b.first];
      // In order to make sensible comparisons, we require the number of
      // entries in both histograms to be >= 10.
      if (vanilla < 10 || centerd < 10) {
        continue;
      }
      // margin we allow. The statistical error at each point is
      // sqrt(N) (= N/sqrt(N)); the statistical error for the difference
      // between two points is hence sqrt(N1 + N2).
      double margin = std::sqrt(vanilla + centerd);
      // we'll make another histogram from the errors divided by the
      // standard deviation (stored in margin).
      int diffbin = std::abs(vanilla - centerd) / margin;
      // everything with a larger deviation than sigmabins*margin is
      // collected in highest bin.
      diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
      ++diffbad[diffbin];
      ++total;
      // maybe we want to print the distributions:
      if (PRINT) {
        std::printf("%d %7.3f %8d %8d\n", c, b.first * dx, vanilla, centerd);
      }
    }
    // if we don't want to print the distributions, then we want to
    // verify them.
    if (!PRINT) {
      // this is the integral over diffbad.
      int totalbad = 0;
      // now, we check diffbad one by one (don't need to check the last
      // one)
      for (int unit = 0; unit < sigmabins - 1; ++unit) {
        totalbad += diffbad[unit];
        double fraction = (totalbad + 0.0) / (total + 0.0);
        VERIFY(fraction > allowed[unit])
            << "\ntoo few entries have less than " << unit + 1
            << " sigma deviation\n(" << totalbad << "/" << total << "="
            << fraction << ", required minimal fraction: " << allowed[unit]
            << ")\nat component " << name[c];
      }
    }
  }
}
