/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vector>

// class for generating a 3 dimensional look up table. 
template <typename T>
class Tabulation3d {
  Tabulation3d(T x0, T x1, T y0, T y1, T z0, T z1, T dx, T dy, T dz,
               std::function<T(T, T, T)> f) : 
               x0_ {x0}, x1_ {x1}, y0_ {y0}, y1_ {y1},x0_ {z0}, z1_ {z1},dx_ {dx}, dy_ {dy}, dz_ {dz}
               {
                   // determine size of vector. reserve space. fill. 
               }

  Tabulation3d(T x0, T x1, T y0, T y1, T z0, T z1, int n_points,
               std::function<T(T, T, T)> f);

  T get_value_step(T x, T y, T z);
  T get_value_linear(T x, T y, T z);

 private:
  const T x0_, x1_, y0_, y1_, z0_, z1_;
  const T dx_, dy_, dz_;
  std::vector<val_type> values_;
};
