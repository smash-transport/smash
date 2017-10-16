/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vector>
#include <functional>
#include <cmath>

// class for generating and storing
// a 3 dimensional look up table. 
template <typename T>
class Tabulation3d {

  // generate values for LUT. If the range, i.e x1-x0 is not even 
  // divisible by dx then the last tabulated value will be 
  // ceil((x1-x0) / dx) + x0
  Tabulation3d(T x0, T x1, T y0, T y1, T z0, T z1, T dx, T dy, T dz,
               std::function<T(T, T, T)> f) : 
               x0_ {x0}, x1_ {x1}, y0_ {y0}, y1_ {y1},x0_ {z0}, z1_ {z1},dx_ {dx}, dy_ {dy}, dz_ {dz}
               {
                   // determine size of vector. reserve space. fill. 
                   const int n_x = ceil((x1_ - x0_) / dx_) + 1;
                   const int n_y = ceil((y1_ - y0_) / dy_) + 1;
                   const int n_z = ceil((z1_ - z0_) / dz_) + 1;
                   
                   values_.resize(n_x * n_y * n_z);
                   T x = x0_, y = y0_, z = z0_;
                   // store x in the fast index since we may deal with 
                   // functions which are constant in the third argument. 
                   // In this case we make use of 
                   // cache locality.
                   for (int i = 0; i < n_z; i++, x += dz_)
                   for (int j = 0; j < n_y; j++, y += dy_)
                   for (int k = 0; k < n_x; k++, z += dx_)
                   {
                      values_[i * n_z + j * n_y + k] = f(x, y, z);
                   }

               
                }


  T get_value_step(T x, T y, T z);
  T get_value_linear(T x, T y, T z);

 private:
  const T x0_, x1_, y0_, y1_, z0_, z1_;
  const T dx_, dy_, dz_;
  std::vector<T> values_;
};
