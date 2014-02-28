/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARAMETERS_H_
#define SRC_INCLUDE_PARAMETERS_H_

#include <cstring>
#include <cstdlib>

namespace Smash {

class Parameters {
 public:
  /* default constructor */
  Parameters(char *k, char *v) : key_(strdup(k)), value_(strdup(v)) {}
  /* copy constructor */
  Parameters(const Parameters &p)
      : key_(strdup(p.key_)), value_(strdup(p.value_)) {}
  /* destructor */
  ~Parameters() {
    if (key_) {
      free(key_);
    }
    if (value_) {
      free(value_);
    }
  }
  /* accessors */
  inline void set_key(char *key_string);
  inline char *key(void) const;
  inline void set_value(char *value_string);
  inline char *value(void) const;

 private:
  char *key_;
  char *value_;
};

inline char *Parameters::key(void) const {
  return key_;
}

inline void Parameters::set_key(char *key_string) {
  if (key_) {
    free(key_);
  }
  key_ = strdup(key_string);
}

inline char *Parameters::value(void) const {
  return value_;
}

inline void Parameters::set_value(char *value_string) {
  if (value_) {
    free(value_);
  }
  value_ = strdup(value_string);
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PARAMETERS_H_
