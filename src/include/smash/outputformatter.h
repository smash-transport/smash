#ifndef OUTPUT_FORMATTER_H
#define OUTPUT_FORMATTER_H

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "smash/particledata.h"  // Assuming ParticleData is defined in the smash namespace

namespace smash {

// Function object to convert to ASCII
struct ASCII {
  using ReturnType = std::string;

  template <typename T>
  std::string operator()(const T& value) const {
    return std::to_string(value);
  }
  template <typename T>
  std::string operator()(const T& value, const char* format) const {
    char buffer[20];
    std::sprintf(buffer, format, value);
    return buffer;
  }
  std::string operator()(const std::string& value) const { return value; }
};

// Function object to convert to Binary
struct Binary {
  using ReturnType = const void*;

  template <typename T>
  const void* operator()(const T& value) const {
    return static_cast<const void*>(&value);
  }
};

template <typename Converter>
class OutputFormatter {
 public:
  OutputFormatter(const std::vector<std::string>& in_quantities)
      : quantities_(in_quantities) {
    for (const std::string& quantity : quantities_) {
      if (quantity == "t") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[0], "%g");
        });
      } else if (quantity == "x") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[1], "%g");
        });
      } else if (quantity == "y") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[2], "%g");
        });
      } else if (quantity == "z") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[3], "%g");
        });
      } else if (quantity == "mass") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.effective_mass(), "%g");
        });
      } else if (quantity == "p0") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[0], "%.9g");
        });
      } else if (quantity == "px") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[1], "%.9g");
        });
      } else if (quantity == "py") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[2], "%.9g");
        });
      } else if (quantity == "pz") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[3], "%.9g");
        });
      } else if (quantity == "pdg") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().string().c_str(), "%s");
        });
      } else if (quantity == "ID") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.id(), "%i");
        });
      } else if (quantity == "charge") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.type().charge(), "%i");
        });
      } else if (quantity == "ncoll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().collisions_per_particle,
                                  "%i");
        });
      } else if (quantity == "form_time") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.formation_time(), "%g");
        });
      } else if (quantity == "xsecfac") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.xsec_scaling_factor(), "%g");
        });
      } else if (quantity == "proc_id_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().id_process, "%i");
        });
      } else if (quantity == "proc_type_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(
              static_cast<int>(in.get_history().process_type), "%i");
        });
      } else if (quantity == "t_last_coll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().time_last_collision, "%g");
        });
      } else if (quantity == "pdg_mother1") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().p1.string().c_str(), "%s");
        });
      } else if (quantity == "pdg_mother2") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().p2.string().c_str(), "%s");
        });
      } else if (quantity == "baryon_number") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().baryon_number(), "%i");
        });
      } else if (quantity == "strangeness") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().strangeness(), "%i");
        });
      } else if (quantity == "spin_projection") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.spin_projection(), "%i");
        });
      } else {
        throw std::invalid_argument("OutputFormatter: Unknown quantity: " +
                                    quantity);
      }
    }
  }

  std::string data_line(const ParticleData& data) {
    std::string line;
    for (const auto& getter : getters_) {
      line += getter(data) + " ";
    }
    if (!line.empty()) {
      line.pop_back();
      line += "\n";
    }
    return line;
  }

  std::string header() {
    std::string header;
    for (const std::string& q : quantities_) {
      header += q + " ";
    }
    header.pop_back();
    return header;
  }

  std::string unit_line() {
    std::string line;
    for (const std::string& q : quantities_) {
      line += units_.at(q) + " ";
    }
    line.pop_back();
    return line;
  }
  // unit line (use map?)
 private:
  Converter converter_;
  std::vector<std::string> quantities_;
  std::vector<
      std::function<typename Converter::ReturnType(const ParticleData&)>>
      getters_;
  const std::map<std::string, std::string> units_ = {
      {"t", "fm"},
      {"x", "fm"},
      {"y", "fm"},
      {"z", "fm"},
      {"mass", "GeV"},
      {"p0", "GeV"},
      {"px", "GeV"},
      {"py", "GeV"},
      {"pz", "GeV"},
      {"pdg", "none"},
      {"ID", "none"},
      {"charge", "e"},
      {"ncoll", "none"},
      {"form_time", "fm"},
      {"xsecfac", "none"},
      {"proc_id_origin", "none"},
      {"proc_type_origin", "none"},
      {"time_last_coll", "fm"},
      {"pdg_mother1", "none"},
      {"pdg_mother2", "none"},
      {"baryon_number", "none"},
      {"strangeness", "none"},
      {"spin_projection", "none"},
  };
  // REN: To be filled
  const std::vector<std::string> OSCAR2013_quantities_ = {
      "t",  "x",  "y",  "z",   "mass", "p0",
      "px", "py", "pz", "pdg", "ID",   "charge"};
  const std::vector<std::string> OSCAR2013Extended_quantities_ = {
      "t",
      "x",
      "y",
      "z",
      "mass",
      "p0",
      "px",
      "py",
      "pz",
      "pdg",
      "ID",
      "charge",
      "ncoll",
      "form_time",
      "xsecfac",
      "proc_id_origin",
      "proc_type_origin",
      "time_last_coll",
      "pdg_mother1",
      "pdg_mother2",
      "baryon_number",
      "strangeness",
      "spin_projection"};
  // not sure what to do with "0"
  // const std::vector<std::string> OSCAR1999_quantities_ =
  // {"id","pdg","0","px","py","pz","p0","mass","x","y","z","t"};
};

}  // namespace smash

#endif  // OUTPUT_FORMATTER_H
