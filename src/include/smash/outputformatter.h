#ifndef OUTPUT_FORMATTER_H
#define OUTPUT_FORMATTER_H

#include <vector>
#include <string>
#include <functional>
#include "smash/particledata.h"  // Assuming ParticleData is defined in the smash namespace

namespace smash {

    // Function object to convert to ASCII
    struct ASCII {
        using ReturnType = std::string;

        template <typename T>
        std::string operator()(const T& value) const {
            return std::to_string(value);
        }
        std::string operator()(const std::string& value) const {
            return value;
        }
    };

    // Function object to convert to Binary
    struct Binary {
        using ReturnType = const void*;

        template <typename T>
        const void* operator()(const T& value) const {
            return static_cast<const void*>(&value);
        }
    };

    template<typename Converter>
    class OutputFormatter {
    public:
        OutputFormatter(const std::vector<std::string>& in_quantities)
        : quantities_(in_quantities) {
        for (const std::string& quantity : quantities_) {
            if (quantity == "id") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.id());
                    }
                );
            } 
            else if (quantity == "pdgcode") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.pdgcode().string());
                    }
                );
            } 
           
             else if (quantity == "formation_time") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.formation_time());
                    }
                );
            } 
            else if (quantity == "charge") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.type().charge());
                    }
                );
            }

         else if (quantity == "mass") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.type().mass());
                    }
                );
            } 
          
             else if (quantity == "spin") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.spin());
                    }
                );
            } 
            else if (quantity == "ncoll") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.get_history().collisions_per_particle);
                    }
                );
            } 
 
            else if (quantity == "proc_id_origin") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.get_history().id_process);
                    }
                );
            } 
 
            else if (quantity == "proc_type_origin") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(static_cast<int>(in.get_history().process_type));
                    }
                );
            }
            else if (quantity == "t_last_coll") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.get_history().time_last_collision);
                    }
                );
            }
            else if (quantity == "pdg_mother1") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.get_history().p1.string());
                    }
                );
            }
            else if (quantity == "pdg_mother2") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.get_history().p2.string());
                    }
                );
            }

            else if (quantity == "baryon_number") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.pdgcode().baryon_number());
                    }
                );
            }
            else if (quantity == "strangeness") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.pdgcode().strangeness());
                    }
                );
            }

            else if (quantity == "x0") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.position()[0]);
                    }
                );
            } 
            else if (quantity == "x1") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.position()[1]);
                    }
                );
            } 
            else if (quantity == "x2") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.position()[2]);
                    }
                );
            } 
            else if (quantity == "x3") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.position()[3]);
                    }
                );
            } 
            else if (quantity == "p0") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.momentum()[0]);
                    }
                );
            } 
            else if (quantity == "p1") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.momentum()[1]);
                    }
                );
            } 
            else if (quantity == "p2") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.momentum()[2]);
                    }
                );
            } 
            else if (quantity == "p3") {
                getters_.push_back(
                    [this](const ParticleData& in) {
                        return this->converter_(in.momentum()[3]);
                    }
                );
            } 
            else {
                throw std::invalid_argument("OutputFormatter: Unknown quantity: " + quantity);
            }
        }
    }



    std::string data_line(const ParticleData& data){
        std::string line;
        for(const auto& getter : getters_){
            line+= getter(data) + ",";

        }
        if (!line.empty()) {
            line.pop_back();
        }
        return line;
        
    }

    private:
        Converter converter_;
        std::vector<std::string> quantities_;
        std::vector<std::function<typename Converter::ReturnType(const ParticleData&)>> getters_;
    };

} // namespace smash

#endif // OUTPUT_FORMATTER_H
