cmake_minimum_required(VERSION 3.28 FATAL_ERROR)
project(SMASH_as_a_library_example)

add_library(project_requirements INTERFACE)
target_compile_features(project_requirements INTERFACE cxx_std_17)
target_compile_options(project_requirements
                       INTERFACE -W
                                 -Wall
                                 -Wextra
                                 -Wmissing-declarations
                                 -Wpointer-arith
                                 -Wshadow
                                 -Wuninitialized
                                 -Winit-self
                                 -Wundef
                                 -Wcast-align
                                 -Wformat=2
                                 -Wold-style-cast
                                 -Werror=switch
                                 -fno-math-errno)

add_executable(example example.cc)
add_executable(example_rate_equations example_rate_equations.cc)
add_executable(example_smash_wrapper example_smash_wrapper.cc)

find_package(SMASH 3.3 REQUIRED COMPONENTS Core) # Possible components: Core, ROOT, HepMC3, Rivet
if(SMASH_FOUND)
    # Finding SMASH sets SMASH_INPUT_FILES_DIR CMake variable
    target_compile_definitions(project_requirements
                               INTERFACE SMASH_INPUT_DIR=\"${SMASH_INPUT_FILES_DIR}\")
    target_link_libraries(example project_requirements SMASH::smash)
    target_link_libraries(example_rate_equations project_requirements SMASH::smash)
    target_link_libraries(example_smash_wrapper project_requirements SMASH::smash)
    # If you need to use the optional components, you can link to them as well. For example, if you
    # want to use ROOT and HepMC3, you can do: target_link_libraries(example project_requirements
    # SMASH::smash_root SMASH::smash_hepmc)
else()
    message(FATAL_ERROR "SMASH library not found!")
endif()
