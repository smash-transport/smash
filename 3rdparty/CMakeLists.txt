########################################################
#
#    Copyright (c) 2014,2016-2018,2021-2024,2026
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# Note that here we do not use smash_project_build_requirements and smash_project_usage_requirements
# interface libraries, because we do not want to force SMASH specific decisions onto vendored
# libraries coming from 3rd-party sources. For example, we do not want to force them to necessarily
# be built with the same C++ standard as SMASH, because this might not be always possible.

message(CHECK_START "Setting up Cuba library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
add_subdirectory(Cuba-4.2.2)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "done")

message(CHECK_START "Setting up Einhard library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
add_subdirectory(einhard)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "done")

message(CHECK_START "Setting up YAML library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
# Avoid YAML to build unnecessary parts of it
set(YAML_CPP_BUILD_TESTS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_TOOLS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE INTERNAL "")
add_subdirectory(yaml-cpp-0.9.0)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "done")

# To make sure that we do not get warnings from YAML, we want to treat its include as SYSTEM. For
# Einhard and Cuba we have changed their CMakeLists files and set it up there already, but for YAML
# we do not want do so and hence we adjust the target here.
set_property(TARGET yaml-cpp PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                                      $<TARGET_PROPERTY:yaml-cpp,INTERFACE_INCLUDE_DIRECTORIES>)

# We also want to make the vir-test library available to the rest of the project, because it
# contains the unit-tests framework that SMASH uses. Since virtest is a header only library, we can
# just create an interface library and attach the include directories to it.
add_library(vir-test INTERFACE)
target_include_directories(vir-test SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/virtest)

# Install third-party libraries relying on targets only. Note that header files are attached to
# targets for Cuba and Einhard and therefore we can refer to that. YAML library is a "vendored
# dependency" and it needs to be installed here as well and it would be wrong to use their shipped
# installation because we want that users using SMASH as a library link to the frozen YAML library
# shipped with SMASH.
install(DIRECTORY yaml-cpp-0.9.0/include/yaml-cpp
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h")
install(TARGETS yaml-cpp EXPORT SMASH-targets)
install(TARGETS cuhre divonne suave vegas
        EXPORT SMASH-targets
        FILE_SET HEADERS
        DESTINATION include/cuba)
install(TARGETS einhard
        EXPORT SMASH-targets
        FILE_SET HEADERS
        DESTINATION include/einhard)
install(EXPORT SMASH-targets NAMESPACE SMASH:: DESTINATION lib/cmake)
