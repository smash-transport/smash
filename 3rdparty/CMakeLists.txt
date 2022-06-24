########################################################
#
#    Copyright (c) 2014,2016-2018,2021-2022
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# Let 3rd-party library use their own C++ standard
unset(CMAKE_CXX_STANDARD)

message(STATUS "Setting up Cuba library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
add_subdirectory(Cuba-4.2.1)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "Setting up Cuba library - done") # From v3.17 one should use message(CHECK_START ...)

message(STATUS "Setting up Einhard library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
add_subdirectory(einhard)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "Setting up Einhard library - done")

# Avoid YAML to build unnecessary parts of it
message(STATUS "Setting up YAML library")
list(APPEND CMAKE_MESSAGE_INDENT "   ")
set(YAML_CPP_BUILD_TESTS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_TOOLS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE INTERNAL "")
add_subdirectory(yaml-cpp-0.7.0)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "Setting up YAML library - done")

# Restore SMASH C++ standard
set(CMAKE_CXX_STANDARD ${CXX_SMASH_STANDARD})
