########################################################
#
#    Copyright (c) 2014,2016-2018,2021-2024
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# Let 3rd-party library use their own C++ standard
unset(CMAKE_CXX_STANDARD)

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
add_subdirectory(yaml-cpp-0.8.0)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "done")

# Restore SMASH C++ standard
set(CMAKE_CXX_STANDARD ${CXX_SMASH_STANDARD})

# To make sure that we do not get warnings from YAML, we want to treat its include as SYSTEM. For
# Einhard and Cuba we have changed their CMakeLists files and set it up there already, but for YAML
# we do not want do so and hence we adjust the target here.
set_property(TARGET yaml-cpp PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                                      $<TARGET_PROPERTY:yaml-cpp,INTERFACE_INCLUDE_DIRECTORIES>)

# Install 3rdparty libraries building their targets first. Note that libraries are installed both as
# static (ARCHIVE) and as shared (LIBRARY) so that whatever built versions get copied.
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND}
                                      --build ${CMAKE_BINARY_DIR}
                                      --target cuhre divonne suave vegas einhard yaml-cpp)")
install(FILES Cuba-4.2.2/cuba.h DESTINATION "include/cuba")
install(FILES einhard/einhard.hpp einhard/stacktrace.h DESTINATION "include/einhard")
install(DIRECTORY yaml-cpp-0.8.0/include/yaml-cpp
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h")
install(TARGETS cuhre
                divonne
                suave
                vegas
                einhard
                yaml-cpp
        ARCHIVE DESTINATION "lib"
        LIBRARY DESTINATION "lib")
