########################################################
#
#    Copyright (c) 2014,2016-2018,2021-2023
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# Let 3rd-party library use their own C++ standard
unset(CMAKE_CXX_STANDARD)

message(STATUS "Setting up Cuba library") # From v3.17 one should use message(CHECK_START ...)
list(APPEND CMAKE_MESSAGE_INDENT "   ")
add_subdirectory(Cuba-4.2.2)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "Setting up Cuba library - done")

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

# Install 3rdparty libraries building their targets first. Note that libraries are installed both as
# static (ARCHIVE) and as shared (LIBRARY) so that whatever built versions get copied.
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND}
                                      --build ${CMAKE_BINARY_DIR}
                                      --target cuhre divonne suave vegas einhard yaml-cpp)")
install(FILES Cuba-4.2.2/cuba.h DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}/cuba")
install(FILES einhard/einhard.hpp einhard/stacktrace.h
        DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}/einhard")
install(DIRECTORY yaml-cpp-0.7.0/include/yaml-cpp
        DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}"
        FILES_MATCHING
        PATTERN "*.h")
install(TARGETS cuhre
                divonne
                suave
                vegas
                einhard
                yaml-cpp
        ARCHIVE DESTINATION "lib/${SMASH_INSTALLATION_SUBFOLDER}"
        LIBRARY DESTINATION "lib/${SMASH_INSTALLATION_SUBFOLDER}")
