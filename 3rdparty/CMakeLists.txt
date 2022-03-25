add_subdirectory(einhard)
# Avoid YAML to build unnecessary parts of it
set(YAML_CPP_BUILD_TESTS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_TOOLS OFF CACHE INTERNAL "")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE INTERNAL "")
add_subdirectory(yaml-cpp-0.7.0)
add_subdirectory(Cuba-4.2.1)
