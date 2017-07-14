include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

add_library(suave STATIC
    Suave.c
    ${COMMON_DEPS}
)
