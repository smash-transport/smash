include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

add_library(divonne STATIC
    Divonne.c
    ${COMMON_DEPS}
)
