include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

add_library(vegas STATIC
    Vegas.c
    ${COMMON_DEPS}
)
