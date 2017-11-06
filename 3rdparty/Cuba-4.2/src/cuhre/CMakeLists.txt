include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

add_library(cuhre STATIC
    Cuhre.c
    ${COMMON_DEPS}
)
