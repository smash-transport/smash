cmake_minimum_required(VERSION 2.6.2)
project(Cuba)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mfpmath=sse")

add_library(pythia
    src/common/Random.c
    src/common/Data.c
    src/common/CSample.c
    src/common/MSample.c
    src/common/Global.c
    src/common/Fork.c
    src/common/ChiSquare.c
    src/common/Parallel.c
    src/common/Erf.c
    src/vegas/Grid.c
    src/vegas/Vegas.c
    src/vegas/common.c
    src/vegas/Integrate.c
    src/cuhre/Cuhre.c
    src/cuhre/common.c
    src/cuhre/Integrate.c
    src/cuhre/Rule.c
    src/suave/Sample.c
    src/suave/Grid.c
    src/suave/common.c
    src/suave/Fluct.c
    src/suave/Integrate.c
    src/suave/Suave.c
    src/divonne/Sample.c
    src/divonne/KorobovCoeff.c
    src/divonne/Iterate.c
    src/divonne/Explore.c
    src/divonne/common.c
    src/divonne/FindMinimum.c
    src/divonne/Integrate.c
    src/divonne/Rule.c
    src/divonne/Divonne.c
    src/divonne/Split.c
)
