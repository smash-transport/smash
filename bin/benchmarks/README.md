# SMASH Benchmarks

In the current state, the benchmarks are very basic. A few common SMASH run
scenarios are tested. More scenarios as well as microbenchmarks
might be added in the future.

## Preparation

Before running, prepare a cleaned and (re)compiled `build` directory with the
version you want to benchmark.

## Running the benchmarks

The benchmarks are run by executing the shell script.

    ./benchmark.sh PREPARED_BUILD_DIR

First the benchmark collects system information and then runs the different
SMASH setups. The measurement is done with `perf` and performed a few times to
obtain an average and measure system fluctuations. The measurement is written
to a markdown formatted file named `bm-results-SMASH-VERSION.md`.

## Results from tagged version

For comparison, results form previous tagged version are kept in the `results`
directory. The directory is reserved for benchmarks done during the tagging
procedure.


## Adding other setups

You may add other common SMASH scenarios. First add the configs to the
respective directory and then modify the shell script accordingly.
