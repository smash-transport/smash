# SMASH Benchmarks

In the current state, the benchmarks are very basic. A few common SMASH run
scenarios are tested. More scenarios as well as microbenchmarks
might be added in the future.

## Preparation

Before running, prepare a cleaned and (re)compiled `build` directory with the
version you want to benchmark.

## Running the benchmarks

The benchmarks are run by executing the shell script.
```console
./benchmark.sh PREPARED_BUILD_DIR
```

First the benchmark collects system information and then runs the different
SMASH setups. The measurement is done with `perf` and performed a few times to
obtain an average and measure system fluctuations. The measurement is written
to a markdown formatted file named `bm-results-SMASH-VERSION.md`.

## Results from tagged version

For comparison, results from previous tagged versions are attached to the
corresponding github releases.

## Comparing benchmarks among versions

A utility script to easily compare the benchmark output is available and can
be used specifying the result files of the version to be compared, e.g.
```console
./compare_benchmarks.bash  bm-results-SMASH-2.1rc.md  bm-results-SMASH-2.2rc.md
```

## Adding other setups

You may add other common SMASH scenarios. First add the configs to the
respective directory and then modify the shell script accordingly.
