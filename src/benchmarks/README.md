# SMASH Benchmarks

In the current state the benchmarks are very basic. A few common SMASH run
scenarios are tested. More scenarios as well as microbenchmarks
might be added in the future.

## Running the benchmarks

The benchmarks are run by executing the shell script.

    ./benchmark.sh

This first builds SMASH, collects are few system infos and runs the different
SMASH setups. The measurement is done with `perf` and performed a few times to
obtain an average. The measurement is written
to a markdown formatted file named `bm-results-SMASH-VERSION.md`.

## Adding other setups

You might add other common SMASH scenarios by adding the configs to the
respective directory. Additionally, the shell script has to be modified manually
at the moment.
