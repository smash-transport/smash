#!/bin/sh

fail() {
  echo "$*" >&2
  exit 1
}

# TODO:
######### "USER" INPUT FOR NOW #################
smash_root_dir="../../.."
configs_dir="${smash_root_dir}/src/benchmarks/configs"
n_cores=8
##############################################

echo
echo "Building SMASH ..."
echo

# TODO: remove it, just have /build
# mkdir -p smash-benchmark
# cd smash-benchmark

mkdir -p build || fail "Failed to create $PWD/build directory."
cd build || fail "Failed to change directory to $PWD/build."

cmake -DCMAKE_BUILD_TYPE=Release $smash_root_dir || fail "Failed to configure the build system."
make -j${n_cores} || fail "Failed to build SMASH."

echo
echo "Running Benchmarks:"
echo

benchmark_run() {
  setup=$1

  exec_time=$( { TIMEFORMAT="%R"; time sh -c "./smash -i ${configs_dir}/${setup}/config.yaml -d ${configs_dir}/decaymodes.txt -p ${configs_dir}/particles.txt &> /dev/null"; } 2>&1 )  # Captures time only. SMASH muted completely.

  echo "$exec_time"
}

echo "Started benchmark for collider ..."
coll_t=$(benchmark_run collider)

echo "Started benchmark for box ..."
box_t=$(benchmark_run box)

echo "Started benchmark for sphere ..."
sphe_t=$(benchmark_run sphere)

echo
echo "Benchmark Result Times ..."
echo

echo "Collider =    ${coll_t} sec"
echo "Box      =    ${box_t} sec"
echo "Sphere   =    ${sphe_t} sec"

cat > bm_results.txt <<EOF
Collider time:    ${coll_t}s
Box time:         ${coll_b}s
Sphere time:      ${sphe_t}s
EOF
echo "\nResults are also written to bm_results.txt in /build."
