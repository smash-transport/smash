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

echo "Building SMASH ..."
echo

# TODO: remove it, just have /build
# mkdir -p smash-benchmark
# cd smash-benchmark

mkdir -p build || fail "Failed to create $PWD/build directory."
cd build || fail "Failed to change directory to $PWD/build."

#cmake -DCMAKE_BUILD_perfYPE=Release $smash_root_dir || fail "Failed to configure the build system."
#make -j${n_cores} || fail "Failed to build SMASH."

echo
echo "Collecting System Information ..."


#TODO: remove tail for master (output is cleaned up there)
SMASH_VER_NUM=$(./smash -v | grep -E 'SMASH-' | tail -n1)
SMASH_VER_OUT=$(./smash -v)
HW_INFO=$(lshw -short 2> /dev/null)
UNAME_INFO=$(uname -a)

echo
echo "Running Benchmarks ..."
echo

benchmark_run() {
  setup=$1

  PERF_OUT=$(perf stat -B -r 3 ./smash -i ${configs_dir}/${setup}/config.yaml -d ${configs_dir}/${setup}/decaymodes.txt -p ${configs_dir}/${setup}/particles.txt 2>&1 >/dev/null)

  echo "$PERF_OUT"
}

echo "Started benchmark for collider ..."
coll_perf=$(benchmark_run collider)
echo "$coll_perf" | grep -E "time elapsed"

echo "Started benchmark for box ..."
box_perf=$(benchmark_run box)
echo "$box_perf" | grep -E "time elapsed"

echo "Started benchmark for sphere ..."
sphere_perf=$(benchmark_run sphere)
echo "$sphere_perf" | grep -E "time elapsed"


OUTPUT_FILE_N=bm-results-${SMASH_VER_NUM}.md

# \` are espaced backticks
cat > ../${OUTPUT_FILE_N}<<EOF
# Benchmark Results for $SMASH_VER_NUM

### \`smash -v\`
\`\`\`
$SMASH_VER_OUT}
\`\`\`

## System information

### \`uname -a\`
\`\`\`
$UNAME_INFO
\`\`\`

### \`lshw -short\`

\`\`\`
$HW_INFO
\`\`\`

## Results

Performed with \`perf\` and the \`stat -B\` options.

### Collider Run
\`\`\`
$coll_perf
\`\`\`

### Box Run
\`\`\`
$box_perf
\`\`\`

### Sphere Run
\`\`\`
$sphere_perf
\`\`\`
EOF
echo
echo "Benchmark results are written to $OUTPUT_FILE_N"
