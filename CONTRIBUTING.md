# Contributing

## Code Documentation

We use "doxygen":http://doxygen.org for generating documentation from our code.

### What to document in the code

Code documentation has two important purposes:
* Documenting how interfaces are supposed to be used. Doxygen creates all the boilerplate for this task by parsing the class inheritance and function signatures into nice HTML pages. Via special comments this can be completed to full API documentation.
* Documenting why things are as they are. Often the how code works is more obvious than why it was done this way and not differently. This information can be very useful to understand design choices and follow along original ideas.

### How to write good doxygen comments

Doxygen is very flexible in the "comments it accepts for documentation generation":http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html. In general documentation should look like this:

    /**
     * \brief Brief description of the class
     * Full description of the class.
     */
    class Something {
     public:
      /**
       * <Explain what this function does>
       * \param x <Explain input parameter>
       * \return <Explain what the function returns>
       * \throws ExceptionClass <The conditions under which the function is expected to throw the named exception>
       */
      int fun(int x);
    };


Doxygen has a lot of commands for markup. In general you can work with "Markdown":http://www.stack.nl/~dimitri/doxygen/manual/markdown.html (known from Wiki syntax).
But for all the connection to code and some special layout commands, refer to the "Special Commands in Doxygen":http://www.stack.nl/~dimitri/doxygen/manual/commands.html.

### How to insert references to publications

In order to refer to a paper inside a doxygen comment, the @iref@ command should be used:

    /** ... this function implements ... as described in \iref{XXX}. */
    int fun int(x);


Here, XXX should be the BibTex key for the paper from Inspire. In order to find it, search for the paper on http://inspirehep.net and then click on 'BibTex', which will show the complete BibTex entry (you only need the key, which is in the first line). Doxygen will automatically translate @\iref{XXX}@ into a link to the paper on Inspire.

After adding a new reference, you should run the script doc/get_bibtex.sh, which will update the file "smash.bib":https://fias.uni-frankfurt.de/pm/projects/smash/repository/revisions/master/entry/doc/smash.bib by fetching the BibTex entries of all @iref@ references from Inspire.

### How to build docs with doxygen yourself

You need to have doxygen version 1.8.4 installed. Then just call:

    make doc

and you'll find doc/html/index.html in your build directory.

Additionally, there are two more targets that can be used to test the completeness of the documentation:

    make undocumented
    make undocumented_count

Both are building the doxygen documentation only for completely documented entities, but the main purpose of both is that all warnings are enabled. So, the first target outputs all doxygen warnings about missing documentation and the second one only counts (and outputs) the number of warnings.

## User Guide

The User Guide will be written in the code base, i.e., documentation of configuration options are described where they are used. Comments that are in normal doxygen format do not appear in the User Guide. Instead only multi-line comments of the form @/*!\Userguide ... */@ will be used. Example:

    /*!\Userguide
     * \if user
     * This text ONLY appears in the User Guide (useful for a headline that the normal documentation already has in the preceeding section)
     * \endif
     * \ifnot user
     * This text will appear only in the Developer Docs
     * \endif
     *
     * Text that appears in both the User Guide and the Developer Docs.
     */

The workflow is that @doc/CMakeLists.txt@ extracts all @/*!\Userguide ... */@ Sections into doc/userguide.dox which is then included into the User Guide Doxygen tree. Currently, the cmake script processes all @src/include/*.h@ and @src/*.cc@ files.

### Creating the User Guide

Call

    make user

to obtain the files in @doc/user/index.html@

## Code formatting

You can run clang-format to automatically format all code in SMASH correctly
using the following command:

    for i in src/*.cc src/include/*.h src/tests/*.cc src/tests/*.h; do clang-format -i $i; done;

If clang-format does changes that don't look good, you can disable it locally
using comments like this:

    // clang-format off
    ...
    // clang-format on


## Choosing a build type

There are different build types available, which compile the SMASH code for
different situations.

For a debug build (with all warnings, but no optimization) use:

      cmake .. -DCMAKE_BUILD_TYPE=Debug

For a release build (full optimization and no warnings) use:

      cmake .. -DCMAKE_BUILD_TYPE=Release

For a profiling build use:

      cmake .. -DCMAKE_BUILD_TYPE=Profiling

The default build type is 'RelWithDebInfo', which provides both optimization
and debug info.


## Enhancing build verbosity

To find cmake build errors (best debugged with full compiler output) use:

      make VERBOSE=1


## Runtime memory checking with valgrind

The SMASH binary is regularly checked by:

      valgrind -v ./smash


## Static code analysis & coding style

The SMASH source code can be checked via:

      make cpplint 2>&1 | grep -v 'Done processing'
      make cppcheck



## Profiling / Benchmarking

### GPROF

You can tell cmake to create a build for profiling with the \c Profiling build
type:

    cmake -DCMAKE_BUILD_TYPE=Profiling ..

This will compile the smash code like \c Release mode, but with the \c -pg flag
to instrument the code and create a \c gmon.out file whenever you run a binary.
You can look at the \c gmon.out information with gprof.

    # first run smash to create the gmon.out file
    ./smash
    # now run gprof to see the profile information
    gprof smash|less

### Perf

A different method for profiling SMASH uses performance counters and the \c perf
tool on Linux. \c perf uses a feature of the CPU (and OS support) to count
performance relevant events in the normal execution of the program. Therefore,
SMASH should execute as "normal as possible", but with meta information about
the program to improve the reporting capabilities of \c perf. To get useful
results compile SMASH in \c RelWithDebInfo mode, i.e. via

    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..

Now you can use \c perf to collect information for subsequent investigation
(similar to \c gprof) or have it collect global stats directly and output them
after execution. Here are some of the common uses:

    # Collect default information about the execution (cycles, IPC, stalls,
    # branching, ...)
    perf stat -B ./smash

    # Collect detailed information about a run (collecting 'cycles')
    perf record --call-graph dwarf ./smash
    perf report --stdio

\c perf also supports to collect different information, but this can be very CPU
specific and also hard to interpret. Therefore, I recommend to focus only on \c
cycles (for \c record) and the default \c stat output.
A short overview what to make of the \c stat numbers:

* \c instructions:
This is the total number of instructions that the CPU had to execute to run the
program. In general, less is better; but the actual run time is determined by how
many instructions the CPU can execute in parallel (current CPUs can issue up to
4 instructions per cycle). You can see this ratio right behind the \c
instructions value, such as: `2,384,060,126 instructions  # 1.59 insns per
cycle`. The bigger that number the better. But it is really hard to judge the
efficiency of a program by this number. (After all there could be unnecessarily
long instruction sequences that have a high IPC, while a shorter sequence would
be faster overall, but with a lower IPC.)

* \c branches and \c branch-misses
Branches are places where the CPU cannot be certain which code needs to be
executed next. But instead of waiting for the decision the CPU predicts the jump
address and starts executing. If it then determines it miss-predicted it has to
roll back its work and start again from the correct address. Therefore
miss-predictions are costly and should be reduced to a minimum. The branch
miss-prediction ratio is shown in the \c perf output.

* \c page-faults
A page fault happens when a memory address is accessed that has not yet been
allocated in physical memory by the operating system. In that case the TLB
(translation lookaside buffer) entry cannot be determined by the CPU and the OS
has to interrupt the execution of the program to allocate the page and fill it
in the page table for use in the TLB. Page-faults must happen after memory
allocations (unless \c malloc is able to reuse previously deallocated memory)
and are therefore an indicator for "irresponsible" memory allocations.

### Flame Graphs

Flame graphs are a useful way to visualize the call graph output of a profiler.
To generate them, do:

    git clone https://github.com/brendangregg/FlameGraph
    export FGPATH=$(pwd)/FlameGraph

This gives you a collection of Perl scripts for flame graph generation.
Different profilers are supported. In the following example we are using perf on
Linux:

    perf record -F 99 -g ./smash
    perf script | $FGPATH/stackcollapse-perf.pl > out.perf-folded
    $FGPATH/flamegraph.pl out.perf-folded > perf.svg
    firefox perf.svg

If you get mangled names of the functions, try this as a second line

    perf script | c++filt | $FGPATH/stackcollapse-perf.pl > out.perf-folded

The x axis represents the total duration that the corresponding stack frame
lived for. The order and the colors are arbitrary and optimized for readability.
The y axis represents the position on the stack.
