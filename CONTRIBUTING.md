# Contributing

If you are reading this, it is probable that you are going to contribute to
SMASH. First of all, we welcome your interest and thank you for your work.
Hopefully, you find what you need in the following. Otherwise feel free to
contact the development team by
[reporting issues](https://github.com/smash-transport/smash/issues)
or contact us [by email](mailto:elfner@itp.uni-frankfurt.de).

Note that any contributions must be licensed under the same terms as SMASH, see
our [LICENSE](https://github.com/smash-transport/smash/blob/main/LICENSE).
In particular, it has been decided that changes to C++ or CMake
code should be reflected in the years of the copyright notice. Therefore, be sure
that the present year is added there in files that you are going to commit.

As an external contributor, go to https://github.com/smash-transport/smash, fork
the repository, work on a topic branch and then create a pull request. Details
on this workflow can be found e.g.
[here](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project).

### Table of content

1. [The most important rule](#consistency)
2. [Testing](#testing)
3. [Choosing a build type](#build-type)
4. [Development tools](#development-tools)
5. [Code documentation](#documentation)
6. [Coding rules](#coding-rules)
7. [General policies](#general-policies)
8. [Profiling and benchmarking](#profiling)
9. [SMASH containers](#containers)

---

<a id="consistency"></a>

## Be consistent!

The first and main guideline beyond all the following ones in this file is:
**Be consistent with anything that is already existing in the codebase.**
This is valid in the most general sense, ranging from style, to notation,
to nomenclature and whatever else you might notice.



<a id="testing"></a>

## Testing

### Unit tests

To run the various unit tests, use the following:

    make test

Another way to do this is to use the CMake test runner:

    ctest

This has the advantage that it can also be used for running tests in parallel on
a multicore machine, e.g. via

    ctest -j4

(on a quad-core machine).

If a test crashes, there might be some leftover in the `test_output` folder,
causing the test to always fail when run again. To fix this problem, just remove
the folder.

### Functional tests

The functional tests require Python3.3 to create a virtual environment where modules can be imported. They are disabled by default, but can be enabled before building the project with

    cmake -DENABLE_FUNCTIONAL_TESTS=ON ..

They can be run collectively with

    ctest -R functional

Notice that there is no executable created for them, and so they cannot be run by calling `make` as for the unit tests.

### Runtime memory checking with valgrind

The SMASH binary memory usage can be checked for the different modi with the the
following cmake targets:

    make memcheck_collider
    make memcheck_box
    make memcheck_sphere

Alternatively, the binary can be checked manually via:

    valgrind -v ./smash

Note: There is known bug with `valgrind-3.11` that leads to an error about an
unrecognized instruction. The memchecks will not run with this version.



<a id="build-type"></a>

## Choosing a build type

There are different build types available, which compile the SMASH code for
different situations.

To build a binary with debug symbols (which is useful to get stacktraces in case
of errors), with support for `TRACE` and `DEBUG` logging output (which costs a
lot of performance) and with internal consistency checks (which can be useful
for finding and fixing bugs, but are expensive), use the `Debug` build type:

    cmake .. -DCMAKE_BUILD_TYPE=Debug

For a build with optimizations, without debug symbols, without support for
`TRACE` or `DEBUG` logging and without internal consistency checks, use the
`Release` build type:

    cmake .. -DCMAKE_BUILD_TYPE=Release

For a profiling build, use the `Profiling` build type:

    cmake .. -DCMAKE_BUILD_TYPE=Profiling

The default build type is 'RelWithDebInfo', which is the same as the `Release`
build type, except that it enables debug info for better stacktraces when SMASH
crashes. The debug info makes the binaries larger, but only has a marginal
performance impact.


### Enhancing build verbosity

To find cmake build errors (best debugged with full compiler output) use:

    VERBOSE=1 make

### Find deprecated C library function calls

If a call to a C function without explicit namespace qualification is made,
a warning appears (if enabled) during compilation, if some code tries to use a C
library function call. The only correct way to write C++ code is to use the C++
interface by specfing the correct std namespace.

This option is OFF by default and is checked from time to time. Enable it as
follows

    cmake .. -DDEPRECATE_C_FNS=ON
    make

Note: Do not use this as a default, since it increases compilation time. Use
clang for compilation here, since gcc runs into an internal compiler error,
warnings are however still reported.



<a id="development-tools"></a>

## Development tools

The following tools can be helpful for development:
- clang-format = 13.0.x
- doxygen >= 1.9 (not 1.9.4 or 1.9.5)
- valgrind
- cpplint
- cppcheck
- codespell

**NOTE:** The formatting with the above mentioned `clang-format` version is enforced
at each merge to the `main` branch.


### Installing clang-format

clang-format is a part of the clang compiler. The usage of version 13.0.x is
enforced (x stands for any number).
You can download the binaries [here](http://releases.llvm.org/download.html).

Make sure to pick a pre-built binary for your system. For example, for Ubuntu
you could run:

    $ lsb_release -a
    No LSB modules are available.
    Distributor ID: Ubuntu
    Description:    Ubuntu 16.04.2 LTS
    Release:        16.04
    Codename:       xenial

This tells you to download "Clang for x86_64 Ubuntu 16.04". (You might have to
look for an older version to get pre-built binaries.)

It is sufficient to unpack the archive with `tar xf` and to copy only the
binary you need (`clang-format` in the `bin` folder of the archive), see
"Installing binaries as a user" below.


**NOTE:** The name of the clang-format executable should not have any suffix,
like e.g. `clang-format-13`. If this is the case, a possible workaround is to
create a symbolic link to this exectuable called `clang-format` in a directory
included in the environment variable PATH.
For example:

```
    mkdir -p "${HOME}/bin"
    ln -s "$(which clang-format-13)" "${HOME}/bin/clang-format"
    export "PATH=${HOME}/bin:${PATH}"
```

Of course, if the user has writing permissions, the symbolic link can be also
directly created in the same directory in which the executable is and, if this
is already in the system PATH (for example `/usr/bin`), no further actions are
needed to use it as a normal command.


### Installing cpplint

We use cpplint to enforce some of our style guide lines as part of our tests.
You can install it like this:

    pip install --user cpplint

You might have to add `~/.local/bin` to your `$PATH`, see "Installing binaries
as a user".


### Installing cppcheck

You can use cppcheck to find some problems in the code, but beware it has quite
a few false positives. Download and compile the latest version:

    git clone git://github.com/danmar/cppcheck.git
    cd cppcheck
    make

You can then copy it to your local binary folder, see "Installing binaries
as a user".


### Installing codespell

If you want to check the spelling in comments, try codespell. You can install
it like this:

    pip install --user codespell

It is the same as installing cpplint.


### Installing binaries as a user

If you do not have administrator privileges on the machine you are using, you can
still install software locally. Just copy the binary (in this example
`clang-format`) to a local folder and update your path, e.g.:

    mkdir ~/bin
    cp ./clang-format ~/bin
    echo 'export PATH=${PATH}:${HOME}/bin' >> ~/.bashrc
    source ~/.bashrc

After this, you can just copy executable files to `~/bin` to install them. This
also works for other executables like cpplint. You might have to set them to be
executable with `chmod u+x ~/bin/my-binary`.



<a id="documentation"></a>

## Code documentation

We use [doxygen](http://doxygen.org) for generating documentation from our code.
The online version of the code documentation is found
[here](https://theory.gsi.de/~smash/doc/current/).


### How to build documentation with Doxygen yourself

You need to have Doxygen installed. Then, from your build directory, just call:

    make doc

and you'll find `doc/html/index.html` in your build directory. Open e.g. with

    firefox doc/html/index.html

Additionally, there are two more targets that can be used to test the
completeness of the documentation:

    make undocumented
    make undocumented_count

The first target outputs all Doxygen warnings about missing documentation and
the second one only counts (and outputs) the number of warnings. Both are
building the Doxygen documentation only for completely documented entities, but
the main purpose of both is that all warnings are displayed when running.


#### Building the user guide

In the build directory, run

    make user

to obtain the files in _**doc/user/**_. Open _doc/user/index.html_ in your
favorite browser.


### Doxygen pages and their ordering

New pages can be created using the `\page` Doxygen command, which takes an anchor
name as first argument and the page title as second argument. It is possible to
use the `\page` command with the **same** anchor name in order to add content
to the same page from different places. Therefore, the second argument is optional
and should be used only at the first `\page` occurrence (to avoid inconsistencies).
In general, this is not trivial, but in SMASH we collect pages "declarations" in
the _doc/index.dox_ file, which is meant to enforce pages ordering (Doxygen picks
up pages in the order it encounters them). So, in this file titles should be
specified, while elsewhere not.

To make searching for documentation pages easier, a `doxypage_` prefix for their
anchors has been decided to be used and you should do the same. For instance, if
you create a new documentation page, you will declare it with title in the
_doc/index.dox_ file where you want it to appear as

    \page doxypage_my_new_page My wonderful title

and then add content to it simply by using

    /**
     * \page doxypage_my_new_page
     *
     * [...]
     */

at the desired place in the codebase.

Sub-pages are pages themselves, but Doxygen classifies pages as sub-pages as soon
as it encounters the `\subpage` command in a given `\page`. At that point the
anchor given to the `\subpage` command tells Doxygen which page should be marked
as sub-page of which other. Also sub-pages ordering is fixed in the _doc/index.dox_
file and you may refer to it for further information. In short, remember to use
the `\subpage` command **only* in the index file and the `\ref` command elsewhere.


### What to document in the code

Code documentation has two important purposes:

* Documenting how interfaces are supposed to be used. Doxygen creates all the
  boilerplate for this task by parsing the class inheritance and function
  signatures into nice HTML pages. Via special comments this can be completed to
  full API documentation.
* Documenting why things are as they are. Often the how code works is more
  obvious than why it was done this way and not differently. This information
  can be very useful to understand design choices and follow along original
  ideas.


### How to write good Doxygen comments

Doxygen is very flexible in the [comments it accepts for documentation
generation](https://www.doxygen.nl/manual/docblocks.html). In
general, it uses both comments and the source code itself to generate the
documentation. In general documentation should look like this:

``` cpp
/**
 * \brief Brief description of the class
 * Full description of the class.
 */
class Something {
 public:
  /**
   * <Explain what this function does>
   * \param[in] x <Explain input parameter>
   * \param[out] y <Explain input parameter that gets modified within fun>
   * \return <Explain what the function returns>
   * \throw ExceptionClass <The conditions under which the function is
   * expected to throw the named exception>
   */
  int fun(int x);
};
```

Doxygen has a lot of commands for markup. Most importantly it uses special
comment formatting. Comments meant for the documentation are either prefixed
with `///` or wrapped with `/** ... */`. For all the connection to code and
some special layout commands, refer to the [Special Commands in
Doxygen](https://doxygen.nl/manual/commands.html).

#### Rules

The code documentation follows a few rules concerning the formatting and the
question what should be documented:

* The general rule is that all functions need to be documented.
* All parameters (if important also the parameter ranges), the return value and
  exceptions that are thrown by the function  need to be documented.
* Use `\param[in]` for input parameters; `\param[out]` for input parameters that
  get modified within the function; `\return` for what is returned; `\throw` for
  exception that are thrown
* Template parameters are documented with `\tparam`
* Syntax for inline math is `\f$ <math expression> \f$`
* Syntax for display math is `\f[ <math expression> \f]`
* Syntax for generic LaTeX equation environments, `\f{<env-name>} <math expression> \f}`
* Use `\mathbf{x}` instead of `\vec{x}` for vectors in formulae and use
  `\boldsymbol{\alpha}` for greek letters or symbols like `\nabla`
* Use square brackets for units `[...]`
* Check if the documentation is properly parsed by doxygen with `make doc` (see below)
* Doxygen comment formatting:
  * 1-line
```
    /// <comment>
```
  * multi-line (pay attention to alignment)
```
    /**
     * <comments>
     */
```
* normal comment formatting:
  * 1-line
```
    // <comment>
```
  * multi-line (more compact)
```
    /* <comment line 1>
     * <comment line n> */
```


### How to insert references to publications

In order to refer to a paper inside a doxygen comment, the `\iref` command should
be used:

```cpp
/** ... this function implements ... as described in \iref{XXX}. */
int fun int(x);
```

Here, `XXX` should be the BibTex key for the paper from Inspire e.g.
`Weil:2016zrk`. In order to find it, search for the paper on
http://inspirehep.net and then click on 'BibTex', which will show the complete
BibTex entry (you only need the key, which is in the first line). Doxygen will
automatically translate `\iref{XXX}` into a link to the paper on Inspire.

After adding a new reference, you should run the script `doc/get_bibtex.sh`,
which will update the file `/doc/inspire.bib` by fetching the BibTex entries of
all `\iref` references from Inspire. It also reports references that are not
found on Inspire.

References that are not contained in the Inspire database can be handled as
follows: A corresponding BibTex entry should be put into `doc/non_inspire.bib`
manually. It can then be referenced via the `\cite` command:

```cpp
/** ... this function implements ... as described in \cite XXX. */
int fun int(x);
```

### User guide

The User Guide will be written in the code base, i.e., documentation of
configuration options are described where they are used. Comments that are in
normal doxygen format do not appear in the User Guide. Instead only multi-line
comments of the form `/*!\Userguide ... */` will be used. Example:

```cpp
/*!\Userguide
 * \if user
 * This text ONLY appears in the User Guide (useful for a headline that the
   normal documentation already has in the preceeding section)
 * \endif
 * \ifnot user
 * This text will appear only in the Developer Docs
 * \endif
 *
 * Text that appears in both the User Guide and the Developer Docs.
 */
```

The workflow is that `doc/CMakeLists.txt` extracts all `/*!\Userguide ... */`
Sections into `doc/userguide.dox` which is then included into the User Guide
Doxygen tree. Currently, Cmake processes all files specified in the
`doc/UserInputFiles.cmake` in that given order.


### Markdown documents

The markdown documents included in the repository follow the Markdown dialect of Github that is specified [here](https://github.github.com/gfm/).


### Changelog

Major code changes are tracked in the [CHANGELOG.md](CHANGELOG.md) file. This
file is meant to inform the users, but also the other developers. The structure
of the sections is explained at the beginning of the file itself. Every version
has its own section with the date it was tagged. Starting from the bottom with
the section for the first public release. On top of these sections an
_Unreleased_ section is kept. The file is meant to be kept up-to-date together
with the code changes itself, so that at all times the code status is reflected.
Most important are (breaking) changes to our in- or output file. Here, all
details have to be specified. For all other code updates only major changes like
new features or bug fixes have to be given.



<a id="coding-rules"></a>

## Coding rules

* [Google Naming & Formatting
Rules](https://theory.gsi.de/~smash/extra/code_guidelines/cppnaming.xml) -
  The rules enforce consistent coding style throughout our code base.

Our Naming & Formatting Rules follow the Google Styleguide. We keep a copy of
the relevant sections to be able to make small adjustments and hide the C++
Coding Guidelines sections.

Note that we use `clang-format` for formatting (installation explained above,
usage explained below), so that you do not have to worry about formatting by
hand to be in accordance with the style guide.

A notable naming rule is that class member variables have to end in an
underscore (`_`).

### Specific naming conventions for SMASH

* for particles going into or coming out of an interaction we use:
  * `incoming_particles`
  * `outgoing_particles`
* for lists of `ParticleData` objects (type `ParticleList`) we use
  * `particle_list`
* for lists of `PdgCode` objects we use
  * `pdg_list`
* to enumerate particles in any context we use:
  * `particle_a, particle_b, particle_c, ...`

### Exceptions to the rules

1. Ensure that there is only one statement per line.

```cpp
if (bla) {
  foo();
}
```

Don't put bla and foo on the same line to better show what the code is doing.

2. Use CamelCase for enum variants.

```cpp
enum class MyEnum {
  FirstVariant,
  SecondVariant,
  ThirdVariant
}
```


### How to add a new input key for SMASH configuration

If a new key is added to the `config.yaml`, it should be done in the
following way:

* If the key has not a C++ built-in type, add a new `enum class` to
  _forwarddeclarations.h_ for the new key type.
* If the key requires a new `enum`, add `Configuration::Value` cast operator overload.
* Declare the new key in the _input_keys.h_ file according to the
  instructions you find in the documentation of the `InputKeys` class.
  Doing so you will also add the needed description for the User Guide.
* If the key requires a new `enum`, add a `to_string(...)` overload to convert the `enum` to a `std::string`.


### Code formatting with used utilities

All C++ code has to be formatted by running [`clang-format`](https://releases.llvm.org/download.html),
(version `13.0.0`) while CMake code requires [`cmake-format`](https://github.com/cheshirekow/cmake_format)
(version `0.6.13`) to be run. Python scripts are formatted using [`autopep8`](https://github.com/hhatto/autopep8)
(version `2.3.1`), instead. These programs automatically format the code in SMASH correctly.
Use the helper script in SMASH's _**bin**_ directory to format the source code via

    ./codebase-format-helper.bash C++ -p
    ./codebase-format-helper.bash CMake -p
    ./codebase-format-helper.bash Python -p

or by simply using

    ./codebase-format-helper.bash -p

to format all languages at once. Review and commit changes afterwards.
You can also use the `-t` option to test whether the code is correctly
formatted (the script has also a `-h` option that you can check out).
`clang-format` does changes that don't look good, you can disable it
locally using comments like this:


    // clang-format off
    ...
    // clang-format on


### Coding style and static code analysis

The SMASH source coding style can be checked via:

    make cpplint 2>&1 | grep -v 'Done processing'

`cpplint` checks the formatting of the code and is also part of the unit tests.

`cppcheck` is a static code analyzer that can be run additionally, but yields a
lot of false positives:

    make cppcheck

The `make` targets will be created by CMake stage only if given versions of the
commands are installed and found, namely version 1.6.0 for `cpplint` and version
2.8 for `cppcheck`.


### Floating-point precision

All floating point numbers are represented using doubles.


### Not-a-number

There is a built-in type `smash::smash_NaN` which represents a templated [STL
`quiet_NaN()`](https://en.cppreference.com/w/cpp/types/numeric_limits/quiet_NaN)
that should be used in the codebase instead of the STL low-level one.


### Guideline to include header files

As a guideline, try to include only those header files in the class which are directly being used by
the file.



<a id="general-policies"></a>

## General policies

### Input and output compatibility

In general, input and output interfaces should be backwards compatible, when
introducing changes. If there are backwards incompatible changes that affect the
YAML input or the binary output, these should be mentioned in the CHANGELOG file. For example, input keys that got deprecated or removed should be listed.

The release notes need to include a prominent mention of **all** changes,
backwards incompatible or not. In particular, newly introduced config parameters
have to be mentioned.

### Third party codes

In general, the usage of third party codes is discouraged. If there is a scientific necessity
or a major performance gain or time saving by using third party libraries, they can be linked
to SMASH. For the common ones, a description on how to install them in the README file
is sufficient whereas for the less common ones, including them in the **_3rdparty_** folder and
shipping them with SMASH is the better solution. Of course, this involves ensuring the proper
copyright.

Some third-party libraries are already shipped within SMASH. More information about them and,
in particular, instruction about how to update them can be found in the README file in the
**_3rdparty_** folder.



<a id="profiling"></a>

## Profiling and benchmarking

This section discusses tools that can be used to measure the
performance.

### Benchmarks

A basic benchmark script for linux machines is included in the `bin` directory.
It runs different common setups of SMASH and measures them with `perf` (see
below). Usage instructions can be found in the corresponding
[README](bin/benchmarks/README.md).

### Nanobenchmarking

In the codebase some nanobenchmarking code has been written in order to let the
developer inspect critical parts of the code, e.g. counting CPU cycles.
This is in general unused and it is by default excluded from compilation.
In order to get access to it and be able to use it, you need to toggle the
correspondent CMake option passing `-DENABLE_NANOBENCHMARKING=ON` command line
option to `cmake`. Please note that nanobenchmarking is not supported on ARM
architectures like Apple machines with M1 chips.

### GPROF

You can tell cmake to create a build for profiling with the `Profiling` build
type:


    cmake -DCMAKE_BUILD_TYPE=Profiling ..


This will compile the smash code like `Release` mode, but with the `-pg` flag
to instrument the code and create a `gmon.out` file whenever you run a binary.
You can look at the `gmon.out` information with `gprof`.

    # first run smash to create the gmon.out file
    ./smash
    # now run gprof to see the profile information
    gprof smash|less


### Perf

A different method for profiling SMASH uses performance counters and the `perf`
tool on Linux. `perf` uses a feature of the CPU (and OS support) to count
performance relevant events in the normal execution of the program. Therefore,
SMASH should execute as "normal as possible", but with meta information about
the program to improve the reporting capabilities of `perf`. To get useful
results compile SMASH in `RelWithDebInfo` mode, i.e. via


    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..


Now you can use `perf` to collect information for subsequent investigation
(similar to `gprof`) or have it collect global stats directly and output them
after execution. Here are some of the common uses:


    # Collect default information about the execution (cycles, IPC, stalls,
    # branching, ...)
    perf stat -B ./smash


    # Collect detailed information about a run (collecting 'cycles')
    perf record --call-graph dwarf ./smash
    perf report --stdio


`perf` also supports to collect different information, but this can be very CPU
specific and also hard to interpret. Therefore, I recommend to focus only on
`cycles` (for `record`) and the default `stat` output.
A short overview what to make of the `stat` numbers:

* `instructions`:
This is the total number of instructions that the CPU had to execute to run the
program. In general, less is better; but the actual run time is determined by how
many instructions the CPU can execute in parallel (current CPUs can issue up to
4 instructions per cycle). You can see this ratio right behind the
`instructions` value, such as: `2,384,060,126 instructions  # 1.59 insns per
cycle`. The bigger that number the better. But it is really hard to judge the
efficiency of a program by this number. (After all there could be unnecessarily
long instruction sequences that have a high IPC, while a shorter sequence would
be faster overall, but with a lower IPC.)

* `branches` and `branch-misses`
Branches are places where the CPU cannot be certain which code needs to be
executed next. But instead of waiting for the decision the CPU predicts the jump
address and starts executing. If it then determines it miss-predicted it has to
roll back its work and start again from the correct address. Therefore
miss-predictions are costly and should be reduced to a minimum. The branch
miss-prediction ratio is shown in the `perf` output.

* `page-faults`
A page fault happens when a memory address is accessed that has not yet been
allocated in physical memory by the operating system. In that case the TLB
(translation lookaside buffer) entry cannot be determined by the CPU and the OS
has to interrupt the execution of the program to allocate the page and fill it
in the page table for use in the TLB. Page-faults must happen after memory
allocations (unless `malloc` is able to reuse previously deallocated memory)
and are therefore an indicator for "irresponsible" memory allocations.


### Callgrind

Valgrind includes a tool that can profile your code, which should compiled with
debug symbols and optimization turned on. This can be easily achieved in SMASH
by specifying the `-DCMAKE_BUILD_TYPE=RelWithDebInfo` CMake option when setting
up the project. Once compiled SMASH in this mode, run

    valgrind --tool=callgrind ./smash

from the ***build*** directory. Note that this tool is **great and accurate**,
but it will make the execution of your code **extremely slow**. For one SMASH
event with the default configuration file, the execution time will pass from few
dozens seconds to the realm of (tens of) minutes, roughly speaking. Once
terminated, the run generates a file called _callgrind.out.X_, where _X_ usually
is the process ID. Use the `kcachegrind` tool to read this file. It will give
you a graphical analysis of the profiling output with results like which lines
cost how much. Alternatively, you can use `gprof2dot`, which is a more general
tool to visualize the output of different profilers (see below).

It is worth remarking that using the `Profiling` build for this type of measurement
will make calls to `mcount` appear in the calling graph, but this is an artefact
of the profiling procedure, which might even affect measured performance and, hence,
should be avoided.


### gprof2dot

This is a very nice way to visualize profilers output. The software is written in
Python and [open-source](https://github.com/jrfonseca/gprof2dot). You can install
it via `pip`, e.g. via

    pip install --user gprof2dot

for a user-only installation. Then you are ready to use it. Check out the [README
examples](https://github.com/jrfonseca/gprof2dot#examples) as quick-start. To
produce a graph out of the Valgrind output, use something like

    gprof2dot --strip -f callgrind callgrind.out.X | dot -Tsvg -o output.svg

where you need to replace `X` by the proper number and you can choose a better
name for the produced SVG output file. The `--strip` option will

> strip function parameters, template parameters, and
> const modifiers from demangled C++ function names

and it is encouraged to be used to obtain a more readable result with SMASH.


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



<a id="containers"></a>

## SMASH containers

Containers for SMASH are built and shipped next to each new public release.
Instructions about how to build them can be found [here](containers/README.md).
