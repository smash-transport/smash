# Contributing to SMASH

If you are reading this, it is probable that you are going to contribute to
SMASH. First of all, we welcome your interest and thank you for your work.
Hopefully, you find what you need in the following. Otherwise feel free to
contact the development team.

Note that any contributions must be licensed under the same terms as SMASH, see
[LICENSE](LICENSE).

## Testing

### Running Tests

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


### Runtime memory checking with valgrind

The SMASH binary memory usage can be checked for the different modi with the the
following cmake targets:

    make memcheck_collider
    make memcheck_box
    make memcheck_sphere

Alternatively, the binary can be checked manually via:

    valgrind -v ./smash

Note: There is known bug with `valgrind-3.11` that leads to an error about an
unrecognized instruction. The memchecks will not run with his version.

## Development Tools

The following tools can be helpful for development:
- clang-format = 6.0
- doxygen >= 1.8.4
- valgrind
- cpplint
- cppcheck
- codespell


### Installing binaries as a user

If you don't have administrator privileges on the machine you are using, you can
still install software locally. Just copy the binary (in this example
`clang-format`) to a local folder and update your path:

    mkdir ~/bin
    cp ./clang-format ~/bin
    echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
    source ~/.bashrc

After this, you can just copy executable files to `~/bin` to install them. This
also works for other exectuables like cpplint. You might have to set them to be
executable with `chmod u+x ~/bin/my-binary`.


### Installing clang-format

clang-format is a part of the clang compiler. You can download the most recent
binaries here:

   http://releases.llvm.org/download.html

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
"Installing binaries as a user" above.


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


## Code Documentation

We use [doxygen](http://doxygen.org) for generating documentation from our code.
The online version of the code documentation is found
[here](https://fias.uni-frankfurt.de/~smash/extra/doc/).


### How to build docs with doxygen yourself

You need to have doxygen installed. Then just call:

    make doc

and you'll find `doc/html/index.html` in your build directory. Open e.g. with

    firefox doc/html/index.html

Additionally, there are two more targets that can be used to test the
completeness of the documentation:

    make undocumented
    make undocumented_count

The first target outputs all doxygen warnings about missing documentation and
the second one only counts (and outputs) the number of warnings. Both are
building the doxygen documentation only for completely documented entities, but
the main purpose of both is that all warnings are displayed when running.


#### Building the User Guide

Call

    make user

to obtain the files in `doc/user/index.html`


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


### How to write good doxygen comments

Doxygen is very flexible in the [comments it accepts for documentation
generation](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html). In
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
with `///` or wraped with `/** ... */`. For all the connection to code and
some special layout commands, refer to the [Special Commands in
Doxygen](http://www.stack.nl/~dimitri/doxygen/manual/commands.html).

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
* Use square brackets for units `[...]`
* Check if the documentation is properly parsed by doxygen with `make doc` (see
  below)
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

In order to refer to a paper inside a doxygen comment, the `iref` command should
be used:

```
/** ... this function implements ... as described in \iref{XXX}. */
int fun int(x);
```

Here, `XXX` should be the BibTex key for the paper from Inspire e.g.
`Weil:2016zrk`. In order to find it, search for the paper on
http://inspirehep.net and then click on 'BibTex', which will show the complete
BibTex entry (you only need the key, which is in the first line). Doxygen will
automatically translate `\iref{XXX}` into a link to the paper on Inspire.

After adding a new reference, you should run the script `doc/get_bibtex.sh`,
which will update the file
[smash.bib](https://fias.uni-frankfurt.de/pm/projects/smash/repository/revisions/master/entry/doc/smash.bib)
by fetching the BibTex entries of all `iref` references from Inspire.


### User Guide

The User Guide will be written in the code base, i.e., documentation of
configuration options are described where they are used. Comments that are in
normal doxygen format do not appear in the User Guide. Instead only multi-line
comments of the form `/*!\Userguide ... */` will be used. Example:

``` cpp
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
Doxygen tree. Currently, the cmake script processes all `src/include/*.h` and
`src/*.cc` files.

### Markdown Documents

The markdown documents included in the repository follow the Markdown dialect of Github that is specified [here](https://github.github.com/gfm/).

## Coding Rules

* [Google Naming & Formatting
Rules](https://fias.uni-frankfurt.de/~smash/extra/code_guidelines/cppnaming.xml) -
  The rules enforce consistent coding style throughout our code base.

It can be very instructive to study the [CERT Coding
Standards](https://www.securecoding.cert.org/confluence/display/seccode/CERT+Coding+Standards)
as well. Secure coding is not as important for SMASH as for most other projects.
However, secure coding basically is a strategy for avoiding program errors, and
as such it is interesting for every C++ project.

Our Naming & Formatting Rules follow the Google Styleguide. We keep a copy of
the relevant sections to be able to make small adjustments and hide the C++
Coding Guidelines sections.

Note that we use `clang-format` for formatting (see [README](README.md) for installation,
usage explained below), so you can just use that instead of having to worry
about formatting by hand in accordance with the style guide.

A notable naming rule is that class member variables have to end in an
underscore (`_`).


### Specific Naming Conventions for SMASH

* for particles going into or coming out of an interaction we use:
  * `incoming_particles`
  * `outgoing_particles`
* for lists of `ParticleData` objects (type `ParticleList`) we use
  * `particle_list`
* for lists of `PdgCode` objects we use
  * `pdg_list`
* to enumerate particles in any context we use:
  * `particle_a, particle_b, particle_c, ...`

### Exceptions to the Rules

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


### How to add an enum for a new configuration value

If a new option is added to the `config.yaml`, it should be done in the
following way:

* Add `enum class` to `forwarddeclarations.h`.
* Add `Configuration::Value` cast operator overload.
* Document possible strings in `config.yaml` and User Guide.


### Code formatting with `clang-format`

All code has to be formatted by running `clang-format`. This automatically
formats the code in SMASH correctly. Use the helper script in SMASH's /bin
directory to format the source code:

    ./clang-format-helper -p

Review and commit changes afterwards. clang-format does changes that 
don't look good, you can disable it locally using comments like this:


    // clang-format off
    ...
    // clang-format on


### Coding style and static code analysis

The SMASH source code can be checked via:

    make cpplint 2>&1 | grep -v 'Done processing'

`cpplint` checks the formatting of the code and is also part of the unit tests.
`cppcheck` is a static code analyzer that can be run additionally, but yields a
lot of false positives:

    make cppcheck


## Profiling and benchmarking

This section discusses tools that can be used to measure the SMASH's
performance.

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
