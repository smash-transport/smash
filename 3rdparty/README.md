# Third-party libraries into SMASH

Few external libraries are included and shipped within SMASH.
Since some of them are still under active development, they are regularly updated.
In general, for different reasons, some adjustments in third-party _CMakeLists.txt_ files has been needed and these modifications have to be maintained whenever an update is performed.
This README file is meant to contain information about the work that has been done on external libraries as well as guidance for future changes.
**It is crucial to follow the advice here while updating any external library.**

## General remarks

All third party libraries are shipped with their requirement about the minimum CMake version and SMASH has one as well.
**Unless the library is needing a larger version of that required by SMASH**, we drop the third library requirement,
```diff
- cmake_minimum_required(VERSION ...)
+ # CMake minimum version inherited from SMASH context
```
and this change must be in place at every update.

When updating either Cuba or YAML,
* the `include_directories` commands in the top-level _CMakeLists.txt_ file needs to be adjusted;
* similarly the `add_subdirectory` and `install` commands in _3rdparty/CMakeLists.txt_ must be changed;
* the _cmake/FindSMASH.cmake_ file contains few occurrences that require the same type of adjustment.

Cuba and Einhard libraries use the `-march=native` flag that is not supported by e.g. Appleclang 13.0 compiler on M1 machines.
Hence, it is important to add this flag using the `add_compiler_flags_if_supported` CMake function.
This is already in place and it has to be maintained.

It is worth adding a last remark about how libraries are to be compiled.
This should not require any maintenance change, although it is not guaranteed to be so in the future.
Since SMASH is going to be shipped as a shared library, we need to offer position independent code.
This is imposed at SMASH CMake top-level for all targets and should propagate automatically for third-party libraries, too.

### Remark about the C++ standard

Each third party library is responsible for itself and, in particular, this includes the choice of the C++ standard, which is not necessary the same of that used by SMASH.
Being each library included as subdirectory, it will have its own variable scope initialized with a copy of the current variable values at the time of the `add_subdirectory()` call.
Therefore, changes to variables like `CMAKE_CXX_FLAGS` do not affect the top-level, SMASH scope.
Using target properties should be preferred, but strictly speaking the obsolete, traditional CMake approach of setting variables instead of properties still works.

### Remark about the CMake `install` target

The installation of SMASH comes along with the frozen third-party libraries frozen in it, so that in this respect using SMASH as external library behaves the same as when using it as software.

If a library provides an installation procedure (i.e. CMake `install` calls), like Einhard, these commands should be removed at every update.
If instead the library, like YAML, offers an installation only when not used in another project, then the CMake files should not be changed, as that installation will not be used by SMASH.
In general, SMASH installation deals with the third-party libraries in the *CMakeLists.txt* file in the ***3rdparty*** folder and those lines should not need changes when updating a library (apart from trivial folder names).

## YAML

We use [this library](https://github.com/jbeder/yaml-cpp) mainly for SMASH configuration files and we update it from time to time, since it is maintained and active.
If you want to include a new version inside SMASH codebase, you need to do something along the following lines (from within the `3rdparty` folder).
```bash
rm -r yaml-cpp-*
YAML_VERSION='X.Y.Z' # Put in the right numbers
wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/"${YAML_VERSION}".zip
unzip "${YAML_VERSION}".zip && rm yaml-cpp-"${YAML_VERSION}".zip
rm -r yaml-cpp-"${YAML_VERSION}"/{test,util}
unset -v 'YAML_VERSION'
```
Removing the `test` and `util` folders is crucial to avoid shipping hundreds of unnecessary files which sum up to some unnecessary MB.
This is connected to the fact we do not build YAML tests or tools (cf. `YAML_CPP_BUILD_TESTS` and `YAML_CPP_BUILD_TOOLS` CMake variables in the _3rdparty/CMakeLists.txt_ file) and it is likely not changing in the future.

**Be sure to drop any CMake version requirement as described in the general remarks above.**

Of course, you should try to compile SMASH with the new version of the library and check that everything works fine.
If third-party warnings pop up, have a look to them and see if something can be done about it.
It is in general not desired to pollute compilation with warnings.
For example, for version `0.7.0` the compiler flag `-Wno-shadow` has been added since some warnings otherwise where given and the YAML developers declared them [as fine](https://github.com/jbeder/yaml-cpp/issues/764) and this was not needed anymore for version `0.8.0`.


## Cuba

Releases and further information are available on [a dedicated webpage](http://www.feynarts.de/cuba/).
The library is maintained and active and we plan to update it from time to time.
Since it is not shipped with CMake code to build it, the needed CMake code has been [added in the past](https://github.com/smash-transport/smash-devel/commit/eef0dd995ced5ff1c54571fa2d296cff58b31739) and all _CMakeLists.txt_ files have to be kept during the update and adjusted as needed.
To update Cuba, proceed along the following lines from within the `3rdparty` folder in a `bash` shell.
```bash
OLD_CUBA_VERSION='X.Y.Z' # Put in the right numbers
NEW_CUBA_VERSION='X.Y.Z' # Put in the right numbers
wget http://www.feynarts.de/cuba/Cuba-"${NEW_CUBA_VERSION}".tar.gz
tar -xf Cuba-"${NEW_CUBA_VERSION}".tar.gz && rm Cuba-"${NEW_CUBA_VERSION}".tar.gz
shopt -s globstar
for f in Cuba-"${OLD_CUBA_VERSION}"/**/CMakeLists.txt; do mv "${f}" "${f/${OLD_CUBA_VERSION}/${NEW_CUBA_VERSION}}"; done
shopt -u globstar
rm -r Cuba-"${OLD_CUBA_VERSION}"
git status Cuba-"${NEW_CUBA_VERSION}"
unset -v 'OLD_CUBA_VERSION' 'NEW_CUBA_VERSION'
```
**Be sure not to forget the steps mentioned in the general remarks above.**


## Virtest

Although not very active, we check for updates of [this library](https://github.com/mattkretz/virtest) from time to time and, in particular, after every release.
To check for updates, you can e.g. follow the procedure here below, from within the `3rdparty` folder.
```bash
rm -r virtest
git clone --depth 1 git@github.com:mattkretz/virtest.git
rm -r virtest/.git
# Apply changes done in commit 97c2b917 to remove CMake minimum version requirement
git status virtest
# If something changed the library evolved from that shipped with SMASH.
```
The changes to drop the CMake minimum requirement of the library can be easily applied by hand to the _3rdparty/virtest/CMakeLists.txt_ file.
Add
```diff
+# CMake minimum version inherited from SMASH context
+
```
at the very top of the file and remove the call to `cmake_minimum_required`, e.g.
```diff
-cmake_minimum_required(VERSION 2.6)
```


## Einhard

[This library](https://gitlab.com/Marix/Einhard) seems inactive, but it is sound and it has never given problems.
Few changes have been done to integrate the library into SMASH.

 * Due to CMake policies, the minimum required CMake version has been implicitly increased as already described.
 * The main _einhard.hpp_ header file was not guarded and header guards have been added.

It is planned to leave this library frozen as untouched as possible, unless C++ problems arise.
If anything is changed at some point in the upstream and such a changes are meant to be pulled into SMASH, **be sure not to forget the steps mentioned in the general remarks above**.
Futhermore, at every update, the `install` lines should be deleted, as we install the libraries differently from the *CMakeLists.txt* file in the ***3rdparty*** library.



