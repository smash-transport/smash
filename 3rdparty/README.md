## Third-party libraries into SMASH

Few external libraries are included and shipped within SMASH.
Since some of them are still under active development, they are regularly updated.
In general, for different reasons, some adjustments in third-party _CMakeLists.txt_ files has been needed and these modifications have to be maintained whenever an update is performed.
This README file is meant to contain information about the work that has been done on external libraries as well as guidance for future changes.
**It is crucial to follow the advice here while updating any external library.**

### General remarks

:warning: :arrow_right: ADD HERE COMMENT ABOUT `-fPIC` compiler flag, see #846.

When updating either Cuba or YAML,
* the `include_directories` commands in the top-level _CMakeLists.txt_ file needs to be adjusted;
* similarly the `add_subdirectory` commands in _3rdparty/CMakeLists.txt_ must be changed;
* the _cmake/FindSMASH.cmake_ file contains few occurrences that require the same type of adjustment.

Cuba and Einhard libraries use the `-march=native` flag that is not supported by e.g. Appleclang 13.0 compiler on M1 machines.
Hence a patch has been added in [commit d064cac0](https://github.com/smash-transport/smash-devel/commit/d064cac0fa5b68089cf0c0c2ac7c11c358a1c15d) and it has to be maintained.

### YAML

We use [this library](https://github.com/jbeder/yaml-cpp) mainly for SMASH configuration files and we update it from time to time, since it is maintained and active.
If you want to include a new version inside SMASH codebase, you need to do something along the following lines (from within the `3rdparty` folder).
```bash
rm -r yaml-cpp-* 
YAML_VERSION='X.Y.Z' # Put in the right numbers
wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-"${YAML_VERSION}".zip
unzip yaml-cpp-"${YAML_VERSION}".zip
mv yaml-cpp-yaml-cpp-"${YAML_VERSION}" yaml-cpp-"${YAML_VERSION}" && rm yaml-cpp-"${YAML_VERSION}".zip
git status yaml-cpp-"${YAML_VERSION}"
unset -v 'YAML_VERSION'
```
**Be sure not to forget the steps mentioned in the general remarks above!**
Of course, you should try to compile SMASH with the new version of the library and check that everything works fine.
If third-party warnings pop up, have a look to them and see if something can be done about it.
It is in general not desired to pollute compilation with warnings.
For example, for version `0.7.0` the compiler flag `-Wno-shadow` has been added since some warnings otherwise where given and the YAML developers declared them [as fine](https://github.com/jbeder/yaml-cpp/issues/764).

### Cuba

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
**Be sure not to forget the steps mentioned in the general remarks above!**

### Virtest

Although not very active, we check for updates of [this library](https://github.com/mattkretz/virtest) from time to time and, in particular, after every release.
To check for updates, you can e.g. follow the procedure here below, from within the `3rdparty` folder.
```bash
rm -r virtest
git clone --depth 1 git@github.com:mattkretz/virtest.git
rm -r virtest/.git
git status virtest
# If something changed the library evolved from that shipped with SMASH.
```

### Einhard

[This library](https://gitlab.com/Marix/Einhard) seems inactive, but it is sound and it has never given problems.
Due to CMake policies, the minimum required CMake version has been increased from `2.6.2` to `3.0`.
It is planned to leave this library frozen, unless C++ problems arise.
If anything will be changed at some point, **be sure not to forget the steps mentioned in the general remarks above**, as well as using a modern CMake version.



