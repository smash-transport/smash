# Changelog


All notable changes to this project will be documented in this file. The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). This project does not adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html), but it uses versioning inspired by it. Given a version number `X.Y.Z`,

* `X` is incremented for major changes in particular large backwards incompatible updates,
* `Y` is incremented for minor changes like an external pull-request that adds one feature and
* `Z` for the indication of a bug fix or other very small changes that are not backwards incompatible.

The major categories to group changes in this log are:

* `Input / Output` for all, in particular breaking, changes, fixes and additions to the in- and output files;
* `Added` for new features;
* `Changed` for changes in existing functionality;
* `Fixed` for any bug fixes;
* `Removed` for now removed features.

Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.

## Unreleased

## SMASH-3.3rc
Date:

### Added
* Added spin 4-vectors (Pauli-Lubanski) to particles including spin interactions. At the moment, spin interactions are treated in elastic scatterings and in inelastic scatterings through `Σ* → Λπ` formation. The spin vector is by default initialized unpolarized.
* Heavy flavour particles can be treated as perturbative particles in inelastic interactions.

### Changed
* Update to Pythia8.316.
* Cross sections for baryon-antibaryon interactions are now shifted w.r.t. the effective mass threshold in AQM rescaling.
* The `ParticleData` class now initially sets the formation time and the time of the last collision to `NAN` and overwrites them with the real values the first time the objects are used. Previously, these values were initially set to `0.0`, which could be misleading if the time evolution started at negative times.
* The procedure for dynamic initial conditions was updated, fluidized hadrons are no longer removed from the evolution, in order to provide energy density to the surrounding particles, being internally named as "core" particles. Core and corona (non-core) particles can only interact elastically.
* The automatic lattice heuristic determination uses a minimum lattice size and has been disabled in `List` modus.
* ⚠️ The `backpropagate_to_the_same_time` method is now a private member of `ListModus`. Its functionality was moved to the free function `backpropagate_straight_line`, which *does not set* the formation time or cross section scaling of particles.
* The `only_res` boolean parameter of `find_final_actions` was removed in the base class `ActionsFindersInterface` and all its derived action classes.

### Fixed
* Properly take into account potentials at threshold in resonances decay channels selection.
* Add better root finding for momentum dependent potentials.

### Input
* The `Sphere` modus now can include thermal heavy flavor particles, via the `Modi: Sphere: Heavy_Flavor_Multiplier` key, a value that multiplies their partial density.
* Back to back jets were introduced to the `Sphere` modus with the  `Modi: Sphere: Jet: Back_To_Back` and `Modi: Sphere: Jet: Back_To_Back_Separation` keys.
* The initial `Sphere` modus velocity profile can now be tuned with  `Modi: Sphere: Add_Radial_Velocity_Exponent`.
* `ListModus` now accepts `Optional_Quantities`, extra to the default values used in Oscar2013.
* Added the `Proper_Time_Scaling` key in the `Initial_Conditions` section under `Modi: Collider:` to scale the switching proper time when using constant tau initial conditions.
* Added a new key `Collision_Term: Spin_Interactions` to enable spin interactions with values `On` and `Off`.
* If `Spin_Interactions` is `On`, the `List` modus requires four additional columns representing the components of the spin 4-vector (s_0, s_1, s_2, s_3) and these shall be specified in the `Optional_Quantities` list.

### Output
* ⚠️ The previous `<content>_binary.bin` output files are now called `<content>_custom.bin`.
* ⚠️ The previous `ASCII` option for `Output: Initial_Conditions: Format` key was renamed to `For_vHLLE`, and the corresponding file to `SMASH_IC_For_vHLLE.dat`. **This breaks the workflows for hybrid models using SMASH as initial conditions!**
* `Dileptons`, `Photons`, and `Initial_Conditions` now accept the `ASCII` and `Binary` custom format, and require setting the desired `Quantities`.
* Added `tau`, `eta_s`, `mt`, `y_rap`, `spin0`, `spinx`, `spiny`, `spinz`, and `perturbative_weight` as custom ASCII/binary quantity.
* Improve writing performance of both ASCII and binary output formats by buffering information event by event.
* The `formation_time` and `time_last_collision` TBranches in the ROOT output have been moved from the extended to the ordinary (default) output.

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-3.2.2...SMASH-3.3)


## SMASH-3.2.2
Date: 2025-09-19

### Fixed
* Fix too strict energy-momentum conservation leading to failures at high energies

### Changed
* The `almost_equal` and `almost_equal_physics` functions to compare floating point numbers have been improved to perform a relative comparison à la Knuth (the latter considers numbers below a given threshold as equal)

### Added
* Offer conversion to string for all input keys so that, when using SMASH as library, setting `Configuration` key values does not require to implement specializations of the `YAML::convert` class template for key types not convertible to strings

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-3.2.1...SMASH-3.2.2)


## SMASH-3.2.1
Date: 2025-07-08

### Changed
* Upgraded to Pythia 8.315, which contains a bug-fix (see [here](https://pythia.org/pdfdoc/aextra.pdf) for further information).
  This hot-fix is meant to provide a SMASH version including the updated Pythia version.

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-3.2...SMASH-3.2.1)


## SMASH-3.2
Date: 2025-02-25

### Added
* Strangeness production by resonances tuned to elementary production cross sections using a genetic algorithm
* Added charmed particles with corresponding interactions (AQM for elastic scatterings and resonance production and decay for inelastic interactions)
* Added alpha-clustered oxygen nuclei that can be configured in a new `Alpha_Clustered` section in the `Projectile`/`Target` section
* Implemented dynamic initial conditions for hydrodynamics, available under `Modi: Collider: Initial_Conditions: Type: "Dynamic"`
* Infrastructure for functional tests, which requires at least Python3.3, and first functional test

### Input
* Moved the `Orientation` section from the `Deformed` to the `Projectile`/`Target` section and changed the default value for `Theta` from `pi/2` to `0`
* The `Include_Weak_And_EM_Decays_At_The_End` was renamed to `Ignore_Minimum_Decay_Width_For_Decays_At_The_End`

### Output
* Implemented a new `ASCII` value for the `Format` key of `Output: Particles` and `Output: Collisions`  which creates a *.dat* file containing columns based on a list of user-input `Quantities` (the existing `Oscar1999` and `Oscar2013` formats are simply a convenient alias for specifying `ASCII` and a predefined list of `Quantities`)
* ⚠️ The ensemble number is now included in the output next to the event number - **this is a potentially breaking change for analysis software, especially if using the binary output**; refer to the documentation of each output type for further details
* ⚠️ The former `Binary` format has been renamed to `Oscar2013_bin` and the `Binary` format now exists as a user-customizable output w.r.t. the output quantities (it works analogously to the `ASCII` format)
* Bump binary output format version from 9 to 10
* Some binary output files have been renamed to better track their content
* Removed the information on the number of particles from the initial conditions output
* Initial conditions, HepMC and Rivet outputs disabled when SMASH is run with multiple parallel ensembles

### Changed
* Upgraded to Pythia 8.312
* The interface of the `Configuration` class has been totally changed making it use `Key` objects instead of low-level strings
* The `List` modus now validates particles from all events at the very beginning and aborts if more than 2 particles are at the same identical position in any event

### Fixed
* Fix a small error and added a more precise documentation for the orientation section for deformed nuclei
* Make few keys really mandatory in code as it was already stated in the documentation
* Oscar2013 format for Initial Conditions no longer show the number of particles
* Fix the number of each event in output files when using multiple parallel ensembles

### Deprecated
* Some physics-related input keys in `Output: Initial_Conditions` were deprecated and duplicated under `Modi: Collider: Initial_Conditions`

### Removed
* The deprecated `Version` key in the configuration file is not accepted anymore
* The `Include_Weak_And_EM_Decays_At_The_End` key has been renamed and hence is not accepted anymore

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-3.1...SMASH-3.2)


## SMASH-3.1
Date: 2024-02-26

### Added
* Added momentum dependence of the Skyrme potential using the cold nuclear matter approximation
* Option `Decay_Initial_Particles` which blocks the decay of resonances created at initialization
* The `List` and `ListBox` modi can now be used specifying a single particles input file using the `Filename` key
* Command-line option to run smash with `-n` to avoid caching integrals in the ***tabulations*** folder on the disk
* ⚠️  The `Total_Cross_Section_Strategy` key determines how the total cross sections for collision finding are evaluated, whether from summing the partial contributions or using a parametrization; previously the default behavior was the former strategy (now accessible with the value `BottomUp`), but it has changed to the latter (accessible with `TopDownMeasured`)
* ⚠️  Option `Pseudoresonance` which adds *ad hoc* inelastic 2→1 processes close to the energy regime around the transition to strings, where few resonances have been measured; by default, this is done for interactions where at least one of the incoming particles is unstable, and the heaviest possible resonance producible from the incoming pair is chosen
* New `Use_Potentials_Outside_Lattice` option in the config file to specify where potentials are applied

### Fixed
* Fix bug in the time step calculation of SMASH output in events successive to the first one when the `End_Time` is not a multiple of the `Output_Interval`
* SMASH now aborts with a user-friendly message if any input key among `End_Time`, `Delta_Time` or `Output_Interval` is zero
* Fix bug in thermodynamics lattice ASCII output when both ASCII and binary formats are requested
* Fix bug about SMASH getting stuck with some particular values of `End_Time` and `Delta_Time`
* Fix value of particles mass in SMASH IC ASCII output (the pole mass was printed instead of the real mass)
* SMASH does not crash if `rho_eckart` is not specified in the `Thermodynamics` output
* SMASH does not crash if the same output format is repeated
* Fix bug when `Power_Particle_Formation` is explicitly specified
* SMASH can now be built with GNU compiler `v13` that requires explicit inclusion of `cstdint` header
* Fix bug when using SMASH as a library and particles should be removed

### Changed
* Use updated gold radius of 6.55 fm instead of 6.38 fm and a diffusiveness of 0.523 fm instead of 0.535 fm
* Upgraded to Pythia 8.310
* SMASH installation procedure was improved and allows now to use SMASH as library relying on installation folder only
* CMake module to find SMASH in an external project offers now a better approach based on the new installation procedure
* Make `Clock` class comparison operators refer to internally represented time and not to `current_time()` return value
* Add strangeness to all extended output formats in the particle line

### Deprecated
* Locate SMASH in external projects using the `SMASH_DIR` and not relying on an installation

### Removed
* Some public setters of the `ListModus` class have been removed

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-3.0...SMASH-3.1)


## SMASH-3.0
Date: 2023-04-27

### New major version of SMASH

Summary of major changes from `SMASH-2.2.1`:
* ⚠️ The `master` branch has been renamed to `main`
* Updated to Pythia 8.309 and introduced automatic usage of Monash tune for high beam energies to match multiplicities up to LHC energies
* SMASH has now more functionality when used as an external library
* Optimized SMASH at low energies achieving a ~25% speedup
* Added new validation mechanism of input configuration file to prevent misuse and improve corresponding documentation
* Upgraded to C++17 standard and get rid of Boost dependency
* Minimum versions requirements of prerequisites have increased and officially supported systems were made explicit
* Backward incompatible changes have occurred, mainly on the input configuration side (few keys changed and the internal parsing class interface was substantially rewritten)

### Added
* New CMake build configuration `MinSizeRel` to optimize in executable size
* New examples for usage of SMASH as a library: Full SMASH wrapper and rate equations
* New validation mechanism of SMASH input (configuration file and `-c` command line option)
* Add an input key to employ the Monash tune for all Pythia processes (this is automatically enabled in `Collider` mode for `Sqrtsnn` ≥ 200 GeV, unless explicitly disabled)
* Cross section parameters for the transition to strings are now optional inputs
* Add a new functionality to add or remove particles when using SMASH as a library

### Fixed
* Fix bug in collider modus when setting the total energy per nucleon of _each_ beam
* The box modus can only be used with the fixed time step mode
* Avoid integer overflow to silently happen in grid construction
* Correct tau lepton mass in particles list files
* Fix code behavior in list modi setup to be coherent with given warnings in case of input discrepancies
* Fix falsely applied smearing factor in thermodynamic output of the charge currents `j_BQS` in the case of no smearing

### Changed
* The `master` branch has been renamed to `main`
* SMASH makes now use of C++17 standard and minimum compiler requirements have changed
* Upgraded to Pythia 8.309
* Minimum supported CMake version is now version 3.16
* Minimum Doxygen version to build documentation is now 1.9
* Only 64-bit UNIX-like operating systems (e.g. Linux, MacOS) are officially supported
* Only GCC, Clang and Apple clang compilers are officially supported
* Default build configuration is now `Release` (instead of `RelWithDebInfo`)
* Use CMake default compiler flags for `RelWithDebInfo` build configuration, hence having `-O2` optimization level instead of `-03`
* Changed interface of the `Configuration` class, which is now not copyable and offers more explicit methods, making it harder to be misused
* The `Version` configuration key is not used anymore and has been deprecated
* SMASH now aborts if any content output format in the configuration file is invalid or absent
* `Format: ["None"]` can be used in the configuration file to suppress any content output
* Added two quantities (baryon number, strangeness) in the particle line in ASCII initial conditions output (SMASH_IC.dat)
* Clebsch-Gordan coefficients are now tabulated resulting into a performance increase, especially at low energies
* Made the input `Jet_PDG` key in the `Jet` section in sphere and box modus mandatory, when the `Jet` section is provided
* The `Lattice` section in the input file cannot be anymore empty; the key `Automatic: True` must be now used to automatically generate the lattice and `Automatic: False` is needed if the geometry is fully specified
* All input mass values that need to match the SMASH internally used ones are now validated
* Updated 2 → 2 angular distribution parametrization for high energies

### Removed
* Boost is not used anymore, since the C++17 standard library is sufficient now
* Removed `Rest_Frame_Density_Derivatives_Mode` input key which was not used in the code

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.2.1...SMASH-3.0)


## SMASH-2.2.1
Date: 2022-05-18

### Fixed
* Properly boost in radial direction in sphere mode

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.2...SMASH-2.2.1)


## SMASH-2.2
Date: 2022-05-10

### Input / Output
* Added option `Minimum_Nonempty_Ensembles` to define a minimum number of ensembles in which interactions took place
* Added `HepMC_treeroot` output
* Added option to `Include_Weak_And_EM_Decays_At_The_End`

### Added
* Light nuclei (A=3) production via 4-to-2 reactions
* Possibility to perform non-strong decays at the end of the calculation
* Possibility to impose transverse momentum or rapidity cut when extracting initial conditions for hydrodynamics
* Option to add a velocity field of the form `u_r = u_0 * r / R` in radial direction to the sphere mode
* Interface functions to initialize SMASH conveniently (see new `library.h`)
* Extend benchmarks to include collisions at high energy to test strings

### Fixed
* The `--version` option of `smash` now returns a version also when SMASH is downloaded as source package without VCS information
* The saturation density of the woods-saxon distribution is calculated such that its integral returns the number of nucleons in the nucleus
* Initial conditions output now includes unformed particles

### Changed
* Minimum supported CMake version is now version 3.9
* Renamed `USE_ROOT`, `USE_HEPMC` and `USE_RIVET` CMake options to `TRY_USE_ROOT`, `TRY_USE_HEPMC` and `TRY_USE_RIVET`, respectively
* Upgraded to Pythia 8.307
* Renamed version macro to `SMASH_VERSION` (from `MAJOR_VERSION`)
* Pythia installation directory can also be specified with `-DPYTHIA_ROOT_DIR` and by the environment variables `PYTHIA_ROOT_DIR` or `PYTHIA8`, with `/usr` as default
* Big endian architectures are not officially supported and compilation on them is, by default, disabled
* Specific parameters for Xenon nuclei are used instead of default values

### Removed
* A possibly available system installation of YAML library is not considered anymore
* The singularity container definition has been dropped, since a singularity container can be obtained from the shipped docker

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.1.4...SMASH-2.2)


## SMASH-2.1.4
Date: 2022-02-24

### Fixed
* Properly reset for each event whether projectile and target interacted

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.1.3...SMASH-2.1.4)


## SMASH-2.1.3
Date: 2022-02-15

### Changed
* Include unformed particles in initial conditions output
* Add definition of endianness macro to FindSMASH.cmake module

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.1.2...SMASH-2.1.3)

### Deprecated
 * Output option `HepMC` to enable HepMC ASCII output (equivalent to the newer `HepMC_asciiv3`)


## SMASH-2.1.2
Date: 2022-02-04

### Changed
* Reflect new Cuba version in FindSMASH.cmake

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.1.1...SMASH-2.1.2)


## SMASH-2.1.1
Date: 2022-01-31

### Changed
* Fixed typo in PgdCode class

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.1...SMASH-2.1.1)


## SMASH-2.1
Date: 2021-12-21

### Input / Output
* Added option to specify a `Fixed_Min_Cell_Length` to control the grid size for the stochastic criterion
* Computation of thermodynamic quantities optionally restricted to participants only using the option `Only_Participants`

### Added
* 5-to-2 reactions for NNbar annihilations via the stochastic collision criterion
* 2-to-5 reactions for NNbar annihilations to 5 pions
* Thermodynamic lattice output
* Parallel ensembles
* VDF potentials
* New smearing modes: discrete + triangular
* New option: finite difference derivatives for calculation of density gradients on lattice
* New option: sample particles according to quantum statistics
* Non-relativistic electromagnetic potentials in a simple form meant for low collision energies
* New option: start a box with an external particle list
* Support of Apple silicon M1 chips
* Optional N(1520) Dalitz decay with constant form factor
* Pre-built docker container on Github

### Changed
* Evaluation of failed string processes. BBbar pairs are now forced to annihilate
* Updated Dockerfile and Singularity definition file (matching pre-built container on Github)

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.0.2...SMASH-2.1)


## SMASH-2.0.2
Date: 2021-06-23

### Input / Output
* Added interface with Rivet (the particle-physics MC analysis toolkit) that can process the data of the simulations at runtime and produce YODA output files with the analysis results
* Refactoring of the HempMC3 output now available for final particles and for the collision history in asciiv3 format

### Changed
* Minimum supported cmake version is now version 3.1
* Continuous integration now uses Github Actions

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.0.1...SMASH-2.0.2)


## SMASH-2.0.1
Date: 2020-12-21

### Input / Output
* Fixed event number counting for intermediate OSCAR output

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-2.0...SMASH-2.0.1)


## SMASH-2.0
Date: 2020-12-17

### New major version of SMASH

Major changes:
* Multi-particle interactions: There is infrastructure for 3 ↔ 2 and 3 ↔ 1 interactions and example processes are implemented
* The [SMASH hybrid](https://github.com/smash-transport/smash-vhlle-hybrid) is available including SMASH in initial and final state as well as the vHLLE viscous hydrodynamics code
* Updated to Pythia 8.303 and optimized the function calls to allow SMASH runs at LHC and high RHIC energies
* No large backwards incompatible updates for this major version

### Input / Output
* Unify event number counting: The event number counting now starts at 0 for all output formats
* Added charge chemical potential to the box modus as an option
* Added option to specify the `Maximum_Cross_Section` that is considered from the config file
* Clarify naming of empty events output flag: changing `empty` to `scattering_projectile_target`
* Newly available HepMC3 output
* Option to scale all cross sections by a global `Cross_Section_Scaling` factor from the config file
* New `Additional_Elastic_Cross_Section` option to add an additional constant contribution to the elastic cross section in the config file
* Option to `Only_Warn_For_High_Probability` in case of long production runs with the stochastic criterion

### Added
* 3-to-1 reactions for mesons via the stochastic collision criterion
* 3-to-2 reactions for deuterons via the stochastic collision criterion
* Parallel ensembles technique to run simulations with mean field potentials faster
* New covariant collision criterion
* Hadron Gas EoS extended by nQ and muQ
* Various tests for photons
* Various tests for the stochastic criterion
* Various tests for multi-particle reactions
* Command-line option to enable completely silent output
* Travis CI check to ensure zero doxygen warning about undocumented instances
* Added a new Process Type to tag for failed string processes
* Option to initialize the box with Bose or Fermi distributions
* More extensive documentation for the different collision criteria
* A map of Pythia objects to reduce number of Pythia initializations
* f-resonance to pi-pi-photon Bremsstrahlung cross-sections
* Tests for nucleon density normalization

### Fixed
* Consider cross section scaling factor of incoming particles for photon production. This factor was previously neglected, resulting in exploding weights and overestimated photon production
* Use form factors for binary scattering photons also in the case of Nfrac=1
* The evaluation of the interaction point in the boxmodus is considering interactions through the wall properly
* Enforced a small time step if the box modus is used
* Bug in displaying the total energy per particle

### Changed
* The Pythia version is increased to 8.303
* In collisions of unformed particles with equal formation time the outgoing particles now always inherit the smaller scaling factor
* The default collision criterion changed from "geometric" to "covariant"
* The interaction point for the stochastic criterion is the coordinate of a random incoming particle instead of the middle point. This prevents density artifacts in the center of grid cells

### Removed
* Integrator1dMonte as it was not used anymore
* Integrator2d with the GSL Monte-Carlo integration functions. Replaced by the Integrator2dCuhre which performs 2D integrations according to the Cuhre algorithm. The latter was now renamed to Integrator2d

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.8...SMASH-2.0)


## SMASH-1.8.1
Date: 2020-08-13

### Changed
* Improve version determination for Pythia

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.8...SMASH-1.8.1)


## SMASH-1.8
Date: 2020-04-07

### Input / Output
* Some column names in the extended ROOT output have changed to be agreement with those in the extended OSCAR output:
  * `coll_per_part` → `ncoll`
  * `formation_time` → `form_time`
  * `xsec_factor` → `xsecfac`
* Add strangeness and baryon number to VTK particles output
* The `Only_Final` option for the Particles output has the new possibility if `IfNotEmpty` is entered. Then only output is written if the event is not empty i.e. there was a collision between projectile and target to save disk space. The other two options are now `Yes` or `No` (previously `True` or `False`)
* Restructure photon configuration:
  * Add new `Photon` subsection to `Collision_Term` section, with options `2to2_Scatterings: True/False`, `Bremsstrahlung: True/False` and `Fractional_Photons: Nfrac`, where `Nfrac` is an arbitrary integer
  * Add new `Dilepton` subsection to `Collision_Term`, with option `Decays: True/False`
  * From now on, only the format and whether or not it shall be extended, can be set in the Output subsection `Photons`

### Added
* Photon production from pion-pion bremsstrahlung
* Option to randomize the reaction plane (rotate all particles by a random angle around the z axis)
* Test to verify that cross sections do not depend on particle order
* Various tests for deformed nuclei
* Various tests for string processes
* Proper CHANGELOG.md file
* Option to set up an equilibration time for the box modus, after which the output is written
* Most resonance integrals are now cached on disk, reducing the time until the simulation starts

### Changed
* Extend box examples in user guide to directly run

### Fixed
* Fix user guide for baryon density dependent symmetry potential
* ASCII initial conditions output:
  * Bugfix: write PDG ID instead of unique particle ID
  * Exclude spectators
* Add lower bound for initial conditions proper time
* Fix floating point exceptions that were raised if initial conditions output is enabled
* Fix light nuclei not being affected by potentials

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.7...SMASH-1.8)


## SMASH-1.7
Date: 2019-10-14

### Input / Output
- New output content: Initial conditions output in ASCII, Oscar, Binary and ROOT format
- New option to explicitly specify the output times in terms of a list
- Input option `Included_2to2` now has a new possible value `NNbar`. This is a breaking change since this value now has to be in the included list of 2to2 reactions in order to use the `resonances` option for `NNbar_treatment`
- Input configuration of deformed nuclei has changed. Additional level of `Orientation` provided where the angles or a random rotation can be specified. This change is backwards incompatible to config files for SMASH-1.6 or older
- Most outputs now specify whether the projectile collided with the target
- ROOT output now includes charge and optionally the extended output
- Added particle ID and number of collisions to VTK particles output
- Thermodynamic output can now also output electric, baryonic and strange currents

### Added
- New stochastic collision criterion (as in A. Lang, H. Babovsky, W. Cassing, U. Mosel, H. G. Reusch, and K. Weber, J. Comp. Phys. 106, 391 (1993)) for 2-to-2 reactions of one particle species with a fixed elastic cross-section
- Initial conditions for hydrodynamic simulations can be extracted and are provided in the initial conditions output. The proper time of the produced hypersurface can, if desired, be specified in the configuration file
- Fermi motion for deformed nuclei
- Random rotation of custom nuclei
- Automatic deformation parameters for ruthenium and zirconium
- Option to add net-baryon density dependence to symmetry potential
- Option to add high momentum particle into a box calculation

### Changed
- Changed the NN transition region to 3.5-4.5 GeV to better describe the experimental data
- Box now has minimum length given by the grid size
- The final-state cross-section output provided by `smash -S pdgID1,pdgID2` no longer includes resonances

### Fixed
- Collision finding with frozen Fermi motion: Use beam momentum for collision finding if nuclei have not interacted
- Fix output time for output at the event end on the command line to always be correct

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.6...SMASH-1.7)


## SMASH-1.6
Date: 2019-05-24

### In- and Output
- Fix thermodynamic VTK output of the Landau velocity

### Added
- Particle properties updated to use PDG 2018 data (this includes 5 new N* and a new Delta*)
- Lund fragmentation function for leading baryons in soft non-diffractive string processes
- Support for customized nuclei: optionally reading in initial nucleon positions in collider modus
- New density type for charge/isospin and possibility to output the charge currents without smearing
- Option to place a single high-momentum particle in the center of the sphere modus (~jet)
- SMASH compiles now also with GCC >= 9, Clang >= 7

### Changed
- No debug output in the default build type, which results in much better performance

### Fixed
- Documentation of formation time
- Correct declaration of Woods-Saxon distribution for deformed nuclei

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.5.2...SMASH-1.6)


## SMASH-1.5.2
Date: 2019-03-07

**Version for JETSCAPE run**

### Added
- New tests for density and potential
- New density type for charge/isospin
- Added option to propagate a high momentum hadron through the sphere

### Changed
- Exact Pythia version is now required
- Improve third party infrastructure
- Default build type amended with one flag
- Unify functions for EM and hadronic widths of resonances

### Fixed
- Bug fixes related to formation of particles including
    - Correcting the boost
    - Adjust resonance decay times
    - Add unformed particles to grid for collision finding, if they form during the time step

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.5.1...SMASH-1.5.2)


## SMASH-1.5.1
Date: 2018-12-04

### Changed
- Documentation of formation time
- Continuous particle formation for unformed particles from string
fragmentation processes

### Fixed
- Correct generic radius for nuclei with A > 16, except for Au, Pb, U, Cu as they are/were exactly specified

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-1.5...SMASH-1.5.1)


## SMASH-1.5
Date: 2018-11-27

**[First public version of SMASH](https://github.com/smash-transport/smash/releases/tag/SMASH-1.5)**

Known issues:

- Simulations at high &radic;<span style="text-decoration: overline">s</span><sub>NN</sub> above ~30 GeV are slow due to the initialization of Pythia for hard scatterings
- Strangeness production at intermediate energies is off due to the string excitation and fragmentation (only tuned to pp elementary data)
- No Coulomb potential is available
