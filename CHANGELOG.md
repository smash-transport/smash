# Changelog


All notable changes to this project will be documented in this file. The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). This project does not adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html), but it uses versioning inspired by it. Given a version number `X.Y.Z`,

* `X` is incremented for major changes in particular large backwards incompatible updates,
* `Y` is incremented for minor changes like an external pull-request that adds one feature and
* `Z` for the indication of a bug fix or other very small changes that are not backwards incompatible.

The major categories to group changes in this log are:

* `Input / Output` for all, in particular breaking, changes, fixes and additions to the in- and output files.
* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Fixed` for any bug fixes.
* `Removed` for now removed features.

Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.

## [SMASH-2.0.2](https://github.com/smash-transport/smash/compare/SMASH-2.0.1...2.0.2)

### Input / Output
* Added interface with Rivet (the particle-physics MC analysis toolkit) that can process the data of the simulations at runtime and produce in output YODA files with the results of the analysis (tested with Rivet version 3.1.4)
* Refactoring of the HempMC3 output, in asciiv3 format, both for final particles and for the collision history (tested with HepMC3 version 3.2.3)

### Changed
* Enabling the Rivet output requires a compiler supporting c++14 
* c++ standard version flags removed from CMakeLists.txt, they can (or must) be supplied to cmake as command line arguments
* Minimum supported cmake version: 3.1

## [SMASH-2.0.1](https://github.com/smash-transport/smash/compare/SMASH-2.0...2.0.1)

### Input / Output
* Fixed event number counting for intermediate OSCAR output

## [SMASH-2.0](https://github.com/smash-transport/smash/compare/SMASH-1.8...2.0)

**New major version of SMASH**

Major changes:

* Multi-particle interactions: There is infrastructure for 3<->2 and 3<->1 interactions and example processes are implemented.
* The [SMASH hybrid](https://github.com/smash-transport/smash-vhlle-hybrid) is available including SMASH in initial and final state as well as the vHLLE viscous hydrodynamics code.
* Updated to Pythia 8.303 and optimised the function calls to allow SMASH runs at LHC and high RHIC energies.
* No large backwards incompatible updates for this major version.

### Input / Output
* Unify event number counting: The event number counting now starts at 0 for all output formats.
* Added charge chemical potential to the box modus as an option.
* Added option to specify the `Maximum_Cross_Section` that is considered from the config file.
* Clarify naming of empty events output flag: changing `empty` to `scattering_projectile_target`.
* Newly available HepMC3 output.
* Option to scale all cross sections by a global `Cross_Section_Scaling` factor from the config file.
* New `Additional_Elastic_Cross_Section` option to add an additional constant contribution to the elastic cross section in the config file.
* Option to `Only_Warn_For_High_Probability` in case of long production runs with the stochastic criterion.

### Added
* 3-to-1 reactions for mesons via the stochastic collision criterion.
* 3-to-2 reactions for deuterons via the stochastic collision criterion.
* New covariant collision criterion.
* Hadron Gas EoS extended by nQ and muQ.
* Various tests for photons.
* Various tests for the stochastic criterion.
* Various tests for multi-particle reactions.
* Command-line option to enable completely silent output.
* Travis CI check to ensure zero doxygen warning about undocumented instances.
* Added a new Process Type to tag for failed string processes.
* Option to initialize the box with Bose or Fermi distributions.
* More extensive documentation for the different collision criteria.
* A map of Pythia objects to reduce number of Pythia initializations.
* f-resonance to pi-pi-photon Bremsstrahlung cross-sections.
* Tests for nucleon density normalization.

### Fixed
* Consider cross section scaling factor of incoming particles for photon production. This factor was previously neglected, resulting in exploding weights and overestimated photon production.
* Use form factors for binary scattering photons also in the case of Nfrac=1.
* The evaluation of the interaction point in the boxmodus is considering interactions through the wall properly.
* Enforced a small time step if the box modus is used.
* Bug in displaying the total energy per particle.

### Changed
* The Pythia version is increased to 8.303.
* In collisions of unformed particles with equal formation time the outgoing particles now always inherit the smaller scaling factor.
* The default collision criterion changed from "geometric" to "covariant".
* The interaction point for the stochastic criterion is the coordinate of a random incoming particle instead of the middle point. This prevents density artifacts in the center of grid cells.

### Removed
* Integrator1dMonte as it was not used anymore.
* Integrator2d with the GSL Monte-Carlo integration functions. Replaced by the Integrator2dCuhre which performs 2D integrations according to the Cuhre algorithm. The latter was now renamed to Integrator2d.



## [SMASH-1.8.1](https://github.com/smash-transport/smash/compare/SMASH-1.8...SMASH-1.8.1)
Date: 2020-08-13

### Changed
* Improve version determination for Pythia.

## [SMASH-1.8](https://github.com/smash-transport/smash/compare/SMASH-1.7...SMASH-1.8)
Date: 2020-04-07

### Input / Output
* Some column names in the extended ROOT output have changed to be agreement with those in the extended OSCAR output:
  * `coll_per_part` -> `ncoll`
  * `formation_time` -> `form_time`
  * `xsec_factor` -> `xsecfac`
* Add strangeness and baryon number to VTK particles output.
* The `Only_Final` option for the Particles output has the new possibility if `IfNotEmpty` is entered. Then only output is written if the event is not empty i.e. there was a collision between projectile and target to save disk space. The other two options are now `Yes` or `No` (previously `True` or `False`).
* Restructure photon configuration:
  * Add new `Photon` subsection to `Collision_Term` section, with options `2to2_Scatterings: True/False`, `Bremsstrahlung: True/False` and `Fractional_Photons: Nfrac`, where `Nfrac` is an arbitrary integer
  * Add new `Dilepton` subsection to `Collision_Term`, with option `Decays: True/False`
  * From now on, only the format and whether or not it shall be extended, can be set in the Output subsection `Photons`.

### Added
* Photon production from pion-pion bremsstrahlung
* Option to randomize the reaction plane (rotate all particles by a random angle around the z axis).
* Test to verify that cross sections do not depend on particle order.
* Various tests for deformed nuclei
* Various tests for string processes
* Proper CHANGELOG.md file
* Option to set up an equilibration time for the box modus, after which the output is written.
* Most resonance integrals are now cached on disk, reducing the time until the simulation starts.

### Changed
* Extend box examples in user guide to directly run.

### Fixed
* Fix user guide for baryon density dependent symmetry potential.
* ASCII initial conditions output:
  * Bugfix: write PDG ID instead of unique particle ID.
  * Exclude spectators.
* Add lower bound for initial conditions proper time.
* Fix floating point exceptions that were raised if initial conditions output is enabled
* Fix light nuclei not being affected by potentials



## [SMASH-1.7](https://github.com/smash-transport/smash/compare/SMASH-1.6...SMASH-1.7)
Date: 2019-10-14

### Input / Output
- New output content: Initial conditions output in ASCII, Oscar, Binary and ROOT format
- New option to explicitly specify the output times in terms of a list
- Input option `Included_2to2` now has a new possible value `NNbar`. This is a breaking change since this value now has to be in the included list of 2to2 reactions in order to use the `resonances` option for `NNbar_treatment`.
- Input configuration of deformed nuclei has changed. Additional level of `Orientation` provided where the angles or a random rotation can be specified. This change is backwards incompatible to config files for SMASH-1.6 or older.
- Most outputs now specify whether the projectile collided with the target
- ROOT output now includes charge and optionally the extended output
- Added particle ID and number of collisions to VTK particles output
- Thermodynamic output can now also output electric, baryonic and strange currents

### Added
- New stochastic collision criterion (as in A. Lang, H. Babovsky, W. Cassing, U. Mosel, H. G. Reusch, and K. Weber, J. Comp. Phys. 106, 391 (1993)) for 2-to-2 reactions of one particle species with a fixed elastic cross-section.
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
- Fix output time for output at the event end on the command line to always be correct.



## [SMASH-1.6](https://github.com/smash-transport/smash/compare/SMASH-1.5.2...SMASH-1.6)
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



## [SMASH-1.5.2](https://github.com/smash-transport/smash/compare/SMASH-1.5.1...SMASH-1.5.2)
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



## [SMASH-1.5.1](https://github.com/smash-transport/smash/compare/SMASH-1.5...SMASH-1.5.1)
Date: 2018-12-04

### Changed
- Documentation of formation time
- Continuous particle formation for unformed particles from string
fragmentation processes

### Fixed
- Correct generic radius for nuclei with A > 16, except for Au, Pb, U, Cu as they are/were exactly specified



## [SMASH-1.5](https://github.com/smash-transport/smash/releases/tag/SMASH-1.5)
Date: 2018-11-27

**First public version of SMASH**

Known issues:

- Simulations at high &radic;<span style="text-decoration: overline">s</span><sub>NN</sub> above ~30 GeV are slow due to the initialization of Pythia for hard scatterings
- Strangeness production at intermediate energies is off due to the string excitation and fragmentation (only tuned to pp elementary data)
- No Coulomb potential is available
