# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

## [0.3.0] - 2026-09-12

### Changed

- the first iteration of PrimalProximalHeur is not penalized, so that
  the bound it produces is the one of the original problem: get_lb() /
  get_ub() report it and it is what the accuracy the heuristic is asked
  for is measured against. intUseWarmStartPSol asks for the previous
  behaviour, i.e. that the primal solution of the warm start be the
  proximal center of the first iteration

- the parameters of PrimalProximalHeur named after the base Solver refer
  to the heuristic itself: intMaxIter is the number of its iterations and
  dblMaxTime the time it is given as a whole, out of which each call to
  the inner Solver gets what is left; the inner Solver is given
  intInnerMaxIter and dblInnerRelAcc (was dbl_LDSRelAcc), and the
  configuration files of the Solver of the warm start and of the primal
  recovery are strWarmStartBSC (empty, i.e. no warm start, by default)
  and strRecoveryBSC

- the trace of PrimalProximalHeur is compiled in and silent by default,
  intLogVerb turning it on at runtime, rather than the other way around

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- the warm start of PrimalProximalHeur asked the auxiliary Solver for the
  dual solution whatever happened to it: a relaxation stopped by the time
  limit has none, and the query threw. What the Solver has is asked for,
  and if the duals are not there the multipliers stay where set_Block()
  put them, exactly as with no warm start

- the package configuration file finds the libraries the module links, so that
  a project using the installed module needs nothing more than find_package()

## [0.2.0] - 2025-12-12

### Added 

- [huge] PrimalProximalHeur Lagrangian-based math-heuristic

- put_State() and get_State() methods

- Configuration for get\_[var/dual]\_solution() of the InnerSolver can now be set

### Changed 

- [big] managing of intPushCostToOwner parameter of LagBFunction

- adapted to new un\_any\_count thing

### Fixed 

- added missing parameter initialization

- flawed handling of inner Configuration

- missing actual setting of int\_InnerS\_W*Cfg in set\_par()

- right solver called in get\_dual\_solution()

- avoided static vectors prone to static initialization fiasco 

- a bunch of stupid bugs

## [0.1.3] - 2024-02-28

### Changed 

- adapted to new CMake / makefile organisation

### Fixed 

- exploiting the new "father of LagBFunction" mechanism to make Modification
  from sub-Bloch to reach their original father (instead of UpdateSolver)

- tentative solution for Modification issues (just ignore them)

## [0.1.2] - 2022-06-28

### Fixed 

- locking the Solver inside compute()

- bunch of minor fixes

## [0.1.1] - 2021-05-02

Minor point release to avoid the master branch to become too stale:

- significant improvements in handling Configurations

- fixed a number of issues (get_lb/ub exchanged, computation of solutions)

- Lagrangian variables now initialized with dual variables from the FRowConstraint

## [0.1.0] - 2021-05-02

Initial release

### Added

- Initial release.

[Unreleased]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.3...0.2.0
[0.1.3]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.2...0.1.3
[0.1.2]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/tags/0.1.0
