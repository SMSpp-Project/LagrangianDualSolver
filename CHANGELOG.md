# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

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

### Fixed 

- get_ub() of a minimization problem published the bound the inner Solver
  has on the Lagrangian Dual, which lies on the same side of the optimum as
  the dual value itself: whenever the inner Solver proved its own
  optimality, LagrangianDualSolver claimed to have solved the Block, which
  is false as soon as there is a duality gap. It now publishes only the
  bound the relaxation gives, and compute() returns kLowPrecision unless
  the two bounds close, so that what is promised is in the return code

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

[Unreleased]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.2.0...develop
[0.2.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.3...0.2.0
[0.1.3]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.2...0.1.3
[0.1.2]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/tags/0.1.0
