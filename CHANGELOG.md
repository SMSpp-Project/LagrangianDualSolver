# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the list of the sub-Block this Solver is told to skip, which
  `get_excluded_blocks()` gives, is handed to the inner Solver as it is: any
  sub-Block excluded here is a sub-Block of the Lagrangian dual or of one of
  its descendants, so the inner Solver has to skip it as well

- `intRecursive`, with which the decomposition does not stop at the children
  of the Block: a child having the shape the Block must have, i.e., no
  Variable and no Objective of its own and sub-Block of its own, is
  decomposed in turn, its linking Constraint relaxed together with those of
  the Block and its own children taking its place as components. The descent
  goes on as deep as the shape holds, so that the components are the leaves
  of the decomposable part of the tree and the multipliers those of every
  level

- `PrimalProximalHeur` gives a feasible solution when its relaxed Constraint
  tie copies of a decision, x_a - x_b = 0, as the non-anticipativity ones
  of a two-stage problem do: the copies are fixed to their mean, rounded if
  integer, and the components, independent then, are solved alone with the
  Solver of `strRecoveryBSC` by `intRecoveryThreads` threads

### Changed

- the check that the Block has no Variable of its own asks it for its groups,
  and the dictionaries of the relaxed Constraint are filled one run at a
  time, the vectors of `boost::any` they used to read not being there any
  more

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

### Fixed

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- a change to a relaxed constraint, such as a coefficient rewritten by the
  scaling of a unit, reaches the right Lagrangian term also when
  `intSparseLagPairs` is on (the default): each `LagBFunction` then holds only
  the dual pairs of the constraints its sub-Block appears in, numbered among
  themselves, while the term was looked up by the index of the constraint
  among all the relaxed ones, which rewrote the wrong term or crashed; the
  term is now found through its multiplier, and a sub-Block that enters a
  relaxed constraint for the first time gets the dual pair it lacked

- `has_var_solution()` answers false after a `compute()` that returned
  `kError` or `kBlockLocked`: it asked the inner Solver whether it had a
  dual solution, which a bundle has (the multipliers of its master problem)
  also when the evaluation of a component failed, and the caller that
  believed it had `get_var_solution()` throw "no coefficients stored" on the
  component whose linearizations were never computed

- `PrimalProximalHeur` records a recovered point only if it is one: a point
  is discarded when its value beats the bound the Lagrangian Dual gives, no
  feasible point being able to do that, and when it violates the rows of the
  Block by more than a relative 1e-6. The branch that had none of these
  checks is the one without static binary Variable to put the proximal term
  on, where the heuristic skips its loop altogether and records whatever the
  consensus recovery leaves in the Block; on a unit commitment instance
  translated from PyPSA it reported an upper bound below its own lower bound
  and handed over a point violating the dualised rows by 5781, which the
  scaling of the rows of the master problem made visible by moving the
  trajectory. Note that `Block::is_feasible()` is not the question here,
  passing as it does over the relaxed Constraint, which are exactly the ones
  a recovered point can violate

- `int_par_idx2str()` names `intRecursive` too, the array of the names having
  stopped one short of the parameters and answered with the name of another
  one

- `PrimalProximalHeur` counted, as the cost of a point, only the objectives
  of the sub-Block of its Block and not those of the Block nested into them,
  e.g., the HydroUnitBlock of a HydroSystemUnitBlock, so that the value it
  reported was not that of the solution it gave; the objectives of the
  whole subtree of each sub-Block are now copied and evaluated

- the initial multipliers of a relaxed Constraint that is reversed (a >=
  one in a minimization, with `int_LDSlv_NNMult`) were read from its dual
  without changing its sign, while `get_dual_solution()` writes them with
  the sign changed, so that a warm start, and that of `PrimalProximalHeur`
  from the duals of the relaxation, gave these multipliers the wrong sign;
  they are now read back as the multipliers that were written

- `PrimalProximalHeur` looks for the BlockSolverConfig of its recovery after
  the filename prefix of all Configuration, where it then opens it

- with `intRecursive` the components are given back to their own fathers
  when the Solver is detached, rather than to the root, which is not the
  father of a component taken below it and made the detach read past the
  end of its sub-Block

- makefile-c and makefile-s bring in MILPSolver, which PrimalProximalHeur
  needs, rather than leaving $(MILPSINC) to the including makefile

## [0.4.0] - 2026-09-13

### Changed

- `intLogVerb` of `PrimalProximalHeur` is one composite value, v + 4 * w,
  carrying the verbosity of the heuristic and that of its inner Solver

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

### Changed

- significant improvements in handling Configurations

- Lagrangian variables now initialized with dual variables from the FRowConstraint

### Fixed

- fixed a number of issues (get_lb/ub exchanged, computation of solutions)

## [0.1.0] - 2021-05-02

Initial release

### Added

- Initial release.

[Unreleased]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.4.0...develop
[0.4.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.3...0.2.0
[0.1.3]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.2...0.1.3
[0.1.2]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/tags/0.1.0
