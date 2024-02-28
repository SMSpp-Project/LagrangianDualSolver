# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

### Changed 

### Fixed 

## [0.1.3] - 2024-02-28

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

[Unreleased]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.2...develop
[0.1.2]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/lagrangiandualsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/lagrangiandualsolver/-/tags/0.1.0
