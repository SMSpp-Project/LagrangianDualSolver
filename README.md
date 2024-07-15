# LagrangianDualSolver

Definition and implementation of the `LagrangianDualSolver` class, which
implements the `CDASolver` interface within the SMS++ framework for a "generic"
Lagrangian-based `Solver`.

This can "solve" (see below for the reason of the scare quotes) any Block
(B) with the following structure:

- no `Variable` in (B)

- (B) does not depend on any "external" `Variable`, i.e., a `Variable` that does
  not belong to (B) (or any of its sub-`Block`, recursively)

- (B) has at least one sub-`Block` (necessarily, for otherwise it would be
  "completely empty")

- if there is more than one sub-`Block`, the `Constraint` in (B) are all and only
  the ones that link its sub-`Block`; that is, no sub-`Block` must depend on any
  "external" `Variable`, i.e., a `Variable` that does not belong to the
  sub-`Block` (or any of its sub-sub-`Block`, recursively)

- all the `Constraint` in (B) are "linear constraint", i.e., `FRowConstraint`
  with a `LinearFunction` inside. Note that `OneVarConstraint` are "linear
  constraint" as well, but since they only concern one variable they cannot
  be "linking constraints". Although it may in principle be that one may want
  to deal with them in a Lagrangian fashion, in most of the cases including
  them in the subproblem is better, and therefore `LagrangianDualSolver`
  currently do not support them (although this may change later if a serious
  use case arises)

- each sub-`Block` may never make any assumption on which type (B) is or make
  any direct reference to any of its data

The reason for the last requirement is that `LagrangianDualSolver` "cheats" on
(B): it stealthily constructs a new `Block` corresponding to its Lagrangian Dual,
"physically moving" the sub-`Block` of (B) inside it while not changing the
pointers in (B). That is, the sub-`Block` of (B) (temporarily) change father
`Block` to a new `Block` that remains hidden inside the `LagrangianDualSolver`
(this is undone when the `LagrangianDualSolver` is unregistered from (B)), while
(B) still "believes" that they remain its sub-`Block`. Yet, consistency is kept
in that any `Modification` coming from the sub-`Block` is forwarded to (B) via
the "father of `LagBFunction`" mechanism,

An appropriate `Solver` is then registered to the Lagrangian Dual `Block`, and it
is used to solve it. The solution it used as the dual solution for (B), while
a primal solution is constructed by convexification. Here comes the reason for
the scare quotes: if (B) does not represent a convex program (say, the
sub-`Block` have integer variables), then the Lagrangian Dual `Block` is not
equivalent to (B) but to its "convexified relaxation", and this is what is
solved.

A different issue is that (B) may represent a convex program which is
"nonlinear enough" so that strong duality does not hold; say, the primal
problem may not have finite optimum (and not be unbounded), or the dual
problem may be infeasible even if the primal does have an optimal solution.
We assume that these cases either do not occur or are dealt with by the
user of LagrangianDualSolver.

## Getting started

These instructions will let you build `LagrangianDualSolver`.


### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

It's not a build requirement but you will need a SMS++ `Solver`
capable of solving the Lagrangian Dual, such as
[BundleSolver](https://gitlab.com/smspp/bundlesolver).


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(LagrangianDualSolver)
target_link_libraries(<my_target> SMS++::LagrangianDualSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `LagrangianDualSolver` has to include a "main
makefile" of the module, which typically is either [makefile-c](makefile-c)
including all necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the
[LagrangianDualSolver + BoxSolver tester in the tests/ repo](https://gitlab.com/smspp/tests/-/blob/develop/LagrangianDualSolver_Box/test.cpp?ref_type=heads).
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries are not at their default location, it should only be
necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/lagrangiandualsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Gorgone**  
  Dipartimento di Matematica ed Informatica  
  Università di Cagliari

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
