# LagrangianDualSolver

Definition and implementation of the LagrangianDualSolver class, which
implements the CDASolver interface within the SMS++ framework for a "generic"
Lagrangian-based Solver.

This can "solve" (see below for the reason of the scare quotes) any Block
(B) with the following structure:

- no Variable in (B)

- (B) does not depend on any "external" Variable, i.e., a Variable that does
  not belong to (B) (or any of its sub-Block, recurively)

- (B) has at least one sub-Block (necessarily, for otherwise it would be
  "completely empty")

- if there is more than one sub-Block, the Constraint in (B) are all and only
  the ones that link its sub-Block; that is, no sub-Block must depend on any
  "external" Variable, i.e., a Variable that does not belong to the sub-Block
  (or any of its sub-Block, recurively)

- each sub-Block may never make any assumption on which type (B) is or make
  any reference to any of its data

The reason for the last requirement is that LagrangianDualSolver "cheats" on
(B): it stealthily constructs a new Block corresponding to its Lagrangian Dual,
physically moving the sub-Block of (B) inside it while not changing the
pointers in (B). That is, the sub-Block of (B) (temporarily) change father
Block to a new Block that remains hidden inside the LagrangianDualSolver (this
is undone when the LagrangianDualSolver is unregistered from (B)), while (B)
still "believes" that they remain its sub-Block.

The new Lagrangian Dual Block keeps consistency, in particular by
forwarding to (B) any Modification coming from the sub-Block.

An appropriate Solver is then registered to the Lagrangian Dual Block, and it
is used to solve it. The solution it used as the dual solution for (B), while
a primal solution is constructed by convexification. Here comes the reason for
the scare quotes: if (B) does not represent a convex program (say, the
sub-Block have integer variables), then the Lagrangian Dual Block is not
equivalent to (B) but to its "convexified relaxation", and this is what is
solved.

A different issue is that (B) may represent a convex program which is
"nonlinear enough" so that strong duality does not hold; say, the primal
problem may not have finite optimum (and not be unbounded), or the dual
problem may be infeasible even if the primal does have an optimal solution.
We assume that these cases are either dealt with by the user of
LagrangianDualSolver.

## Getting started

These instructions will let you build LagrangianDualSolver.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

- A SMS++ Solver capable of solving the Lagrangian Dual, such as
  [BundleSolver](https://gitlab.com/smspp/bundlesolver)

### Build and install

Configure and build the library with:
```sh
mkdir build
cd build
cmake ..
make
```

Optionally, install the library in the system with:
```sh
sudo make install
```

## Usage

After the library is configured and built, you can use it in your CMake project with:
```cmake
find_package(LagrangianDualSolver)
target_link_libraries(<my_target> SMS++:: LagrangianDualSolver)
```


## Contributing

This section is not ready yet.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
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
