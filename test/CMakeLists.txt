# --------------------------------------------------------------------------- #
#    CMake file for LagrangianDualSolver/test                                 #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- LagrangianDualSolver_unit_test -------------------------------------- #
add_executable(${modName}_unit_test test.cpp)
target_link_libraries(${modName}_unit_test PRIVATE SMS++::LagrangianDualSolver)

smspp_copy_runtime(${modName}_unit_test)

# The tester builds its Block in memory and reads its Configuration files,
# which are next to it, from the working directory.
add_test(NAME ${modName}_unit_test
         COMMAND ${modName}_unit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- partial_solution_test ----------------------------------------------- #
# The parts of the solution that the "meta" Configuration of
# get_var_solution() and get_dual_solution() asks for, on a Block of three
# levels whose Lagrangian Dual is solved for real by a BundleSolver and a
# GRBMILPSolver; it reads BSPar-sol.txt and the files it names from the
# working directory
if (TARGET SMS++::BundleSolver)
    add_executable(partial_solution_test test_partial_solution.cpp)
    target_link_libraries(partial_solution_test PRIVATE
                          SMS++::LagrangianDualSolver
                          SMS++::BundleSolver)
    add_test(NAME partial_solution_test
             COMMAND partial_solution_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()

# --------------------------------------------------------------------------- #
