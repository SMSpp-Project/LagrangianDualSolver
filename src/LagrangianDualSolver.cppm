/*--------------------------------------------------------------------------*/
/*--------------------- File LagrangianDualSolver.cppm ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.lagrangiandualsolver, which exports
 * every name that the headers of LagrangianDualSolver declare at namespace
 * scope, together with the modules it depends on; a translation unit can then
 * write
 *
 *     import smspp.lagrangiandualsolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * The headers LagrangianDualRelaxationSolverML.h are not in this module: they
 * need a library of their own, or only a source includes them.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "LagrangianDualRelaxationSolver.h"
#include "LagrangianDualSolver.h"
#include "PrimalProximalHeur.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.lagrangiandualsolver;

export import smspp.branchandxsolver;
export import smspp.milpsolver;

export namespace SMSpp_di_unipi_it {

 // LagrangianDualRelaxationSolver.h
 using SMSpp_di_unipi_it::LagrangianChange;
 using SMSpp_di_unipi_it::LagrangianDualRelaxationSolver;

 // LagrangianDualSolver.h
 using SMSpp_di_unipi_it::LagrangianDualSolver;

 // PrimalProximalHeur.h
 using SMSpp_di_unipi_it::PrimalProximalHeur;
}

/*--------------------------------------------------------------------------*/
/*------------------- End File LagrangianDualSolver.cppm -------------------*/
/*--------------------------------------------------------------------------*/
