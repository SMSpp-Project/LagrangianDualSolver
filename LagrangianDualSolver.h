/*--------------------------------------------------------------------------*/
/*---------------------- File LagrangianDualSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the LagrangianDualSolver class, which implements the
 * CDASolver interface within the SMS++ framework for a "generic"
 * Lagrangian-based Solver.
 *
 * \version 0.01
 *
 * \date 11 - 11 - 2020
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * Copyright &copy by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LagrangianDualSolver
 #define __LagrangianDualSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class AbstractBlock;  // forward declaration of AbstractBlock

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup LagrangianDualSolver_CLASSES Classes in LagrangianDualSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BundleSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver soving the Lagrangian Dual of a "generic" Block
/** The LagrangianDualSolver class implements the CDASolver interface within
 * the SMS++ framework for a "generic" Lagrangian-based Solver.
 *
 * This can "solve" (see below for the reason of the scare quotes) any Block
 * (B) with the following structure:
 *
 * - No Variable in (B).
 *
 * - (B) has at least one sub-Block (necessarily, for otherwise it would be
 *   "completely empty").
 *
 * - No Objective in (B) (which makes sense: since all Variable in (B)
 *   actually belong to the sub-Block, recurively, it's them who define
 *   the Objective for these Variable: no need for (B) to do it).
 *
 * - (B) and all its sub-Block, recursively, do not depend on any "external"
 *   Variable, i.e., a Variable that does not belong to (B) (actually, to any
 *   of its sub-Block, recurively, since (B) cannot have any Variable of its
 *   own).
 *
 * - If there is more than one sub-Block, the Constraint in (B) are all and
 *   only the ones that link the sub-Block between them; that is, no sub-Block
 *   must depend on any "external" Variable, i.e., a Variable that does not
 *   belong to the sub-Block (or any of its sub-sub-Block, recurively).
 *
 * - All the Constraint in (B) are "linear constraint", i.e., FRowConstraint
 *   with a LinearFunction inside. Note that OneVarConstraint are "linear
 *   constraint" as well, but since they only concern one variable they
 *   cannot be "linking constraints". Although it may in principle be that
 *   one may want to deal with them in a Lagrangian fashion, in most of the
 *   cases including them in the subproblem is better, and therefore
 *   LagrangianDualSolver currently do not support them (although this may
 *   change later if a serious use case arises).
 *
 * - Each sub-Block of (B) may never make any assumption on which type (B) is
 *   or make any direct reference to any of its data.
 *
 * The reason for the last requirement is that LagrangianDualSolver may
 * "cheat" on (B): it stealthily constructs a new Block corresponding to its
 * Lagrangian Dual, possibly "physically moving" the sub-Block of (B) inside
 * it while not changing the pointers in (B). That is, the sub-Block of (B)
 * temporarily change father Block to a new Block that remains hidden inside
 * the LagrangianDualSolver, while (B) still "believes" that they remain its
 * sub-Block. This is undone when the LagrangianDualSolver is unregistered
 * from (B). Consistency is kept, in that any Modification coming from the
 * sub-Block is also "forwarded" to (B).
 *
 * An appropriate Solver is then registered to the Lagrangian Dual Block, and
 * it is used to solve it. The solution it used as the dual solution for (B),
 * while a primal solution is constructed by convexification. Here comes the
 * reason for the scare quotes: if (B) does not represent a convex program
 * (say, the sub-Block have integer variables), then the Lagrangian Dual
 * Block is not equivalent to (B) but to its "convexified relaxation", and
 * this is what is solved.
 *
 * IMPORTANT NOTE ON TWO-SIDED CONSTRAINTS. The FRowConstraint in (B) in
 * general have the form l <= ax <= u, i.e., they correspond to *two* linear
 * constraints. However, in many cases only *one* Lagrangian multiplier need
 * be defined for them:
 *
 * - if l == u, i.e., the equality constraint ax = u (= l); in this case the
 *   corresponding Lagrangian multiplier is constrained in  sign;
 *
 * - if l == -INF and u < INF, i.e., the less-than constraint ax <= u; in
 *   this case the corresponding Lagrangian multiplier is constrained in
 *   sign (>= 0 if the (B) is max, <= 0 if (B) is min);
 *
 * - if l > -INF and u == INF, i.e., the greater-than constraint ax >= l; in
 *   this case the corresponding Lagrangian multiplier is constrained in
 *   sign (<= 0 if the (B) is max, >= 0 if (B) is min);
 *
 * Save for the degenerate case l == -INF and u == INF, which is not allowed,
 * this leves the case -INF < l < u < INF. One possible approach for this
 * would be to consider the constraint as actually being the two less-than
 * and greater-then (ax <= u, a >= l) with two different Lagrangian
 * multipliers, both constrained in sign (in the right way). However, this
 * would significantly complicate the handling of these constraints since
 * each original one may give rise to either one or two multipliers, which
 * would be very though especially if changing the rhs/lhs of the constraint
 * would change its two-sidedness status (say, an equality constraint becoming
 * a ranged one, or an INF bound becoming finite). A different approach is to
 * reformulate the constraint as
 *
 *   ax - s = 0  ,  l <= s <= u
 *
 * and relax it, with an *unconstrained* multiplier (call it y). This would
 * lead to the same single Lagrangian term y ( ax ), plus the extra "mini
 * Lagrangian subproblem"
 *
 *   min/max { y ( - s ) : l <= s <= u }
 *
 * All these may be gathered into a single very simple LagBFunction, to
 * which possibly a VerySimpleLPSolver could be attached. Alternatively,
 * each of these may be represented as a separate one-variable LagBFunction
 * (again, possibly solved by a VerySimpleLPSolver). In particular for Solver
 * capable of exploiting the structure of the LagBFunction to properly modify
 * the Master Problem, such a reformulation would basically be a 0-cost one
 * and most likely the best approach. If one really wants to handle all the
 * cases of changes in the lhs/rhs, even those hanging the two-sidedness
 * status, this (these) extra LagBFunction(s) still should be dynamic and
 * allow new s variables to be created and destroyed; yet, this is in general
 * possible.
 *
 * Hence, LagrangianDualSolver ASSUMES ONLY ONE MULTIPLIER PER RELAXED
 * CONSTRAINTS IN ALL CASES. However THE CONSTRUCTION OF THE
 * "MINI-LagBFunction" FOR THE s VARIABLE IS NOT SUPPORTED YET, WHICH MEANS
 * THAT TRUE TWO-SIDED FRowConstraint ARE NOT ALLOWED YET. Fortunately, true
 * two-sided FRowConstraint are rare in practice, and they can always be
 * avoided by explicitly modelling them as the less-than and greater-than
 * version if needed. Yet, the mini-LagBFunction will hopefully one day be
 * actually handled.
 *
 * A different issue is that (B) may represent a convex program which is
 * "nonlinear enough" so that strong duality does not hold; say, the primal
 * problem may not have finite optimum (and not be unbounded), or the dual
 * problem may be infeasible even if the primal does have an optimal solution.
 * We assume that these cases either do not occur or are dealt with by the
 * user of LagrangianDualSolver. */

 class LagrangianDualSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from Function and C05Function.
 *
 *  @{ */

 using Index = Block::Index;
 using c_Index = Block::c_Index;

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;

/*----------------------------- CONSTANTS ----------------------------------*/

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that LagrangianDualSolver has in addition to these of CDASolver. The 
  * value intLastLDSSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum int_par_type_LDSlv {

 int_LDSlv_iBCopy = intLastParCDAS ,
 ///< if the R3Block has be used for the father block

 intLastLDSSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_LDSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that LagrangianDualSolver has in addition to these of CDASolver. The
  * value dblLastLDSSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum dbl_par_type_LDSlv {
  dblLastLDSlvPar = dblLastParCDAS ,
   ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set
   * of double algorithmic parameters. */

  };  // end( dbl_par_type_LDSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters
 /** Public enum describing the different algorithmic parameters of string
  * type that LagrangianDualSolver has in addition to these of CDASolver. The
  * value strLastLDSSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum str_par_type_LDSlv {
  str_LDSlv_ISName = strLastParCDAS ,  ///< classname of the inner Solver
  
  strLastLDSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of string parameters. */

  };  // end( str_par_type_LDSlv )

/*@} -----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING LagrangianDualSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing LagrangianDualSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 LagrangianDualSolver( void ) : CDASolver() , NumVar( 0 ) , f_nsb( 0 ) ,
  f_convex( false ) , LagrDual( nullptr ) , f_LDBConfig( nullptr ) ,
  static_cons( 0 )
 {
  // ensure all parameters are properly given their default value
  iBCopy  = dflt_int_par[ int_LDSlv_iBCopy - intLastParCDAS ];
  ISName  = dflt_str_par[ str_LDSlv_ISName - intLastParCDAS ];

  // ensure that the inner Solver is always well defined
  InnerSolver = new_Solver( ISName );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the LagrangianDualSolver from the Block

 virtual ~LagrangianDualSolver() { set_Block( nullptr ); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the int paramaters of LagrangianDualSolver / the inner Solver
 /** Set the int paramaters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double paramaters of LagrangianDualSolver / the inner Solver
 /** Set the double paramaters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , double value ) override {
  InnerSolver->set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the string paramaters of LagrangianDualSolver / the inner Solver
 /** Set the string paramaters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-int paramaters of LagrangianDualSolver / inner Solver
 /** Set the vector-of-int paramaters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , std::vector< int > && value ) override {
  InnerSolver->set_par( par , std::move( value ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// set the vector-of-double paramaters of LagrangianDualSolver / inner Solver
 /** Set the vector-of-double paramaters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , std::vector< double > && value ) override {
  InnerSolver->set_par( par , std::move( value ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// set the vector-of-string paramaters of LagrangianDualSolver / inner Solver
 /** Set the vector-of-string paramaters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , std::vector< std::string > && value ) override {
  InnerSolver->set_par( par , std::move( value ) );
  }

/*--------------------------------------------------------------------------*/
 /// "translate" an int parameter index of the inner Solver
 /** Takes the index \p par of an int parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type int_par_is( idx_type par ) {
  if( par >= intLastParCDAS )
   par += intLastLDSSlvPar - intLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a double parameter index of the inner Solver
 /** Takes the index \p par of a double parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type dbl_par_is( idx_type par ) {
  if( par >= dblLastParCDAS )
   par += dblLastLDSSlvPar - dvlLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a string parameter index of the inner Solver
 /** Takes the index \p par of a string parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type str_par_is( idx_type par ) {
  if( par >= strLastParCDAS )
   par += strLastLDSSlvPar - strLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-int parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-int parameter of the inner Solver
  * and returns the index that has to be passed to LagrangianDualSolver to
  * have that very same parameter set in the inner Solver; see the comments
  * to set_ComputeConfig() for details. */

 idx_type vint_par_is( idx_type par ) { return( par ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-double parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-double parameter of the inner
  * Solver and returns the index that has to be passed to
  * LagrangianDualSolver to have that very same parameter set in the inner
  * Solver; see the comments to set_ComputeConfig() for details. */

 idx_type vdbl_par_is( idx_type par ) { return( par ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-string parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-string parameter of the inner
  * Solver and returns the index that has to be passed to
  * LagrangianDualSolver to have that very same parameter set in the inner
  * Solver; see the comments to set_ComputeConfig() for details. */

 idx_type vstr_par_is( idx_type par ) { return( par ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a int parameter index of the LagrangianDualSolver
 /** Takes the index \p par of an int parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type int_par_lds( idx_type par ) {
  if( par >= intLastLDSSlvPar )
   par -= intLastLDSSlvPar - intLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// "translate" a double parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a double parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type dbl_par_lds( idx_type par ) {
  if( par >= dblLastLDSSlvPar )
   par -= dblLastLDSSlvPar - sblLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// "translate" a string parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a string parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type str_par_lds( idx_type par ) {
  if( par >= strLastLDSSlvPar )
   par -= strLastLDSSlvPar - strLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-int parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-int parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vint_par_lds( idx_type par ) { return( par ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// translate a vector-of-double parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-double parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vdbl_par_lds( idx_type par ) { return( par ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// translate a vector-of-string parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-string parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vstr_par_lds( idx_type par ) { return( par ); }

/*--------------------------------------------------------------------------*/
 /// set the whole set of parameters in one blow
 /** This method sets the whole set of parameters in one blow using a
  * ComputeConfig object.
  *
  * LagrangianDualSolver is, in some sense, no more than a wrapper of the
  * "true" CDASolver used to solve the Lagrangian Dual (except that this
  * wrapping is nontrivial). As such, it has comparatively few algorithmic
  * parameters, while the inner Solver may have many. It is therefore
  * advantageous to allow to set the algorithmic parameters of the inner
  * 
  *


  *
  * The method of LagrangianDualSolver calls
  * ThinComputeInterface::set_ComputeConfig() to do the bulk of the work, and
  * then manages f_extra_Configuration. The field must be non-nullptr at
  * least once before compute() and any follow up-method is called, and it
  * can contain:
  *
  * - a BlockSolverConfig *
  *
  * - a SimpleConfiguration< std::pair< Configuration * , Configuration * > >
  *   there .first is a BlockSolverConfig * and .second is a
  *   BlockConfig *
  *
  * The BlockSolverConfig is apply()-ed to the Lagrangian Dual Block; it
  * has to register *at least* one (appropriate) CDASolver to it, and it
  * also has the chance to register Solver to the inner Block (inside the
  * LagBFunction inside the FRealObjective of each sub-Block). The *first*
  * CDASolver registered to the Lagrangian Dual Block will be used as the
  * "main" Solver by LagrangianDualSolver, and there typically must be at
  * least one Solver registered to each sub-Block (unless the "main" Solver
  * allows not to because it deals with some LagBFunction in specialised
  * ways), the first of which will be used by the LagBFunction to compute()
  * themselves.
  *
  * If a BlockConfig is provided, it is also apply()-ed to the Lagrangian
  * Dual Block (and therefore it has the chance to BlockConfig-ure also
  * all the inner Block inside the LagBFunction inside the FRealObjective of
  * each sub-Block).
  *
  * The BlockSolverConfig becomes "property" of the LagrangianDualSolver;
  * it is clear()-ed and used when the LagrangianDualSolver is unregistered
  * from the Block to clear away all the Solver that it has attached.
  *
  * Note that, if the LagrangianDualSolver is not registered to a Block
  * when set_ComputeConfig() is called, the pointers are stored and the
  * configuration of the Block is done in set_Block(). For the
  * BlockSolverConfig, this could be repeated any number of times since a
  * BlockSolverConfig is not "consumed" when it is apply()-ed. However, a
  * BlockConfig is "consumed" instead. Thus, for consistency also the
  * BlockSolverConfig is immediately clear()-ed: registering the
  * LagrangianDualSolver to a new Block requires a new call to
  * set_ComputeConfig().
  *
  * Note that if multiple calls to set_ComputeConfig() occur before the
  * LagrangianDualSolver is registered to a Block, all the BlockSolverConfig
  * / BlockConfig of all calls save the last one are lost and have no effect
  * (but at least are properly deleted). */

 void set_ComputeConfig( ComputeConfig * scfg = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the LagrangianDualSolver log

 void set_log( std::ostream *log_stream = nullptr ) override {
  InnerSolver->set_log( f_log );
  }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Set event handlers
 *
 * Since LagrangianDualSolver basically only acts as a "front end" for the
 * "inner Solver" that actually solves the Lagrangian Dual, it does not
 * handle the events itself; rather, it passes them through to the "true"
 * Solver.
 *
 *  @{ */

 EventID set_event_handler( int type , EventHandler && event ) override {
  return( InnerSolver->set_event_handler( type , std::move( event ) ) );
  }

/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/

 void reset_event_handler( int type , EventID id ) override {
  InnerSolver->reset_event_handler( type , id );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Lagrangian Dual of the given Block
 *  @{ */

 /// (try to) solve the Lagrangian Dual of the given Block

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
 /// returns the "inner" CDASolver used to solve the Lagrangian Dual
 /** Returns a pointer to the "inner" CDASolver used to solve the Lagrangian
  * Dual. This should not be necessary since LagrangianDualSolver makes it
  * possible to do most of the useful operations (in particular Configuring
  * it) via the standard interface, but this is still provided for complete
  * generality. */

 CDASolver * get_inner_Solver( void ) { return( InnerSolver ); }

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 OFValue get_lb( void ) override { return( InnerSolver->get_ub() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { return( InnerSolver->get_lb() ); }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override {
  return( InnerSolver->has_dual_solution() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_var_solution( void ) override {
  return( InnerSolver->has_dual_solution() );
  }

/*--------------------------------------------------------------------------*/

 bool is_var_feasible( void ) override {
  return( InnerSolver->is_dual_feasible() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool is_dual_feasible( void ) override {
  return( InnerSolver->is_var_feasible() );
  }

/*--------------------------------------------------------------------------*/
 /// write the "current" solution

 void get_var_solution( Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution

 void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override {
  return( InnerSolver->new_dual_solution() )
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_dual_solution( void ) override {
  return( InnerSolver->new_var_solution() );
  }

/*--------------------------------------------------------------------------*/
/*
  void set_unbounded_threshold( const VarValue thr ) override { }

  bool has_var_direction( void ) override { return( true ); }

  bool has_dual_direction( void ) override { return( true ); }

  void get_var_direction( Configuration *dirc = nullptr ) override {}

  void get_dual_direction( Configuration *dirc = nullptr ) override {}

  virtual bool new_var_direction( void ) override { return( false ); }
  
  virtual bool new_dual_direction( void ) override{ return( false ); }
*/

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the LagrangianDualSolver
 *
 * While LagrangianDualSolver itself has comparatively few parameters, it
 * allows to change via its standard parameter interface all the parameters
 * of the "inner Solver" used to actually solve the Lagrangian Dual. The
 * indices of these parameters are automatically translated (see *_par_is()
 * and the comments to set_ComputeConfig() for details) so that they can be
 * automatically set and queried as if they were "natural" parameters of
 * LagrangianDualSolver itself.
 * @{ */

 idx_type get_num_int_par( void ) const override {
  return( int_par_is( InnerSolver->get_num_int_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( dbl_par_is( InnerSolver->get_num_dbl_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_str_par( void ) const override {
  return( str_par_is( InnerSolver->get_num_str_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_vint_par( void ) const override {
  return( vint_par_is( InnerSolver->get_num_vint_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_vdbl_par( void ) const override {
  return( vdbl_par_is( InnerSolver->get_num_vdbl_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_vstr_par( void ) const override {
  return( vstr_par_is( InnerSolver->get_num_vstr_par() ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_dflt_int_par( idx_type par ) const override {
  if( ( par >= intLastParCDAS ) && ( par < intLastLdsSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );

  return( InnerSolver->get_dflt_int_par( int_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( idx_type par ) const override {
  if( ( par >= dblLastParCDAS ) && ( par < dblLastLdsSlvPar ) )
   return( dflt_dbl_par[ par - dblLastParCDAS ] );

  return( InnerSolver->get_dflt_dbl_par( dbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 const std::string & get_dflt_str_par( idx_type par ) const override {
  if( ( par >= strLastParCDAS ) && ( par < strLastLdsSlvPar ) )
   return( dflt_str_par[ par - strLastParCDAS ] );

  return( InnerSolver->get_dflt_str_par( str_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< int > & get_dflt_vint_par( idx_type par )
  const override {
  return( InnerSolver->get_dflt_vint_par( vint_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_dflt_vdbl_par( idx_type par )
  const override {
  return( InnerSolver->get_dflt_vdbl_par( vdbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< std::string > & get_dflt_vstr_par( idx_type par )
  const override {
  return( InnerSolver->get_dflt_vstr_par( vstr_par_lds( par ) ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & get_str_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< int > & get_vint_par( idx_type par ) const override {
  return( InnerSolver->get_vint_par( vint_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_vdbl_par( idx_type par ) const override {
  return( InnerSolver->get_vdbl_par( vdbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< std::string > & get_vstr_par( idx_type par )
  const override {
  return( InnerSolver->get_vstr_par( vstr_par_lds( par ) ) );
  }

/*--------------------------------------------------------------------------*/

 idx_type int_par_str2idx( const std::string & name ) const override {
  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );

  return( int_par_is( InnerSolver->int_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );

  return( dbl_par_is( InnerSolver->dbl_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type str_par_str2idx( const std::string & name ) const override {
  const auto it = str_pars_map.find( name );
  if( it != str_pars_map.end() )
   return( it->second );

  return( str_par_is( InnerSolver->str_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type vint_par_str2idx( const std::string & name ) const override {
  return( vint_par_is( InnerSolver->vint_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type vdbl_par_str2idx( const std::string & name ) const override {
  return( vdbl_par_is( InnerSolver->vdbl_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type vstr_par_str2idx( const std::string & name ) const override {
  return( vstr_par_is( InnerSolver->vstr_par_str2idx( name ) ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( idx_type idx ) const override {
  if( ( idx >= intLastParCDAS ) && ( idx < intLastLdsSlvPar ) )
   return( int_pars_str[ idx - intLastParCDAS ] );

  return( InnerSolver->int_par_idx2str( int_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( idx_type idx ) const override {
  if( ( idx >= dblLastParCDAS ) && ( idx < dblLastLdsSlvPar ) )
   return( dbl_pars_str[ idx - dblLastParCDAS ] );

  return( InnerSolver->dbl_par_idx2str( dbl_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & str_par_idx2str( idx_type idx ) const override {
  if( ( idx >= strLastParCDAS ) && ( idx < strLastLdsSlvPar ) )
   return( str_pars_str[ idx - strLastParCDAS ] );

  return( InnerSolver->str_par_idx2str( str_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & vint_par_idx2str( idx_type idx ) const override {
  return( InnerSolver->vint_par_idx2str( vint_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & vdbl_par_idx2str( idx_type idx ) const override {
  return( InnerSolver->vdbl_par_idx2str( vdbl_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & vstr_par_idx2str( idx_type idx ) const override {
  return( InnerSolver->vstr_par_idx2str( vstr_par_lds( idx ) ) );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void register_inner_Solver( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void unregister_inner_Solver( void );

/*--------------------------------------------------------------------------*/

 void set_default_inner_BlockSolverConfig( void );

/*--------------------------------------------------------------------------*/

 void configure_LagrangianDualBlock( void );

/*--------------------------------------------------------------------------*/
 /** Returns the index as active variable of the LagBFunction of the given
  * FRowConstraint, be it static or dynamic.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variablee, Int< Index >() if
  *         \p con does not correspond to any FRowConstraint */

 Index index_of_constraint( const FRowConstraint * con ) {
  auto i = index_of_static_constraint( con );
  if( i < Inf< Index >() )
   return( i );

  return( index_of_dynamic_constraint( con ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the index as active variable of the LagBFunction of the given
  * static FRowConstraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variable, Int< Index >() if
  *         \p con does not correspond to any static FRowConstraint */

 Index index_of_static_constraint( const FRowConstraint * con );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the index as active variable of the LagBFunction of the given
  * dynamic FRowConstraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variable, Int< Index >() if
  *         \p con does not correspond to any dynamic FRowConstraint */

 Index index_of_dynamic_constraint( const FRowConstraint * con );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to) FRowConstraint corresponding to a the active
  * variable of the LagBFunction with the given index, be it static or dynamic.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint */
 
FRowConstraint * constraint_with_index( Index i ) {
 #ifdef NDEBUG
  if( i >= NumVar )
   throw( std::invalid_argument(
		     "LagrangianDualSolver::invalid index of constraint" ) );
 #endif

 if( i < static_cons )
  return( static_constraint_with_index( i ) );

 return( dynamic_constraint_with_index( i ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to the) static FRowConstraint corresponding to a the
  * active variable of the LagBFunction with the given index.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint */

 FRowConstraint * static_constraint_with_index( Index i );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to the) static FRowConstraint corresponding to a the
  * active variable of the LagBFunction with the given index.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint */

 FRowConstraint * dynamic_constraint_with_index( Index i ) {
  return( idx_to_dcon[ i - static_cons ] );
  }

/*--------------------------------------------------------------------------*/

 Index Block2Index( Block * blck ) {
  for( ;; ) {
   auto fb = blck->get_f_Block();
   if( fb == LagrDual ) {
    auto it = std::lower_bound( blck_to_idx.begin() , blck_to_idx.end() ,
				blck_int_pf( blck , 0 , nullptr ) ,
				[]( const auto & a , const auto & b ) {
				 return( std::get< 0 >( a ) <
					 std::get< 0 >( b ) ); } );
    return( std::get< 1 >( *it ) );
    }
   blck = fb;
   if( ! blck )
    throw( std::invalid_argument(
	        "LagrangianDualSolver: Variable belonging to wrong Block" ) );
   }
  }

/*--------------------------------------------------------------------------*/

 void split_constraint( const FRowConstraint & con ,
			std::vector< LinearFunction::v_coeff_pair > & split );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorthmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool iBCopy;         ///< if the R3Block conversion has to be done

 std::string ISName;  ///< classname of the inner Solver

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumVar;      ///< (current) number of variables

 Index f_nsb;       ///< number of sub-Block

 bool f_convex;     ///< true if (B) was a max problem, false otherwise

 AbstractBlock * LagrDual;  ///< the automatically constructed Lagrangian Dual

 CDASolver * InnerSolver;   ///< the Solver attached to LagrDual

 Configuration * f_LDBConfig;   ///< the Configuration for LagrDual

 // dictionaries- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




 Index static_cons;  ///< number of static constraints
                     /** Total number of static constraints in the Block: the
		      * first static_cons active variables in the LagBFunction
		      * will never change. */

 /* The following vectors are used in order to keep track between the
  * Constraints of the Block and the Lagrangian variables of all the
  * corresponding Lagrangian functions.
  *
  *  - scon_to_idx: vector of tuples that store 1) the address of the first
  *    element of each group of static constraints, respectively, 2) the
  *    corresponding index in the set of active variables of the Lagrangian
  *    functions, and 3) the number of elements in the group. the vector is
  *    kept sorted in ascending order by address
  *
  *  - dcon_to_idx: vector of pairs that store the address of each dynamic
  *    constraints the corresponding index in the set of active variables of
  *    the Lagrangian function; the vector is kept sorted in ascending order
  *    by address
  *
  *  - idx_to_scon: vector of pairs that store the index of the first active
  *    variable of the Lagrangian functions corresponding to each group of
  *    static constraints and the address of the first constraint of the
  *    group; the vector is kept sorted in ascending order by index
  *
  *  - idx_to_dcon: vector of pairs that store index of each active variable
  *    of the Lagrangian functions corresponding to a dynamic constraints and 
  *    the address of that constraint; the vector is kept sorted in ascending
  *    order by index
  *
  * Using these dictionaries we can efficiently find the index of each
  * constraint as a Lagrangian variable and vice-versa. */

 typedef std::pair< FRowConstraint * , Index > const_int;
 typedef std::pair< Index , FRowConstraint * > int_const;
 typedef std::tuple< FRowConstraint * , Index , Index > con_int_int;

 std::vector< con_int_int > scon_to_idx; ///< from static constraint to index
 std::vector< int_const > idx_to_scon;   ///< from index to static constraint

 std::vector< const_int > dcon_to_idx;   ///< from dynamic constraint to index
 std::vector< FRowConstraint * > idx_to_dcon;
 ///< From index to dynamic constraint

 typedef std::tuple< AbstractBlock * , Index , LagBFunction * > blck_int_pf;

 /** This is a vector of triples that contain the pointer of a sub-Block of
  * the LagrangianDual, its index into the vector od sub-Block, and the
  * pointer to the corresponding LagBFunction in the FRealObjective; this is
  * ordered by Block. */

 std::vector< blck_int_pf > blck_to_idx;  ///< from Block * to index


 
/*--------------------------------------------------------------------------*/

 const static std::vector< int > dflt_int_par;
 ///< the (static const) vector of int parameters default values

 const static std::vector< double > dflt_dbl_par;
 ///< the (static const) vector of double parameters default values

 const static std::vector< double > dflt_str_par;
 ///< the (static const) vector of string parameters default values

 const static std::vector< std::string > int_pars_str;
 ///< the (static const) vector of int parameters names

 const static std::vector< std::string > dbl_pars_str;
 ///< the (static const) vector of double parameters names

 const static std::vector< std::string > str_pars_str;
 ///< the (static const) vector of string parameters names

 const static std::map< std::string , idx_type > int_pars_map;
  ///< the (static const) map for int parameters names

 const static std::map< std::string , idx_type > dbl_pars_map;
 ///< the (static const) map for double parameters names

 const static std::map< std::string , idx_type > str_pars_map;
 ///< the (static const) map for string parameters names

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class LagrangianDualSolver )

/*@}  end( group( LagrangianDualSolver_CLASSES ) ) -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* LagrangianDualSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File LagrangianDualSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
