/*--------------------------------------------------------------------------*/
/*---------------------- File LagrangianDualSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the LagrangianDualSolver class, which implements the
 * CDASolver interface within the SMS++ framework for a "generic"
 * Lagrangian-based Solver.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LagrangianDualSolver
 #define __LagrangianDualSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"

#include "CDASolver.h"

#include "LagBFunction.h"

#include "LinearFunction.h"

#include "UpdateSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class FRowConstraint;  // forward definition of FRowConstraint
  
/*--------------------------------------------------------------------------*/
/*-------------------- CLASS LagrangianDualSolver --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver solving the Lagrangian Dual of a "generic" Block
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
 *   actually belong to the sub-Block, recursively, it's them who define
 *   the Objective for these Variable: no need for (B) to do it).
 *
 * - (B) and all its sub-Block, recursively, do not depend on any "external"
 *   Variable, i.e., a Variable that does not belong to (B) (actually, to any
 *   of its sub-Block, recursively, since (B) cannot have any Variable of its
 *   own).
 *
 * - If there is more than one sub-Block, the Constraint in (B) are all and
 *   only the ones that link the sub-Block between them; that is, no sub-Block
 *   must depend on any "external" Variable, i.e., a Variable that does not
 *   belong to the sub-Block (or any of its sub-sub-Block, recursively).
 *
 * - All the Constraint in (B) are "linear constraint", i.e., FRowConstraint
 *   with a LinearFunction inside. Note that OneVarConstraint are "linear
 *   constraint" as well, but since they only concern one variable they
 *   cannot be "linking constraints". Although it may in principle be that
 *   one may want to deal with them in a Lagrangian fashion, in most of the
 *   cases including them in the subproblem is better, and therefore
 *   LagrangianDualSolver currently does not support them (although this may
 *   change later if a serious use case arises).
 *
 * - Each sub-Block of (B) may never make any assumption on which type (B) is
 *   or make any direct reference to any of its data.
 *
 * The reason for the last requirement is that LagrangianDualSolver may
 * "cheat" on (B): it stealthily constructs a new Block LB corresponding to
 * its Lagrangian Dual, possibly "physically moving" the sub-Block of (B)
 * inside it while not changing the pointers in (B). In fact, LB is
 * constructed with as many sub-Block as those of (B); for each sub-Block
 * (B_i) of (B) a sub-Block (LB_i) of (LB) is constructed. (LB_i) is empty
 * save for a FRealObjective that contains a LagBFunction; in turn, (B_i)
 * is set as the Block inside that LagBFunction. That is, the father of
 * (B_i) is set to the LagBFunction in (LB_i) (which is both a Function and
 * a Block, hence can be a father), while (B) still "believes" that (B_i)
 * remains its sub-Block. This is undone when the LagrangianDualSolver is
 * unregistered from (B). Consistency is kept, in that any Modification
 * coming from (B_i) is also "forwarded" to (B) by the trick of setting (B)
 * as the father Block of the LagBFunction, exploiting the support that
 * LagBFunction provides for this (somehow tricky) use case. Note that these
 * Modification also reach the LagBFunction, where they are "translated" into
 * FunctionMod; these are the, forwarded to (LB_i) (which Observe-s the
 * LagBFunction via the FRealObjective) and then to (LB) and all the Solver
 * registered to it (not to any Solver registered to the LagBFunction, since
 * there must not be any). Thus, the Modification reach both (LB) (having
 * been properly "translated") and the original father (B) of (B_i). This also
 * has the advantage of partly reconstructing the two-way link between (B)
 * and the (B_i): not only "going down from (B) one reaches the (B_i)", but
 * also "going up from the (B_i) one eventually gets to (B)". The mapping is
 * not perfectly kept since now (B_i) is a "grandson" of (B) rather than a
 * "son" (the LagBFunction having been slotted in the middle), but this is
 * still good enough for any operation that just requires to be able to
 * eventually reach up (B) from (B_i) up the father Block chain. An important
 * example of such an operation is MILPSolver::is_mine(), which checks if a
 * given Block is a sub-Block of the one the MILPSolver is registered to. By
 * the above trick the operation still works even if a LagrangianDualSolver
 * has been registered to (B) before the MILPSolver is (while otherwise it
 * would break).
 *
 * In case this is still not enough, i.e., something in (B) or in any other
 * Solver attached to it requires to keep the original father-son relationship
 * between (B) and its (B_i), the parameter int_LDSlv_iBCopy allows to
 * instruct LagrangianDualSolver to rather build a copy (B'_i) of (B_i) and
 * insert (B'_i) into the LagBFunction, with an UpdateSolver forwarding all
 * Modification from (B_i) to (B'_i). However
 *
 *     THIS REQUIRES get_R3_Block( nullptr ) AND map_back_solution() TO
 *     WORK FOR ALL (B_i)
 *
 * Mathematically speaking, the original (B) can be seen as
 *
 *  (B)   max / min 0
 *                  l <= \sum_{k \in K} g^k( x^k ) <= u
 *                  B^k    k \in K
 *
 * where
 *
 *   (B^k)   max / min c^k(x^k) : x^k \in X^k
 *
 * for each k \in K are all its sub-Block. There is basically no requirement
 * on X^k, save that the corresponding Solver be able to "efficiently" solve
 * the corresponding problem (albeit this can be done only approximately);
 * the same would be true for c^k(), except that we currently need an "easy"
 * function whose abstract representation we can manipulate since this is how
 * LagBFunction works. However, note that arbitrarily complex Objective could
 * be present in any sub-Block of B^k, again provided that the Solver can
 * handle them.
 *
 * It is important to remark that LagrangianDualSolver does not solve (B)
 * but rather its Lagrangian Dual. This is obtained by relaxing the linking
 * constraints, and here immediately comes some caveat. Indeed, the linking
 * FRowConstraint in (B) in general have the form l <= g( x ) <= u, i.e.,
 * they correspond to *two* linear constraints. However, in many cases only
 * *one* Lagrangian multiplier need be defined for them:
 *
 * - if l == u, i.e., the equality constraint g( x ) = u (= l), in which case
 *   the corresponding Lagrangian multiplier is unconstrained in sign;
 *
 * - if l == -INF and u < INF, i.e., the less-than constraint g( x ) <= u, in
 *   which case the corresponding Lagrangian multiplier is constrained in
 *   sign;
 *
 * - if l > -INF and u == INF, i.e., the greater-than constraint g( x ) >= l,
 *   in which case the corresponding Lagrangian multiplier is constrained in
 *   sign.
 *
 * Save for the degenerate case l == -INF and u == INF, which is not allowed,
 * this leaves the case -INF < l < u < INF. One possible approach for this
 * would be to consider the constraint as actually being the two less-than
 * and greater-than [g( x ) <= u, g( x ) >= l] with two different Lagrangian
 * multipliers, both constrained in sign (in the right way). However, this
 * would significantly complicate the handling of these constraints since
 * each original one may give rise to either one or two multipliers, which
 * would be very though especially if changing the rhs/lhs of the constraint
 * would change its two-sidedness status (say, an equality constraint becoming
 * a ranged one, or an INF bound becoming finite). A different approach is to
 * reformulate the constraint as
 *
 *   g( x ) - s = 0  ,  l <= s <= u
 *
 * and relax it, with an *unconstrained* multiplier (call it y). This would
 * lead to the same single Lagrangian term y g( x ), plus the extra "mini
 * Lagrangian subproblem"
 *
 *   max / min { y ( - s ) : l <= s <= u }
 *
 * All these may be gathered as a very simple inner Block into a single very
 * simple LagBFunction, to which possibly a BoxSolver could be attached.
 * Alternatively, each of these may be represented as a separate one-variable
 * LagBFunction (again, possibly solved by a BoxSolver). Even better, Solver
 * capable of exploiting the structure of the LagBFunction to properly modify
 * the Master Problem could see this to reformulate the Master Problem at
 * basically 0-cost, which would most likely be the best approach. If one
 * really wants to handle all the cases of changes in the lhs/rhs, even those
 * changing the two-sidedness status, this (these) extra LagBFunction(s)
 * would need to be dynamic and allow new s variables to be created and
 * destroyed; yet, this is in general possible.
 *
 * Hence, LagrangianDualSolver ASSUMES ONLY ONE MULTIPLIER PER RELAXED
 * CONSTRAINTS IN ALL CASES. However, THE CONSTRUCTION OF THE
 * "MINI-LagBFunction" FOR THE s VARIABLE IS NOT SUPPORTED YET, WHICH MEANS
 * THAT TRUE TWO-SIDED FRowConstraint ARE NOT ALLOWED YET. Fortunately, true
 * two-sided FRowConstraint are rare in practice, and they can always be
 * avoided by explicitly modeling them as the less-than and greater-than
 * version if needed. Yet, the mini-LagBFunction will hopefully one day be
 * actually handled.
 *
 * A different issue is that (B) may represent a convex program which is
 * "nonlinear enough" so that strong duality does not hold; say, the primal
 * problem may not have finite optimum (and not be unbounded), or the dual
 * problem may be infeasible even if the primal does have an optimal solution.
 * We assume that these cases either do not occur or are dealt with by the
 * user of LagrangianDualSolver.
 *
 * A final important note regards the choice of the "sign" in the relaxed
 * constraints, since this impact on the sign constraints of the Lagrangian
 * multipliers. Rewriting (B) for simplicity as
 *
 *   (B)   max / min c(x) : l <= g( x ) <= u , x \in X
 *
 * the exact form of the Lagrangian Dual, and of its Variable y, can actually
 * be done in different ways. Since y will be used as the dual value for the
 * constraints, we need to be consistent with the standards set by
 * RowConstraint [see RowConstraint.h], which are now recalled and discussed.
 * Let us start with the minimization case; then, the Lagrangian function
 * associated with (B)
 *
 *   L-( w , z ) = min c( x ) + w ( l - g( x ) ) + z ( g( x ) - u ) : x \in X
 *               = w l - z  u + min c( x ) + ( z - w ) g( x ) : x \in X   ,
 *
 * where z >= 0 is the Lagrangian multiplier of the constraint g( x ) <= u,
 * while w >= 0 is the Lagrangian multiplier of the constraint l <= g( x ).
 * The Lagrangian dual of (B) is then
 *
 *   (D-)   max  w l - z u + min { c( x ) + ( z - w ) g( x ) : x \in X }
 *               w >= 0  ,  z >= 0
 *
 * Note that in the maximization case things are analogous but different:
 *
 *   L+( w , z ) = max c( x ) + w ( g( x ) - l ) + z ( u - g( x ) ) : x \in X
 *               = z u - w l + max { c( x ) + ( w - z ) g( x ) : x \in X }  ,
 *
 *   (D+)   min  z u - w l + max { c( x ) + ( w - z ) g( x ) : x \in X } 
 *               w >= 0  ,  z >= 0
 *
 * As previously discussed, we will only consider *one* Lagrangian multiplier.
 * This is due to the fact that, at optimality, only one between z and w need
 * be nonzero. The standard set out by RowConstraint is that
 *
 * - for a minimization problem,  y = z - w   (-)
 *
 * - for a maximization problem,  y = w - z   (+)
 *
 * is the value to be written in the dual value of the RowConstraint. In the
 * case l == u (an equality constraint) then (D-) and (D+) become the same
 *
 *   (D-)   max  - y l + min { c( x ) + y g( x ) : x \in X }
 *
 *   (D+)   min  - y l + max { c( x ) + y g( x ) : x \in X } 
 *
 * where y is unconstrained in sign (and beware of the initial "-"). However,
 * in the inequality case things are different. Indeed, let us consider the
 * case of the single inequality constraint g( x ) <= u (l == -INF), whereby
 * then only z is defined: we have
 *
 *   L-( z ) = - z u + min c( x ) + z g( x ) : x \in X   ,
 *
 *   (D-)   max  - z u + min { c( x ) + z g( x ) : x \in X } : z >= 0
 *
 *   L+( z ) = z u + max { c( x ) - z g( x ) : x \in X }  ,
 *
 *   (D+)   min  z u + max { c( x ) - z g( x ) : x \in X } : z >= 0
 *
 * Yet, because of the two different choices (-) and (+), where we take
 * w == 0  ==>  y = z in (-) and y = -z in (+), whence
 *
 *   (D-)   max  - y u + min { c( x ) + y g( x ) : x \in X } : y >= 0
 *
 *   (D+)   min  - y u + max { c( x ) + y g( x ) : x \in X } : y <= 0
 *
 * In the opposite case of the single inequality constraint g( x ) >= l
 * (u == INF) only w is defined and we rather have
 *
 *   L-( w ) = w l + min c( x ) - w g( x ) : x \in X   ,
 *
 *   (D-)   max  w l + min { c( x ) - w g( x ) : x \in X } : w >= 0
 *
 *   L+( w ) = - w l + max { c( x ) + w g( x ) : x \in X }  ,
 *
 *   (D+)   min  - w l + max { c( x ) + w g( x ) : x \in X } : w >= 0
 *
 * Again, due to the difference between (-) and (+), we end up with
 *
 *   (D-)   max  - y l + min { c( x ) + y g( x ) : x \in X } : y <= 0
 *
 *   (D+)   min  - y l + max { c( x ) + y g( x ) : x \in X } : y >= 0
 *
 * To summarise, the Lagrangian Dual always has the form
 *
 *   (D)   max  - y r + min { c( x ) + y g( x ) : x \in X }
 *
 * where r is the non-INF between the two bounds. Then:
 *
 * - for a original maximization problem, corresponding to a minimization
 *   Lagrangian dual (D+):
 *   = for a g( x ) <= u constraint, y = - z <= 0
 *   = for a g( x ) >= l constraint, y =   w >= 0
 *
 * - for a original minimization problem, corresponding to a maximization
 *   Lagrangian dual (D-):
 *   = for a g( x ) <= u constraint, y =   z >= 0
 *   = for a g( x ) >= l constraint, y = - w <= 0
 *
 * Note that LagrangianDualSolver provides a mechanism whereby the inner
 * Solver is only presented with y >= 0 constraints by appropriately changing
 * sign on the data, but this is kept completely hidden from the outside user.
 *
 * After all is said and done, an appropriate Solver is then registered to
 * the Lagrangian Dual Block, and it is used to solve it. The solution is
 * used as the dual solution for (B), while a primal solution is constructed
 * by convexification. Hence, if (B) is *not* a convex program (say, some
 * sub-Block (B_i) has integer variables), then the Lagrangian Dual Block is
 * *not* equivalent to (B) but to its "convexified relaxation", and this is
 * what is solved. */

class LagrangianDualSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 // "import" basic types from Block
 using Index = Block::Index;
 using c_Index = Block::c_Index;

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that LagrangianDualSolver has in addition to these of CDASolver. The 
  * value intLastLDSSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum int_par_type_LDSlv {

 int_LDSlv_iBCopy = intLastParCDAS ,
 ///< if the R3Block has be used for the father block
 
 int_LDSlv_NNMult ,  ///< if Lagrangian multipliers are all >= 0

 int_LDSlv_CloneCfg ,  ///< if BlockSolverConfig need be clone()-d

 int_InnerS_WVarSCfg ,  ///< the Configuration for InnerSolver->get_var_sol

 int_InnerS_WDualSCfg ,  ///< the Configuration for InnerSolver->get_dual_sol

 intPushCostToOwner ,  ///< where the Objective is changed in sub-Block

 intSparseLagPairs ,
 ///< build sparse v_dual_pair for the LagBFunctions
 /**< If nonzero (default 1), in set_Block() each LagBFunction[ h ] is
  * given only the subset of dual pairs whose Lagrangian term
  * LagTerms[ i ][ h ] is non-empty, instead of NumVar dense pairs with
  * empty LinearFunctions for the missing entries. On problems with very
  * sparse coupling (e.g. AC OPF over many time steps) this dramatically
  * reduces setup time, peak memory, and the work the inner Solver has to
  * do per iteration. As a consequence the LagBFunctions may expose
  * different sets of "active" Variables: it is then the responsibility
  * of the inner Solver to handle that — say, by considering each
  * "active" Variable of each LagBFunction as a subset of some "global
  * variable space", which is the union of all of them. The inner Solver
  * may auto-detect the sparse case and switch to a sparse code path, or
  * fall back to the dense one if every LagBFunction happens to expose
  * the full union: in the latter case the sparse setting has no
  * measurable cost.
  *
  * Set to 0 to force the legacy "every LagBFunction sees all multipliers
  * as dense active vars" construction (useful for an inner Solver that
  * cannot handle heterogeneous active sets, or for reproducing the
  * pre-Phase-A behavior for debugging). */

 intLastLDSlvPar   ///< first allowed new int parameter for derived classes
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

  str_LagBF_BCfg ,
  ///< filename of the "default" BlockConfig of the LagBFunction(s)

  str_LagBF_BSCfg ,
  ///< filename of the  "default" BlockSolverConfig of the LagBFunction(s)

  str_LDBlck_BCfg ,   ///< filename of the BlockConfig of the LD

  str_LDBlck_BSCfg ,  ///< filename of the BlockSolverConfig of the LD

  strLastLDSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of string parameters. */

  };  // end( str_par_type_LDSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-int parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-int type that LagrangianDualSolver has in addition to these of
  * CDASolver. The value vintLastLDSlvPar is provided so that the list can
  * be easily further extended by derived classes. */

 enum vint_par_type_LDSlv {
  vint_LDSl_WBCfg = vintLastParCDAS ,
  ///< parameter for associating sub-Block to BlockConfig
  /**< The vector vint_WBCfg, if nonempty, maps the BlockConfig constructed
   * with the filenames out of the parameter vstrCfg [see] into the sub-Block
   * of the Lagrangian Dual; see set_par( std::vector< int > ) for details. */

  vint_LDSl_W2BCfg ,
  ///< parameter telling which sub-Block need be BlockConfig-ured
  /**< The vector vint_W2BCfg, if nonempty, tells which sub-Block of the
   * Lagrangian Dual need be BlockConfig-ured, basically changing the meaning
   * of vint_WBCfg [see] from "dense" to "sparse"; see
   * set_par( std::vector< int > ) for details. */

  vint_LDSl_WBSCfg ,
  ///< parameter for associating sub-Block to BlockSolverConfig
  /**< The vector vint_WBSCfg, if nonempty, maps the BlockSolverConfig
   * constructed with the filenames out of the parameter vstrSCfg [see] into
   * the sub-Block of the Lagrangian Dual; see set_par( std::vector< int > )
   * for details. */

  vint_LDSl_W2BSCfg ,
  ///< parameter telling which sub-Block need be BlockSolverConfig-ured
  /**< The vector vint_W2BCfg, if nonempty, tells which sub-Block of the
   * Lagrangian Dual need be BlockSolverConfig-ured, basically changing the
   * meaning of vint_WBSCfg [see] from "dense" to "sparse"; see
   * set_par( std::vector< int > ) for details. */

  vintWhichPushCost ,  ///< in which sub-Block the Objective is changed
                       /** The vector vintWhichPushCost works in tandem
			* with intPushCostToOwner to decide for which
   * LagBFunction their intPushCostToOwner is (not) activated. */

  vintLastLDSlvPar  ///< first allowed new vector-of-int parameter
                    /**< Convenience value for easily allow derived classes
		     * to extend the set of vector-of-int parameters. */

  };  // end( vint_par_type_LDSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-string parameters
 /** Public enum describing the different parameters of vector-of-string type
  * that LagrangianDualSolver has in addition to these of CDASolver. The
  * value vstrLastLDSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum vstr_par_type_LDSlv {
  vstr_LDSl_Cfg = vstrLastParCDAS ,
  ///< parameter for "the cache of Configurations"

  vstrLastLDSlvPar  ///< first allowed new vector-of-string parameter
                    /**< Convenience value for easily allow derived classes
		     * to extend the set of vector-of-string parameters. */

  };  // end( vstr_par_type_LDSlv )

/** @} ---------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING LagrangianDualSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing LagrangianDualSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 LagrangianDualSolver( void ) : CDASolver() , NumVar( 0 ) , f_nsb( 0 ) ,
  f_max( false ) , LagrDual( nullptr ) , f_BCfg( nullptr ) ,
  f_BSCfg( nullptr ) ,  f_DBCfg( nullptr ) , f_DBSCfg( nullptr ) ,
  f_DBSCfg_map( nullptr ) , static_cons( 0 ) {
  // ensure all parameters are properly given their default value
  iBCopy          = get_dflt_int_par( int_LDSlv_iBCopy );
  NNMult          = get_dflt_int_par( int_LDSlv_NNMult );
  CloneCfg        = get_dflt_int_par( int_LDSlv_CloneCfg );
  WVarSCfg        = get_dflt_int_par( int_InnerS_WVarSCfg );
  WDualSCfg       = get_dflt_int_par( int_InnerS_WDualSCfg );
  PushCostToOwner = get_dflt_int_par( intPushCostToOwner );
  SparseLagPairs  = get_dflt_int_par( intSparseLagPairs );
  ISName          = get_dflt_str_par( str_LDSlv_ISName );
  // all the other string parameters are empty by default, which corresponds
  // to f_BCfg == f_BSCfg == f_DBCfg == f_DBSCfg == nullptr

  // being lazy and not redefining all the other vint parameters, as their
  // default value is empty

  // being lazy and not redefining all the other vdbl parameters, as their
  // default value is empty

  // ensure that the inner Solver is always well defined
  auto ts = new_Solver( ISName );
  InnerSolver = dynamic_cast< CDASolver * >( ts );
  if( ! InnerSolver ) {
   delete ts;
   throw( std::logic_error( ISName + " not a CDASolver" ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the LagrangianDualSolver from the Block

 virtual ~LagrangianDualSolver() { set_Block( nullptr ); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the LagrangianDualSolver has to solve
 /** Set (or changes) the Block that the LagrangianDualSolver has to solve. As
  * customary, set_Block() is a very important method where \p block is
  * thoroughly scanned (after being lock()-ed) and all relevant information is
  * extracted that the LagrangianDualSolver uses; this is even more crucial
  * here, since
  *
  *     INSIDE set_Block() THE LAGRANGIAN DUAL Block IS CONSTRUCTED, POSSIBLY
  *     "EVICTICTING" THE SUB-Block FROM WITHIN THE ORIGINAL Block
  *
  * As a consequence, as customary, generate_abstract_*() are immediately
  * called on \p block since LagrangianDualSolver relies on the existence of
  * a proper Objective in the sub-Block (as this is in turn a requirement for
  * LagBFunction to work) and on proper FRowConstraint to relax, as well as
  * on the *non*-existence of Variable and an Objective in \p block. Hence
  *
  *     THE Block PASSED TO set_Block() MUST BE PROPERLY BlockConfig-URED
  *
  * However, there is a caveat to this. If the sub-Block are *not* evicted
  * but rather R3Block-copied, the copied sub-Block will *not* be
  * BlockConfig-ured even if the original ones were (Block are always born
  * "naked"). Thus, if BlockConfig-uration is needed for the sub-Block, this
  * must be performed via the various ways that LagrangianDualSolver allows
  * for it. However, note that this means
  *
  *     THE BlockConfig-URATION OF THE INDIVIDUAL SUB-Block, RATHER THAN THAT
  *     OF THE WHOLE LAGRANGIAN DUAL Block
  *
  * In fact, the sub-Block are BlockConfig-ured *before* the Lagrangian Dual
  * Block is formed, and generate_abstract_*() is called for them,
  * precisely because LagBFunction needs the relevant pieces of abstract
  * representation to work. All this is irrelevant is sub-Block are evicted,
  * since they will need to be BlockConfig-ured from the start and
  * generate_abstract_*() is (supposedly) called for them when it is for the
  * father \p block. */
 
 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the int parameters of LagrangianDualSolver / the inner Solver
 /** Set the int parameters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details.
  *
  * The parameters handled here are:
  *
  * - int_LDSlv_iBCopy [0]: true (nonzero) if the inner Block of (B) will be
  *   copied (via R3B) when the Lagrangian Dual is created, false (zero) if
  *   they will be moved
  *
  * - int_LDSlv_NNMult [1]: true (nonzero) if the Lagrangian multipliers of
  *   inequality constraints are all constructed as to be non-negative in the
  *   inner Solver and then changed sign, if necessary, when the dual solution
  *   is written in the Block
  *
  * - int_LDSlv_CloneCfg [0]: true (nonzero) if each time a BlockSolverConfig
  *   is apply()-ed to a Block (either the inner Block in a LagBFunction or
  *   the Lagrangian Dual Block itself) it needs to be clone()-d. this is only
  *   necessary if the BlockSolverConfig contains any component (typically,
  *   something in the "extra" Configuration of a ComputeConfig) that gets
  *   "consumed" when apply()-ed, which can happen, but it is not frequent.
  *   it is therefore in general necessary to foresee the possibility of
  *   cloning, but this is not done by default unless this parameter is
  *   properly set (in which case it will apply to *all* BlockSolverConfig,
  *   which may be overkill in some cases but a balance needs to be had).
  *
  * - int_InnerS_WVarSCfg [-1]: the index in the "cache of Configurations"
  *   created with vstr_LDSl_Cfg of the Configuration that is used in the
  *   call of InnerSolver->get_var_solution() to retrieve the dual solution
  *   of the Lagrangian Dual.
  *
  * - int_InnerS_WDualSCfg [-1]: the index in the "cache of Configurations"
  *   created with vstr_LDSl_Cfg of the Configuration that is used in the
  *   call of InnerSolver->get_dual_solution() to retrieve the var solution
  *   of the Lagrangian Dual.
  *
  * - intPushCostToOwner [1]: whether the LagBFunction are instructed (via
  *   the same-named parameter) to change the Objective of all the sub-Block
  *   of the inner Block (in particular, for each variable change the
  *   Objective in the Block in which the variable is defined) as opposed to
  *   changing only the Objective of the "root" inner Block of the
  *   LagBFunction. The parameter works in tandem with vintWhichPushCost
  *   (see) to allow setting intPushCostToOwner differently for each
  *   LagBFunction. If neither this parameter nor vintWhichPushCost are set,
  *   the default behaviour of the LagBFunction (intPushCostToOwner == 1) is
  *   maintained. */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of LagrangianDualSolver / the inner Solver
 /** Set the double parameters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , double value ) override {
  InnerSolver->set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the string parameters of LagrangianDualSolver / the inner Solver
 /** Set the string parameters specific of LagrangianDualSolver, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details.
  *
  * The parameters handled here are:
  *
  * - str_LDSlv_ISName [FakeCDASolver]: the classname used in the Solver
  *   factory to create the inner Solver that actually solves the Lagrangian
  *   Dual. the default FakeCDASolver corresponds to a placeholder that is not
  *   really a suitable choice, and it has to be replaced with a functional one
  *   for LagrangianDualSolver to work, but at least it ensures that
  *   LagrangianDualSolver is not dependent on any other SMS++ module except
  *   the "core" SMS++.
  *
  * - str_LagBF_BCfg [""]: the filename of the "default" BlockConfig of the
  *   inner Block of the LagBFunction(s). If non-empty(), this parameter is
  *   used to create a BlockConfig that is apply()-ed to the inner Block of
  *   all LagBFunction unless specific BlockConfig are provided for that
  *   specific component [see vintWBCfg and vstrBCfg]. if left empty(), no
  *   BlockConfig is apply()-ed unless for those LagBFunction for which
  *   specific ones are provided
  *
  * - str_LagBF_BSCfg [""]: the filename of the "default" BlockSolverConfig
  *   of the inner Block of the LagBFunction(s). If non-empty(), this
  *   parameter is used to create a BlockSolverConfig that is apply()-ed to
  *   the inner Block of all LagBFunction unless specific BlockSolverConfig
  *   are provided for that specific component [see vintWBSCfg and
  *   vstrBSCfg]. if left empty(), no BlockSolverConfig is apply()-ed unless
  *   for those LagBFunction for which specific ones are provided. note that
  *   each LagBFunction does require a working Solver attached to its inner
  *   Block (unless the inner Solver can avoid it for some specially
  *   structured inner Block), so this will have to be provided in some way
  *   (but there are plenty of: besides this parameter and vstrBSCfg, it can
  *   come from the BlockSolverConfig of the whole Lagrangian Dual, see
  *   str_LDBlck_BSCfg, or even from the "extra" Configuration in the
  *   ComputeConfig of LagrangianDualSolver, see set_ComputeConfig())
  *
  * - str_LDBlck_BCfg [""]: the filename of the BlockConfig that is apply()-ed
  *   to the whole Lagrangian Dual Block. This can be used to set the
  *   BlockConfig of the inner Block in the LagBFunction (since the
  *   Lagrangian Dual Block itself has no significant BlockConfig) but it has
  *   to be structured properly, i.e., knowing that the original sub-Block of
  *   the original Block (either itself or a copy) is now set as the inner
  *   Block of a LagBFunction inside the FRealObjective of the corresponding
  *   sub-Block of the Lagrangian Dual Block. doing so allows complete and
  *   fine control on what is done at the cost of being a bit more complex
  *   to handle, which is why the "simpler" str_LagBF_BCfg and vstr_LDSl_Cfg
  *   parameters are provided to spare the user the need to constructing such
  *   a more complex BlockConfig. note that:
  *
  *   = this BlockConfig is ignored if a BlockConfig for the Lagrangian Dual
  *     Block is passed directly as the "extra" Configuration of the
  *     LagrangianDualSolver, see set_ComputeConfig();
  *
  *   = this BlockConfig is only apply()-ed *after* that these created via
  *     str_LagBF_BCfg and vstr_LDSl_Cfg (if any) have been apply()-ed
 *
  * - str_LDBlck_BSCfg [""]: the filename of the BlockSolverConfig that is
  *   apply()-ed to the whole Lagrangian Dual Block. This can be used to set
  *   the BlockSolverConfig of the inner Block in the LagBFunction and
  *   possibly to configure the inner Solver of the LagrangianDualSolver,
  *   although the ComputeConfig of the LagrangianDualSolver itself provides
  *   a more convenient way to do this. in order to BlockSolverConfig-ure
  *   the inner Block such a BlockSolverConfig has to be structured properly,
  *   i.e., knowing that the original sub-Block of the original Block (either
  *   itself or a copy) is now set as the inner Block of a LagBFunction
  *   inside the FRealObjective of the corresponding sub-Block of the
  *   Lagrangian Dual Block. doing so allows complete and fine control on
  *   what is done at the cost of being a bit more complex to handle, which
  *   is why the "simpler" str_LagBF_BSCfg and vstr_LDSl_Cfg parameters
  *   are provided to spare the user the need to constructing such a more
  *   complex BlockSolverConfig. note that:
  *
  *   = this BlockSolverConfig is ignored if a BlockSolverConfig for the
  *     Lagrangian Dual Block is passed directly as the "extra" Configuration
  *     of the LagrangianDualSolver, see set_ComputeConfig();
  *
  *   = this BlockSolverConfig is only apply()-ed *after* that these created
  *     via str_LagBF_BSCfg and vstr_LDSl_Cfg (if any) have been apply()-ed
  *
  *   note that each LagBFunction does require a working Solver attached to
  *   its inner Block (unless the inner Solver can avoid it for some specially
  *   structured inner Block), so this will have to be provided in some way,
  *   this parameter being one of the many */

 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-int parameters of LagrangianDualSolver / inner Solver
 /** Set the vector-of-int parameters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details.
  *
  * The parameters handled here are:
  *
  * - vint_LDSl_WBCfg [empty]: if non-empty(), this maps the BlockConfig
  *   created according to vstr_LDSl_Cfg [see] into the actual inner Block
  *   of the LagBFunction of the sub-Block in the Lagrangian Dual Block. The
  *   exact meaning of this parameter, however, depends on vint_LDSl_W2BCfg
  *   [see]. If the latter is empty, the correspondence is positional:
  *   vint_LDSl_WBCfg[ h ] = k means that the (inner Block of the 
  *   LagBFunction of the) h-th sub-Block will be BlockConfig-ured with the
  *   BlockConfig created using the k-th position in vstr_LDSl_Cfg. If k is
  *   not a valid position into vstr_LDSl_Cfg (it is negative or >=
  *   vstr_LDSl_Cfg.size()), or the corresponding Configuration is not a
  *   BlockConfig, then the h-th sub-Block is not BlockConfig-ured by this
  *   mechanism (but there are plenty of other mechanisms that allows this to
  *   happen, see e.g. str_LagBF_BCfg and str_LDBlck_BCfg). If
  *   vint_LDSl_WBCfg.size() < get_number_nested_Blocks() the sub-Block whose
  *   index exceeds vint_LDSl_WBCfg.size() are ignored, if rather
  *   vint_LDSl_WBCfg.size() > get_number_nested_Blocks() then the excess
  *   Configuration are ignored. If, rather, vint_LDSl_W2BCfg is not empty,
  *   then the correspondence is "sparse":  vint_LDSl_WBCfg[ h ] = k means
  *   that the (inner Block ...) sub-Block vint_LDSl_W2BCfg[ h ] will be
  *   BlockConfig-ured with the BlockConfig created using the k-th position in
  *   vstr_LDSl_Cfg. If k is not a valid position into vstr_LDSl_Cfg (it is
  *   negative or >= vstr_LDSl_Cfg.size()), or the corresponding Configuration
  *   is not a BlockConfig, then the h-th sub-Block is not BlockConfig-ured.
  *   This is done only when information is available in both vectors, i.e.,
  *   for h < min( vint_LDSl_WBCfg.size() , vint_LDSl_W2BCfg.size() ).
  *
  * - vint_LDSl_W2BCfg [empty]: if non-empty(), this provides the list of
  *   sub-Block whose inner Block of the LagBFunction will be BlockConfig-ured
  *   using the Configuration created using vstr_LDSl_Cfg indicated by
  *   vint_LDSl_WBCfg; see the comments to that parameter for details. The
  *   vector needs be ordered in increasing sense and without replications.
  *
  * - vint_LDSl_WBSCfg [empty]: if non-empty(), this maps the
  *   BlockSolverConfig created according to vstr_LDSl_Cfg [see] into the
  *   actual inner Block of the LagBFunction of the sub-Block in the
  *   Lagrangian Dual Block. The exact meaning of this parameter, however,
  *   depends on vint_LDSl_W2BSCfg [see]. If the latter is empty, the
  *   correspondence is positional: vint_LDSl_WBSCfg[ h ] = k means that the
  *   (inner Block of the LagBFunction of the) h-th sub-Block will be
  *   BlockSolverConfig-ured with the BlockSolverConfig created using the
  *   k-th position in vstr_LDSl_Cfg. If k is not a valid position into
  *   vstr_LDSl_Cfg (it is negative or >= vstr_LDSl_Cfg.size()), or the
  *   corresponding Configuration is not a BlockSolverConfig, then the h-th
  *   sub-Block is not BlockSolverConfig-ured by this mechanism (but there
  *   are plenty of other mechanisms that allows this to happen, see e.g.
  *   str_LagBF_BSCfg and str_LDBlck_BSCfg). If vint_LDSl_WBSCfg.size() <
  *   get_number_nested_Blocks() the sub-Block whose index exceeds
  *   vint_LDSl_WBSCfg.size() are ignored, if rather vint_LDSl_WBSCfg.size()
  *   > get_number_nested_Blocks() then the excess Configuration are ignored.
  *   If, rather, vint_LDSl_W2BSCfg is not empty, then the correspondence is
  *   "sparse":  vint_LDSl_WBSCfg[ h ] = k means that the (inner Block ...)
  *   sub-Block vint_LDSl_W2BSCfg[ h ] will be BlockSolverConfig-ured with
  *   the BlockSolverConfig created using the k-th position in vstr_LDSl_Cfg.
  *   If k is not a valid position into vstr_LDSl_Cfg (it is negative or >=
  *   vstr_LDSl_Cfg.size()), or the corresponding Configuration is not a
  *   BlockSolverConfig, then the h-th sub-Block is not
  *   BlockSolverConfig-ured (by this mechanism ...). This is done only when
  *   information is available in both vectors, i.e., for
  *   h < min( vint_LDSl_WBSCfg.size() , vint_LDSl_W2BSCfg.size() ).
  *
  * - vint_LDSl_W2BSCfg [empty]: if non-empty(), this provides the list of
  *   sub-Block whose inner Block of the LagBFunction will be
  *   BlockSolverConfig-ured using the Configuration created using
  *   vstr_LDSl_Cfg indicated by vint_LDSl_WBSCfg; see the comments to that
  *   parameter for details. The vector needs be ordered in increasing sense
  *   and without replications.
  *
  * - vintWhichPushCost [empty]: this works in tandem with intPushCostToOwner
  *   to decide for which of the LagBFunction their own intPushCostToOwner
  *   parameter is set (to a non-default value, i.e., changing only the
  *   objective of the "root" sub-Block of LagBFunction). If not empty(),
  *   vintWhichPushCost[] must contain valid indices of sub-Block: all the
  *   LagBFunction in the given sub-Block are set to the value specified by
  *   intPushCostToOwner, all the others to the opposite value. If
  *   vintWhichPushCost is empty, then *all* the LagBFunction are set to the
  *   value specified by intPushCostToOwner. If neither this parameter nor
  *   intPushCostToOwner are set, the default behaviour of the LagBFunction
  *   (intPushCostToOwner == 1) is maintained. */

 void set_par( idx_type par , std::vector< int > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-double params of LagrangianDualSolver / inner Solver
 /** Set the vector-of-double parameters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , std::vector< double > && value ) override {
  InnerSolver->set_par( par , std::move( value ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-string params of LagrangianDualSolver / inner Solver
 /** Set the vector-of-string parameters specific of LagrangianDualSolver, and
  * allow to directly set those of the inner Solver used to solve the
  * Lagrangian Dual; see the comments to set_ComputeConfig() for details.
  *
  * The parameters handled here are:
  *
  * - vstr_LDSl_Cfg [empty]: if non-empty(), this has to contain a vector of
  *   filenames out of which a vector of Configuration is created; this is a
  *   "chace" of Configurations, that may be BlockConfig, BlockSolverConfig
  *   or any other Configuration, whose elements can then be used for
  *   various purposes; see all times where vstr_LDSl_Cfg is mentioned, such
  *   as str_LagBF_BCfg, str_LDBlck_BCfg, str_LagBF_BSCfg and
  *   str_LDBlck_BSCfg. */

 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*--------------------------------------------------------------------------*/
 /// "translate" an int parameter index of the inner Solver
 /** Takes the index \p par of an int parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type int_par_is( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= intLastParCDAS )
   par += intLastLDSlvPar - intLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a double parameter index of the inner Solver
 /** Takes the index \p par of a double parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type dbl_par_is( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= dblLastParCDAS )
   par += dblLastLDSlvPar - dblLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a string parameter index of the inner Solver
 /** Takes the index \p par of a string parameter of the inner Solver and
  * returns the index that has to be passed to LagrangianDualSolver to have
  * that very same parameter set in the inner Solver; see the comments to
  * set_ComputeConfig() for details. */

 idx_type str_par_is( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= strLastParCDAS )
   par += strLastLDSlvPar - strLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-int parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-int parameter of the inner Solver
  * and returns the index that has to be passed to LagrangianDualSolver to
  * have that very same parameter set in the inner Solver; see the comments
  * to set_ComputeConfig() for details. */

 idx_type vint_par_is( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= vintLastParCDAS )
   par += vintLastLDSlvPar - vintLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-double parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-double parameter of the inner
  * Solver and returns the index that has to be passed to
  * LagrangianDualSolver to have that very same parameter set in the inner
  * Solver; see the comments to set_ComputeConfig() for details. */

 idx_type vdbl_par_is( idx_type par ) const { return( par ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-string parameter index of the inner Solver
 /** Takes the index \p par of a vector-of-string parameter of the inner
  * Solver and returns the index that has to be passed to
  * LagrangianDualSolver to have that very same parameter set in the inner
  * Solver; see the comments to set_ComputeConfig() for details. */

 idx_type vstr_par_is( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= vstrLastParCDAS )
   par += vstrLastLDSlvPar - vstrLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a int parameter index of the LagrangianDualSolver
 /** Takes the index \p par of an int parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type int_par_lds( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= intLastLDSlvPar )
   par -= intLastLDSlvPar - intLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// "translate" a double parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a double parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type dbl_par_lds( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= dblLastLDSlvPar )
   par -= dblLastLDSlvPar - dblLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// "translate" a string parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a string parameter of the LagrangianDualSolver
  * that is actually meant for the inner Solver and returns the value that
  * it would have to be used to set directly in there; see the comments to
  * set_ComputeConfig() for details. */

 idx_type str_par_lds( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= strLastLDSlvPar )
   par -= strLastLDSlvPar - strLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// "translate" a vector-of-int parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-int parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vint_par_lds( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= vintLastLDSlvPar )
   par -= vintLastLDSlvPar - vintLastParCDAS;
  return( par );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// translate a vector-of-double parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-double parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vdbl_par_lds( idx_type par ) const { return( par ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// translate a vector-of-string parameter index of the LagrangianDualSolver
 /** Takes the index \p par of a vector-of-string parameter of the
  * LagrangianDualSolver that is actually meant for the inner Solver and
  * returns the value that it would have to be used to set directly in there;
  * see the comments to set_ComputeConfig() for details. */

 idx_type vstr_par_lds( idx_type par ) const {
  if( par == Inf< idx_type >() )
   return( par );
  if( par >= vstrLastLDSlvPar )
   par -= vstrLastLDSlvPar - vstrLastParCDAS;
  return( par );
  }

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
  * Solver to be directly set by the set_par() of LagrangianDualSolver.
  * 
  * This is done by "translating" all the indices of the parameters of the
  * inner Solver, apart from the standard ones that any CDASolver has, so
  * that they are > than the indices of LagrangianDualSolver parameters (of 
  * the same type). If necessary *_par_is() and *_par_lds() are provided to
  * translate an index of the inner Solver into one of LagrangianDualSolver
  * and vice-versa, respectively, but it is not necessary to use them
  * directly if the Configuration is made via set_ComputeConfig(), as the
  * translation is automatically done by the parameters setting / getting
  * methods. Basically
  *
  *     ALL PARAMETERS OF THE INNER Solver BEHAVE AS IF THEY WERE NATIVE
  *     PARAMETERS OF LagrangianDualSolver WHEN ACCESSED VIA THEIR STRING
  *     NAME
  *
  * Hence, to set the int parameter "intMyParam" of the iner Solver, it is
  * possible to just use
  *
  *     LDS->set_par( LDS->int_par_str2idx( "intMyParam" ) , value );
  *
  * as if one would be accessing the inner Solver directly; this allows all
  * the standard Configuration stuff to work unchanged. When using indices
  * to access parameters, instead, one has to do either
  *
  *     LDS->set_par( LDS->int_par_is( intMyParam ) , value );
  *
  * (assuming the enum value intMyParam to correspond to the string name
  * "intMyParam" as customary) or
  *
  *     LDS->get_inner_Solver()->set_par( intMyParam , value );
  *
  * In both cases one has to know that it is using a LagrangianDualSolver,
  * but this is not a big deal since using of an explicit index (intMyParam)
  * implies compile-time knowledge of the specific solver one is using (in
  * this case, a specific inner Solver inside a LagrangianDualSolver).
  *
  * However, this mechanism has a consequence:
  *
  *     THE PARAMETER INDICES CHANGE MEANING IF THE TYPE OF THE INNER Solver
  *     CHANGES
  *
  * which happens changing its classname ("str_LDSlv_ISName"). Hence, if one
  * wants to change the inner Solver and configure it, it should first do
  * the change and then set the parameters (which is logically required
  * anyway). This is why in set_ComputeConfig() first it is checked if
  * "str_LDSlv_ISName" changes, and only after this is acted upon the
  * standard ThinComputeInterface::set_ComputeConfig() is called to do the
  * bulk of the work.
  *
  * However, set_ComputeConfig() also manages f_extra_Configuration. If not
  * nullptr, the field can contain any amongst
  *
  * - a BlockSolverConfig *
  *
  * - a BlockConfig *
  *
  * - a SimpleConfiguration< std::pair< Configuration * , Configuration * > >
  *   where .first is a BlockSolverConfig * and .second is a
  *   BlockConfig *
  *
  * If any of these are passed, they are apply()-ed to the automatically
  * constructed Lagrangian Dual Block at the moment it is built (which is
  * when set_Block() is called). Note that the structure of the Lagrangian
  * Dual Block is:
  *
  * - a LagrangianDualBlock having a linear Objective (a FRealObjective with
  *   a LinearFunction inside)
  *
  * - The LagrangianDualBlock has exactly as many sub-Block as the original
  *   Block, each of them being an AbstractBlock with the Lagrangian function
  *   of the corresponding sub-Block (a FRealObjective with a LagBFunction
  *   inside, the LagBFunction containing either the original sub-Block or
  *   ita R3B copy).
  *
  * Setting a BlockSolverConfig is a way to register Solver to each
  * sub-Block inside each LagBFunction, which is necessary for the
  * LagBFunction to be able to compute() itself and therefore likely
  * necessary for the inner Solver to work (unless it deals with
  * LagBFunction in specialised ways).
  *
  * Note that the BlockConfig and BlockSolverConfig become property of the
  * LagrangianDualSolver, that clear()-s them and uses them to properly
  * cleanup the Lagrangian Dual Block. If multiple calls to
  * set_ComputeConfig() occur before the LagrangianDualSolver is registered
  * to a Block, all the BlockSolverConfig / BlockConfig of all calls save
  * the last one are lost and have no effect (but at least are properly
  * deleted).
  *
  * This way of setting the Block*Config of the Lagrangian Dual Block takes
  * precedence over doing the same via the str_LDBlck_BCfg and
  * str_LDBlck_BSCfg parameters. */

 void set_ComputeConfig( const ComputeConfig * scfg = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the LagrangianDualSolver log

 void set_log( std::ostream * log_stream = nullptr ) override {
  f_log = log_stream;
  InnerSolver->set_log( f_log );
  /*!!
  for( auto lbf : v_LBF )
   for( auto s : lbf->get_inner_block()->get_registered_solvers() )
    s->set_log( f_log );
    !!*/
  }

/** @} ---------------------------------------------------------------------*/
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void reset_event_handler( int type , EventID id ) override {
  InnerSolver->reset_event_handler( type , id );
  }

/** @} ---------------------------------------------------------------------*/
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

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 double get_elapsed_time( void ) const override {
  return( InnerSolver->get_elapsed_time() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 long get_elapsed_iterations( void ) const override {
  return( InnerSolver->get_elapsed_iterations() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 long get_elapsed_calls( void ) const override {
  return( InnerSolver->get_elapsed_calls() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_lb( void ) override { return( InnerSolver->get_lb() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { return( InnerSolver->get_ub() ); }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override {
  return( InnerSolver->has_dual_solution() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_dual_solution( void ) override {
  return( InnerSolver->has_var_solution() );
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
 /** Write the current solution in the Variable of the Block. This is the
  * solution of the "convexified" problem corresponding to the Lagrangian
  * Dual, and therefore the *dual* solution of the  Lagrangian Dual itself;
  * basically, this is the "important linearization" of the LagBFunction.
  *
  * Since LagrangianDualSolver hinges on the InnerSolver to actually solve
  * the Lagrangian Dual, it must first of all call get_dual_solution() to
  * retrieve the optimal dual solution there. This may require an
  * appropriate Configuration to be passed, which is handled via the
  * int_InnerS_WDualSCfg parameter, giving the index into the "cache of
  * Configuration" constructed via vstr_LDSl_Cfg; if the index provided in
  * this way is invalid (< 0, or larger than the number of Configuration
  * available in the cache), nullptr is used.
  *
  * Note that it is possible for each sub-Block to only generate "a part" of
  * the solution in case this is all the user wants and time/memory can be
  * saved by not dealing with the whole solution. This can be controlled via
  * the BlockConfiguration of each sub-Block, and therefore via the parameters
  * of LagrangianDualSolver that allow to set it. As a consequence, \p solc
  * does not control this aspect.
  *
  * However, \p solc can be used to control which of the sub-Block
  * actually have their solution computed, in case the user is not interested
  * in them all. If *solc is a SimpleConfiguration< std::vector< int > >,
  * it is supposed to be the vector of indices of sub-Block whose variable
  * solution needs be computed (ordered in increasing sense and without
  * replication); the solution of the i-th sub-Block is computed only if i
  * is found anywhere in solc->value.
  *
  * If solc == nullptr, the variable solution is constructed for all the
  * sub-Block. */

 void get_var_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution
 /** Write the current dual solution in the dual values of the RowConstraint
  * of the Block.
  *
  * Since LagrangianDualSolver hinges on the InnerSolver to actually solve
  * the Lagrangian Dual, it must first of all call get_var_solution() to
  * retrieve the optimal var solution there (a primal solution for the
  * original Block, since the Lagrangia Dual is its dual). This may require
  * an appropriate Configuration to be passed, which is handled via the
  * int_InnerS_WVarSCfg parameter, giving the index into the "cache of
  * Configuration" constructed via vstr_LDSl_Cfg; if the index provided in
  * this way is invalid (< 0, or larger than the number of Configuration
  * available in the cache), nullptr is used.
  *
  * The dual solution is actually made by different "pieces":
  *
  * - the dual solution of the relaxed constraints in the father Block;
  *
  * - the dual solution of the constraints inside each individual sub-Block.
  *
  * The former is immediately available from the ColVariable of the Lagrangian
  * Dual, whereas the latter can be written in the sub-Block directly by the
  * Solver that is used to compute() the LagBFunction, provided this is a
  * CDASolver; then, if the sub-Block are a R3B copy, they must be map_back-ed
  * to the original sub-Block.
  *
  * The \p solc Configuration controls which of these pieces is written:
  *
  * - If \p solc is nullptr, then both the dual solution of the relaxed
  *   constraints in the father Block and that of all the sub-Block is
  *   written; the calls to get_dual_solution() of the Solver of the
  *   sub-Block happen with nullptr Configuration (meaning, all of it).
  *
  * - If \p solc is not nullptr, then it can be:
  *
  *   = a pointer to a SimpleConfiguration< std::vector< std::pair< int ,
  *     Configuration * > > >. In this case, the dual solution is only
  *     written for those sub-Block i for which solc->f_value[ h ].first == i
  *     for some h, in which case solc->f_value[ h ].second is passed as the
  *     argument to get_dual_solution() (which can be nullptr). However, the
  *     vector is also used to decide whether the dual solution of the relaxed
  *     constraints in the father Block is written: this happens if and only
  *     if there exist any h such that solc->f_value[ h ].first is invalid
  *     (i.e., either negative or >= get_number_nested_Blocks()), in which
  *     case the Configuration * is ignored. That is, if one only wants the
  *     dual solution of the relaxed constraints, then setting
  *     solc->f_value = { { -1 , nullptr } } does the job.
  *
  *   = a pointer to a SimpleConfiguration< std::vector< std::pair< int ,
  *     int > > >. In this case, the dual solution is only written for those
  *     sub-Block i for which solc->f_value[ h ].first == i for some h, in
  *     which case the argument to get_dual_solution() is the Configuration
  *     in position solc->f_value[ h ].second in the "global cache of
  *     Configuration" created using vstr_LDSl_Cfg. If
  *     solc->f_value[ h ].second is not a valid index in that vector (i.e.,
  *     it is negative or >= vstr_LDSl_Cfg.size()) then nullptr is used.
  *     However, the vector is also used to decide whether the dual solution
  *     of the relaxed constraints in the father Block is written: this
  *     happens if and only if there exist any h such that
  *     solc->f_value[ h ].first is invalid (i.e., either negative or >=
  *     get_number_nested_Blocks()), in which case solc->f_value[ h ].second
  *     is ignored. That is, if one only wants the dual solution of the
  *     relaxed constraints, then setting solc->f_value = { { -1 , 0 } }
  *     does the job.
  *
  *   = a pointer to a SimpleConfiguration< std::vector< Configuration * > >.
  *     In this case, the dual solution is only written for those sub-Block i
  *     with i < solc->f_value.size(), with solc->f_value[ i ] being passed as
  *     the argument to get_dual_solution() (which can be nullptr). However,
  *     the vector is also used to decide whether the dual solution of the
  *     relaxed constraints in the father Block is written: this happens if
  *     and only solc->f_value.size() > get_number_nested_Blocks().
  *
  *   = a pointer to a SimpleConfiguration< std::vector< int > >. In this
  *     case, the dual solution is only written for those sub-Block i with
  *     i < solc->f_value.size(), the argument to get_dual_solution() is the
  *     Configuration in position solc->f_value[ i ] in the "global cache of
  *     Configuration" created using vstr_LDSl_Cfg. If solc->f_value[ i ] is
  *     not a valid index in that vector (i.e., it is negative or >=
  *     vstr_LDSl_Cfg.size()) then nullptr is used. However, the
  *     vector is also used to decide whether the dual solution of the relaxed
  *     constraints in the father Block is written: this happens if and only
  *      solc->f_value.size() > get_number_nested_Blocks().
  *  
  * Note that in order to get the dual solution of the constraints inside a
  * specific sub-Block a CDASolver need be registered there and having been
  * used to compute() the LagBFunction. If there are no Solver registered to
  * the sub-Block, or the one that is used to compute() it (the first one)
  * is not a CDASolver, then the dual solution of the constraints inside that
  * specific sub-Bloc will not be written, and no warning will be issued. Also
  *
  *     IT IS ASSUMED THAT THE LAST TIME THE CDASolver HAS compute()-d THE
  *     SUB-Block IS WITH THE LAGRANGIAN COSTS CORRESPONDING TO THE OPTIMAL
  *     SOLUTION OF THE LAGRANGIAN DUAL.
  *
  * This is automatic if the InnerSolver reports as optimal value the one
  * computed in the very last iteration before it stops, which however may not
  * always be the case because the dual function is nondifferentiable. This
  * issue must be dealt with by properly instructing the InnerSolver to
  * re-compute the dual function in the point it reports as the optimal dual
  * solution in case this does not happen automatically.
  *
  * Note that if the sub-Block are copies, the dual solution will have to be
  * map_back-ed to the originals. For this
  *
  *     THE SAME Configuration IN solc->value[ i ] WILL BE PASSED TO
  *     map_back_solution(), AS THIS SHOULD REASONABLY IDENTIFY THE SAME
  *     SUBSET OF THE DUAL SOLUTION IT DOES IN get_dual_solution().
  *
  * In fact, the "solution Configuration" is the same for both variable and
  * dual solutions. */

 void get_dual_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override {
  return( InnerSolver->new_dual_solution() );
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

/** @} ---------------------------------------------------------------------*/
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

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( int_par_is( InnerSolver->get_num_int_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( dbl_par_is( InnerSolver->get_num_dbl_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( str_par_is( InnerSolver->get_num_str_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vint_par( void ) const override {
  return( vint_par_is( InnerSolver->get_num_vint_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vdbl_par( void ) const override {
  return( vdbl_par_is( InnerSolver->get_num_vdbl_par() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vstr_par( void ) const override {
  return( vstr_par_is( InnerSolver->get_num_vstr_par() ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  static const std::array dflt_int_par = {
    0 , // int_LDSlv_iBCopy
    1 , // int_LDSlv_NNMult
    0 , // int_LDSlv_CloneCfg
   -1 , // int_InnerS_WVarSCfg
   -1 , // int_InnerS_WDualSCfg
    1 , // intPushCostToOwner
    1 , // intSparseLagPairs (default on; sparse path is bit-equivalent
        //                   to dense and unlocks DoEasy=1 in LDS)
   };

  if( ( par >= intLastParCDAS ) && ( par < intLastLDSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );

  return( InnerSolver->get_dflt_int_par( int_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  return( InnerSolver->get_dflt_dbl_par( dbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par ) const
  override {
  static const std::array< std::string , 5 > dflt_str_par = {
   "FakeCDASolver" ,  // str_LDSlv_ISName
   "" ,               // str_LagBF_BCfg
   "" ,               // str_LagBF_BSCfg
   "" ,               // str_LDBlck_BCfg
   ""                 // str_LDBlck_BSCfg
   };

  if( ( par >= strLastParCDAS ) && ( par < strLastLDSlvPar ) )
   return( dflt_str_par[ par - strLastParCDAS ] );

  return( InnerSolver->get_dflt_str_par( str_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< int > & get_dflt_vint_par( idx_type par )
  const override {
  static const std::vector< int > _empty;
  if( ( par == vint_LDSl_WBCfg ) || ( par == vint_LDSl_W2BCfg ) ||
      ( par == vint_LDSl_WBSCfg ) || ( par == vint_LDSl_W2BSCfg ) ||
      ( par == vintWhichPushCost ) )
   return( _empty );
  else
   return( InnerSolver->get_dflt_vint_par( vint_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< double > & get_dflt_vdbl_par( idx_type par )
  const override {
  return( InnerSolver->get_dflt_vdbl_par( vdbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< std::string > & get_dflt_vstr_par(
					     idx_type par ) const override {
  static const std::vector< std::string > _empty;
  if( par == vstr_LDSl_Cfg )
   return( _empty );
  else
  return( InnerSolver->get_dflt_vstr_par( vstr_par_lds( par ) ) );
  }

/*--------------------------------------------------------------------------*/
 
 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( int_LDSlv_iBCopy ):     return( iBCopy );
   case( int_LDSlv_NNMult ):     return( NNMult );
   case( int_LDSlv_CloneCfg ):   return( CloneCfg );
   case( int_InnerS_WVarSCfg ):  return( WVarSCfg );
   case( int_InnerS_WDualSCfg ): return( WDualSCfg );
   case( intPushCostToOwner ):   return( PushCostToOwner );
   case( intSparseLagPairs ):    return( SparseLagPairs );
   }

  return( InnerSolver->get_int_par( int_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  return( InnerSolver->get_dbl_par( dbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override {
  switch( par ) {
   case( str_LDSlv_ISName ): return( ISName );
   case( str_LagBF_BCfg ):   return( LagBF_BCfg );
   case( str_LagBF_BSCfg ):  return( LagBF_BSCfg );
   case( str_LDBlck_BCfg ):  return( LDBlck_BCfg );
   case( str_LDBlck_BSCfg ): return( LDBlck_BSCfg );
   }

  return( InnerSolver->get_str_par( str_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override {
  switch( par ) {
   case( vint_LDSl_WBCfg ):   return( WBCfg );
   case( vint_LDSl_W2BCfg ):  return( W2BCfg );
   case( vint_LDSl_WBSCfg ):  return( WBSCfg );
   case( vint_LDSl_W2BSCfg ): return( W2BSCfg );
   case( vintWhichPushCost ): return( WhichPushCost );
   }

  return( InnerSolver->get_vint_par( vint_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< double > & get_vdbl_par( idx_type par )
  const override {
  return( InnerSolver->get_vdbl_par( vdbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< std::string > & get_vstr_par( idx_type par )
  const override  {
  if( par == vstr_LDSl_Cfg )
   return( FCfg );

  return( InnerSolver->get_vstr_par( vstr_par_lds( par ) ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > int_pars_map = {
   { "int_LDSlv_iBCopy"     , int_LDSlv_iBCopy } ,
   { "int_LDSlv_NNMult"     , int_LDSlv_NNMult } ,
   { "int_LDSlv_CloneCfg"   , int_LDSlv_CloneCfg } ,
   { "int_InnerS_WVarSCfg"  , int_InnerS_WVarSCfg } ,
   { "int_InnerS_WDualSCfg" , int_InnerS_WDualSCfg } ,
   { "intPushCostToOwner"   , intPushCostToOwner } ,
   { "intSparseLagPairs"    , intSparseLagPairs }
   };

  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );

  return( int_par_is( InnerSolver->int_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  return( dbl_par_is( InnerSolver->dbl_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > str_pars_map = {
   { "str_LDSlv_ISName" , LagrangianDualSolver::str_LDSlv_ISName } ,
   { "str_LagBF_BCfg"   , LagrangianDualSolver::str_LagBF_BCfg } ,
   { "str_LagBF_BSCfg"  , LagrangianDualSolver::str_LagBF_BSCfg } ,
   { "str_LDBlck_BCfg"  , LagrangianDualSolver::str_LDBlck_BCfg } ,
   { "str_LDBlck_BSCfg" , LagrangianDualSolver::str_LDBlck_BSCfg }
   };

  const auto it = str_pars_map.find( name );
  if( it != str_pars_map.end() )
   return( it->second );

  return( str_par_is( InnerSolver->str_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > vint_pars_map = {
   { "vint_LDSl_WBCfg"   , vint_LDSl_WBCfg } ,
   { "vint_LDSl_W2BCfg"  , vint_LDSl_W2BCfg } ,
   { "vint_LDSl_WBSCfg"  , vint_LDSl_WBSCfg } ,
   { "vint_LDSl_W2BSCfg" , vint_LDSl_W2BSCfg } ,
   { "vintWhichPushCost" , vintWhichPushCost }
   };

  const auto it = vint_pars_map.find( name );
  if( it != vint_pars_map.end() )
   return( it->second );

  return( vint_par_is( InnerSolver->vint_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vdbl_par_str2idx( const std::string & name )
  const override {
  return( vdbl_par_is( InnerSolver->vdbl_par_str2idx( name ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override {
  if( name == "vstr_LDSl_Cfg" )
   return( vstr_LDSl_Cfg );

  return( vstr_par_is( InnerSolver->vstr_par_str2idx( name ) ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 7 > int_pars_str = {
   "int_LDSlv_iBCopy" , "int_LDSlv_NNMult" , "int_LDSlv_CloneCfg" ,
   "int_InnerS_WVarSCfg" , "int_InnerS_WDualSCfg" , "intPushCostToOwner" ,
   "intSparseLagPairs" };

  if( ( idx >= intLastParCDAS ) && ( idx < intLastLDSlvPar ) )
   return( int_pars_str[ idx - intLastParCDAS ] );

  return( InnerSolver->int_par_idx2str( int_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {
  return( InnerSolver->dbl_par_idx2str( dbl_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 5 > str_pars_str = {
   "str_LDSlv_ISName" , "str_LagBF_BCfg" , "str_LagBF_BSCfg" ,
   "str_LDBlck_BCfg" , "str_LDBlck_BSCfg" };

  if( ( idx >= strLastParCDAS ) && ( idx < strLastLDSlvPar ) )
   return( str_pars_str[ idx - strLastParCDAS ] );

  return( InnerSolver->str_par_idx2str( str_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 5 > vint_pars_str = {
   "vint_LDSl_WBCfg"  , "vint_LDSl_W2BCfg" ,
   "vint_LDSl_WBSCfg" , "vint_LDSl_W2BSCfg" , "vintWhichPushCost" };

  if( ( idx >= vintLastParCDAS ) && ( idx < vintLastLDSlvPar ) )
   return( vint_pars_str[ idx - vintLastParCDAS ] );

  return( InnerSolver->vint_par_idx2str( vint_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vdbl_par_idx2str( idx_type idx )
  const override {
  return( InnerSolver->vdbl_par_idx2str( vdbl_par_lds( idx ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override {
  static const std::string _vstr_LDSl_Cfg = "vstr_LDSl_Cfg";
  if( idx == vstr_LDSl_Cfg )
   return( _vstr_LDSl_Cfg );

  return( InnerSolver->vstr_par_idx2str( vstr_par_lds( idx ) ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING THE State OF THE LagrangianDualSolver -------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the State of the LagrangianDualSolver
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns the current "internal state" of the inner Solver
 /** If this LagrangianDualSolver has an inner Solver, then it simply
  * returns the State of the inner Solver. Otherwise, it returns nullptr. */

 State * get_State( void ) const override {
  if( InnerSolver )
   return( InnerSolver->get_State() );
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// sets the current "internal state" of the inner Solver
 /** If this LagrangianDualSolver has an inner Solver, then this
  * method sets the State of the inner Solver by simply calling
  * put_State() of the inner Solver and passing the given State
  * on. Otherwise, it does nothing. */

 void put_State( const State & state ) override {
  if( InnerSolver )
   InnerSolver->put_State( state );
  }

/*--------------------------------------------------------------------------*/
 /// sets the current "internal state" of the iiner Solver
 /** If this LagrangianDualSolver has an inner Solver, then this
  * method sets the State of the inner Solver by simply calling
  * put_State() of the inner Solver and passing the given State
  * on. Otherwise, it does nothing. */

 void put_State( State && state ) override {
  if( InnerSolver )
   InnerSolver->put_State( std::move( state ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** Small private "copy" of AbstractBlock, for the Block representing the
 * Lagrangian Dual that LagrangianDualSolver constructs.
 *
 * This does almost eactly everything that AbstractBlock does, save for a
 * minor detail: in add_Modifiation(), if the Modification is sent to a
 * channel that is not defined, no exception is raised (as Block does) but
 * the Modification is simply reverted to the default channel.
 *
 * The reason for this is the "eviction mode" in which sub-Block are "stolen"
 * from the original Block that LagrangianDualSolver is registered to. In
 * this case, an issue can occur if a channel is opened in that Block, and
 * a Modification from a sub-Block is sent to the channel: since LagBFunction
 * keeps the channel information, the Modification arrives to the Lagrangian
 * Dual Block on a channel that it does not "know", causing the exception to
 * be raised.
 */
 
 class LagrangianDualBlock : public AbstractBlock {

 public:

  explicit LagrangianDualBlock( Block * father = nullptr ) :
           AbstractBlock( father ) {}

  ~LagrangianDualBlock() override {};

  void add_Modification( sp_Mod mod , ChnlName chnl ) override;
  };


/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void register_inner_Solver( void ) {
  if( ! LagrDual )
   return;

  InnerSolver->set_id( this );
  LagrDual->register_Solver( InnerSolver );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void unregister_inner_Solver( void ) {
  if( InnerSolver && LagrDual ) {
   InnerSolver->set_id();
   LagrDual->unregister_Solver( InnerSolver );
   }
  }

/*--------------------------------------------------------------------------*/

 void clear_LD_BlockSolverConfig( bool keepcfg = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void clear_LD_BlockConfig( bool keepcfg = false );

/*--------------------------------------------------------------------------*/

 void clear_inner_BlockSolverConfig( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void clear_inner_BlockConfig( void );

/*--------------------------------------------------------------------------*/
 /** Returns the index as active variable of the LagBFunction of the given
  * FRowConstraint, be it static or dynamic.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variable, Int< Index >() if
  *         \p con does not correspond to any FRowConstraint */

 Index index_of_constraint( const FRowConstraint * con ) {
  auto i = index_of_static_constraint( con );
  if( i < Inf< Index >() )
   return( i );

  return( index_of_dynamic_constraint( con ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the index as active variable of the LagBFunction of the given
  * static FRowConstraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variable, Int< Index >() if
  *         \p con does not correspond to any static FRowConstraint */

 Index index_of_static_constraint( const FRowConstraint * con );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the index as active variable of the LagBFunction of the given
  * dynamic FRowConstraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding index as active variable, Int< Index >() if
  *         \p con does not correspond to any dynamic FRowConstraint */

 Index index_of_dynamic_constraint( const FRowConstraint * con );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to) FRowConstraint corresponding to a the active
  * variable of the LagBFunction with the given index, be it static or
  * dynamic.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint
  */
 
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to the) static FRowConstraint corresponding to the
  * active variable of the LagBFunction with the given index.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint
  */

 FRowConstraint * static_constraint_with_index( Index i );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Returns the (pointer to the) static FRowConstraint corresponding to the
  * active variable of the LagBFunction with the given index.
  *
  * @param i the index of an active variable of the LagBFunction
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if \p i doesn't correspond to a constraint
  */

 FRowConstraint * dynamic_constraint_with_index( Index i ) {
  return( idx_to_dcon[ i - static_cons ] );
  }

/*--------------------------------------------------------------------------*/
 // return the index of the given sub-Block

 Index Block2Index( Block * blck ) {
  for( ;; ) {
   if( auto lbf = dynamic_cast< LagBFunction * >( blck ) ) {
    blck = lbf->get_Observer()->get_Block();
    continue;
    }
   auto fb = blck->get_f_Block();
   if( fb == LagrDual ) {
    auto it = std::lower_bound( blck_to_idx.begin() , blck_to_idx.end() ,
				blck_int( static_cast< AbstractBlock * >(
							      blck ) , 0 ) ,
				[]( const auto & a , const auto & b ) {
				 return( a.first < b.first ); } );
    return( it->second );
    }
   blck = fb;
   if( ! blck )
    throw( std::invalid_argument(
	        "LagrangianDualSolver: Variable belonging to wrong Block" ) );
   }
  }

/*--------------------------------------------------------------------------*/
/* If NNMult == true, the Lagrangian multipliers of inequality constraints
 * are all constructed as to be non-negative. This means that:
 *
 * - for a maximization original problem (f_max == true) a <= constraint
 *   is passed as it is in the Lagrangian term, while a >= need be changed
 *   sign (both all the coefficients and the RHS/LHS)
 *
 * - for a minimization original problem (f_max == false) a >= constraint
 *   is passed as it is in the Lagrangian term, while a <= need be changed
 *   sign (both all the coefficients and the RHS/LHS) */

 bool to_be_reversed( const FRowConstraint & con );

/*--------------------------------------------------------------------------*/

 double constr2val( const FRowConstraint & con , ColVariable & lvar );

/*--------------------------------------------------------------------------*/
 /// maps the dual box of con, if any, onto the box of the multiplier lvar
 /** If con declares a box for its optimal Lagrangian multiplier, writes in
  * box the corresponding box for the Lagrangian variable lvar and returns
  * true; returns false otherwise. */

 bool constr2box( const FRowConstraint & con , const ColVariable & lvar ,
		  std::array< ColVariable::VarValue , 2 > & box );

/*--------------------------------------------------------------------------*/

 void split_constraint( const FRowConstraint & con ,
			std::vector< LinearFunction::v_coeff_pair > & split );

/*--------------------------------------------------------------------------*/

 void set_PushCostToOwner( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters- - - - - - - - - - - - - - - - - - - - - - - - - -

 bool iBCopy;         ///< true if the R3Block conversion has to be done

 bool NNMult;         ///< true if Lagrangian multipliers are all >= 0

 bool CloneCfg;       ///< true if BlockSolverConfig need be clone()-d

 int WVarSCfg;        ///< the Configuration for IS->get_var_solution()

 int WDualSCfg;       ///< the Configuration for IS->get_dual_solution()

 bool PushCostToOwner;  ///< how to set the same-named LagBFunction parameter

 bool SparseLagPairs;
 ///< true if dual pairs with empty Lagrangian term are skipped in set_Block

 std::string ISName;  ///< classname of the inner Solver

 std::string LagBF_BCfg;
 ///< the filename for the BlockConfig of the individual LagBFunction

 std::string LagBF_BSCfg;
 ///< the filename for the BlockSolverConfig of individual LagBFunction

 std::string LDBlck_BCfg;  ///< the filename for the BlockConfig of the LD

 std::string LDBlck_BSCfg;
 ///< the filename for the BlockSolverConfig of the LD

 std::vector< int > WBCfg;    ///< map between sub-Block and BlockConfig

 std::vector< int > W2BCfg;   ///< which sub-Block to BlockConfig

 std::vector< int > WBSCfg;   ///< map between sub-Block and BlockSolverConfig

 std::vector< int > W2BSCfg;  ///< which sub-Block to BlockSolverConfig

 std::vector< int > WhichPushCost;
 ///< for which sub-Block change PushCostToOwner

 std::vector< std::string > FCfg;  ///< filenames for Configurations
 
 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumVar;      ///< (current) number of variables

 Index f_nsb;       ///< number of sub-Block

 bool f_max;        ///< true if (B) was a max problem, false otherwise

 LagrangianDualBlock * LagrDual;
                    ///< the automatically constructed Lagrangian Dual

 CDASolver * InnerSolver;   ///< the Solver attached to LagrDual

 BlockConfig * f_BCfg;      ///< the BlockConfig for LagrDual

 BlockSolverConfig * f_BSCfg;   ///< the BlockSolverConfig for LagrDual

 BlockConfig * f_DBCfg;      ///< the default individual BlockConfig

 BlockSolverConfig * f_DBSCfg;   ///< the default individual BlockSolverConfig

 /// "meta" form of the default individual BlockSolverConfig
 /** If str_LagBF_BSCfg points to a meta BlockSolverConfig (a
  * SimpleConfiguration< std::map< std::string , Configuration * > > mapping a
  * Block classname() to the BlockSolverConfig for it), it is stored here
  * instead of f_DBSCfg and dispatched per inner Block by classname(). This
  * allows a single str_LagBF_BSCfg to configure inner Blocks of different
  * types (e.g. ThermalUnitBlock and HydroSystemUnitBlock) at once. Exactly one
  * of f_DBSCfg / f_DBSCfg_map is non-null (or both null). */
 SimpleConfiguration< std::map< std::string , Configuration * > > * f_DBSCfg_map;

 /// the default individual BlockSolverConfig for inner Block \p inner
 /** Returns the BlockSolverConfig to use as the "default" for the inner Block
  * \p inner: the per-classname() entry of f_DBSCfg_map if a meta config was
  * given, else f_DBSCfg; nullptr if none applies. Defined out-of-line since it
  * dynamic_cast<>s to the (here incomplete) BlockSolverConfig. */
 BlockSolverConfig * default_BSCfg_for( Block * inner ) const;

 std::vector< Configuration * > v_Cfg;  ///< the "Configuration cache"

 std::vector< UpdateSolver * > v_US;   /// the UpdateSolvers

 std::vector< LagBFunction * > v_LBF;  /// the LagBFunction

 // dictionaries- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index static_cons;  ///< number of static constraints
                     /** Total number of static constraints in the Block: the
		      * first static_cons active variables in the LagBFunction
		      * will never change. */

 /* The following vectors are used in order to keep track between the
  * [FRow]Constraints of the Block and the Lagrangian variables of all the
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

 typedef std::pair< AbstractBlock * , Index > blck_int;

 /** This is a vector of pairs that contain the pointer of a sub-Block of
  * the LagrangianDual and its index into the vector of sub-Block; this is
  * ordered by Block. */

 std::vector< blck_int > blck_to_idx;  ///< from Block * to index

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

 void cleanup_LagrDual( bool keepcfg = false );

/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/

 void flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod );

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class LagrangianDualSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* LagrangianDualSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File LagrangianDualSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
