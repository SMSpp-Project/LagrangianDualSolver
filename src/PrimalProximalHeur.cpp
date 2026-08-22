/*--------------------------------------------------------------------------*/
/*--------------------- File PrimalProximalHeur.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework as a Lagrangian-based
 * primal-proximal heuristic for binary optimization problems.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

#define BIN_VARS 1
/**< Compile-time switch: when defined, PrimalProximalHeur expects all the
 * static "active" Variables of the inner sub-Blocks to be binary (0/1),
 * and exploits this by linearising the quadratic proximal term into a
 * simple update of the linear coefficient. When undefined the proximal
 * term is left quadratic. */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define PrimalProximalHeur_LOG 0
/* If non-zero, enables the verbose trace messages emitted by
 * PrimalProximalHeur to f_log under runtime logVerb control. Default 0
 * keeps the trace silent; set to 1 manually during development to follow
 * the heuristic step by step. */

#if PrimalProximalHeur_LOG
 #define LOG_VERB( lvl ) if( f_log && ( logVerb >= ( lvl ) ) )
#else
 #define LOG_VERB( lvl ) if constexpr ( false )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PrimalProximalHeur.h"

#include "ColVariable.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "BlockSolverConfig.h"

#include <cmath>

#include <filesystem>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using p_FRO = FRealObjective *;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register PrimalProximalHeur to the Solver factory

SMSpp_insert_in_factory_cpp_0( PrimalProximalHeur );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF PrimalProximalHeur -----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_Block( Block * block )
{
 if( f_Block )           // was already attached to a Block
  guts_of_destructor();  // release the PPH-specific state first

 LagrangianDualSolver::set_Block( block );  // call the base method

 if( f_Block ) {         // a new Block is now attached
  initialize();
  best_bound  =   f_max ? - Inf< double >() : Inf< double >();
  worst_bound = - best_bound;
  addterm     = 0;
  }
 }  // end( PrimalProximalHeur::set_Block )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::initialize( void )
{
 // count and check the static "active" binary ColVariable, and meanwhile
 // build the per-sub-Block dictionaries (idx_to_var_sbi1/2) carrying the
 // linear and quadratic coefficients of the inner objective Function of
 // each sub-Block. The dictionaries are later read by add_penalty_terms()
 // and remove_penalty_terms() to apply / strip the proximal term.

 LOG_VERB( 2 )
  *f_log << "PrimalProximalHeur::initialize: maxIter = " << maxIter
         << ", R = " << R << ", logVerb = " << logVerb << std::endl;

 NumStatVar = 0;

 const auto n_sub = f_Block->get_number_nested_Blocks();
 pos_id_sbi.resize( n_sub );
 idx_to_var_sbi1.resize( n_sub );
 idx_to_var_sbi2.resize( n_sub );
 Funct_sbi.resize( n_sub );
 Funct_sbi_quad.resize( n_sub );
 is_linear.resize( n_sub );

 Index index = 0;
 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  // identify the kind (linear / quadratic) of the inner objective Function
  // and cache a copy of it for later evaluation via get_funct_value(); see
  // the C05Function hierarchy in C05Function.h for the actual taxonomy
  const auto fobj = static_cast< Function * >(
                        static_cast< p_FRO >( sbi->get_objective()
                                              )->get_function() );
  if( const auto qobj = dynamic_cast< p_DQF >( fobj ) ) {
   is_linear[ index ]      = false;
   Funct_sbi_quad[ index ] = *qobj;
   }
  else {
   is_linear[ index ] = true;
   Funct_sbi[ index ] = *static_cast< p_LF >( fobj );
   }

  pos_id = 0;

  // helper: register a single binary ColVariable into the sub-Block's
  // dictionaries — read the coefficients off the inner objective Function,
  // zero-clamp tiny linear coefficients to keep the proximal term clean.
  // Variables not currently present in the inner objective are skipped:
  // injecting them in add_penalty_terms() via add_variable() would issue
  // a FunctionModVarsAddd that some leaf Blocks (e.g. ThermalUnitBlock)
  // cannot handle (their obj has a fixed physical layout). The proximal
  // penalty is therefore not applied to those Variables.
  auto register_binvar = [ & ]( ColVariable * pv ) {
   const auto i_in_obj = fobj->is_active( pv );
   if( i_in_obj >= fobj->get_num_active_var() )
    return;  // not in inner obj: skip to avoid unsupported add_variable
   double c1 = 0;
   double c2 = 0;
   if( ! fobj->is_linear() ) {
    auto qf = static_cast< p_DQF >( fobj );
    c1 = qf->get_linear_coefficient( i_in_obj );
    c2 = qf->get_quadratic_coefficient( i_in_obj );
    }
   else
    c1 = static_cast< p_LF >( fobj )->get_coefficient( i_in_obj );
   if( std::abs( c1 ) < 1e-6 )
    c1 = 0;
   idx_to_var_sbi1[ index ].emplace_back( c1 , pv );
   idx_to_var_sbi2[ index ].emplace_back( c2 , pv );
   ++pos_id;
   ++NumStatVar;
   };

  // helper: check the binary nature of a Variable in BIN_VARS mode
  auto is_binary = []( const ColVariable & v ) {
   #ifdef BIN_VARS
    return( v.is_integer() && ( v.get_lb() == 0.0 ) && ( v.get_ub() == 1.0 ) );
   #else
    return( true );
   #endif
   };

  // scan static Variables of sbi (Singles / Vectors / Multiarrays)
  for( const auto & el : sbi->get_static_variables() ) {

   // single ColVariable
   if( un_any_thing_0( ColVariable , el , {
                       if( is_binary( var ) )
                        register_binvar( &var );
                       } ) )
    continue;

   // vector of ColVariable
   if( un_any_thing_1( ColVariable , el , {
                       for( Index j = 0 ; j < var.size() ; ++j ) {
                        const auto pv = var.data() + j;
                        if( is_binary( *pv ) )
                         register_binvar( pv );
                        }
                       } ) )
    continue;

   // multi-array of ColVariable
   if( un_any_thing_K( ColVariable , el , {
                       for( Index j = 0 ; j < var.num_elements() ; ++j ) {
                        const auto pv = var.data() + j;
                        if( is_binary( *pv ) )
                         register_binvar( pv );
                        }
                       } ) )
    continue;
   }

  pos_id_sbi[ index ] = pos_id;
  ++index;
  }

 if( NumStatVar == 0 ) LOG_VERB( 1 )
  *f_log << "PrimalProximalHeur::initialize: no static binary Variable, "
            "compute() will delegate to LagrangianDualSolver" << std::endl;

 }  // end( PrimalProximalHeur::initialize )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIterPPH ): maxIter = value; break;
  case( intMaxSol ):     f_MaxSol = value; break;
  case( intLogVerb ):
   logVerb = value & 3;
   LagrangianDualSolver::set_par( par , std::max( 0 , value >> 2 ) );
   break;
  default:
   LagrangianDualSolver::set_par( par , value );
  }
 }  // end( PrimalProximalHeur::set_par( int ) )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dbl_penaltyFactor ): R = value; break;

  // the accuracy required of the heuristic is kept here, while
  // dbl_LDSRelAcc is the one the inner Solver has to solve the Lagrangian
  // Dual with; before, both were the latter
  case( dblRelAcc ):     RelAcc = value; break;
  case( dbl_LDSRelAcc ): InnerSolver->set_par( dblRelAcc , value ); break;
  default:
   LagrangianDualSolver::set_par( par , value );
  }
 }  // end( PrimalProximalHeur::set_par( double ) )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

bool PrimalProximalHeur::gap_closed( void ) const
{
 // the optimum lies between the valid bound of the unpenalized iteration
 // and the value of the best feasible solution found, so the heuristic has
 // delivered the required accuracy as soon as the two are that close
 const double lb = f_max ? best_bound : valid_bound;
 const double ub = f_max ? valid_bound : best_bound;

 if( ( ! std::isfinite( lb ) ) || ( ! std::isfinite( ub ) ) )
  return( false );

 return( ub - lb <= RelAcc * std::max( double( 1 ) , std::abs( lb ) ) );

 }  // end( PrimalProximalHeur::gap_closed )

/*--------------------------------------------------------------------------*/

int PrimalProximalHeur::compute( bool changedvars )
{
 // no static binary Variable to apply the proximal penalty to: PPH has
 // nothing to add over the inner Lagrangian Dual, fall back to it
 if( NumStatVar == 0 )
  return( LagrangianDualSolver::compute( changedvars ) );

 lock();  // lock the Solver mutex

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // helper: try to insert a new candidate solution (with the given objective
 // value) into the v_best_sol heap, maintaining the heap invariant and
 // updating best_bound / worst_bound. Used both by the inner-Solver event
 // handler (which discovers feasible-integer solutions at a Bundle SS) and
 // by the main loop (which checks at the end of each PPH iteration).

 // comparator making the heap top the best Solution: a is "worse" than b
 // in the sense of the optimization direction
 auto sol_cmp = [ this ]( const sol_value & a , const sol_value & b ) {
  return( f_max ? ( a.second < b.second ) : ( a.second > b.second ) );
  };

 auto record_feasible = [ this , & sol_cmp ]( double value ) {
  // is the new value better than the current best?
  const bool better = f_max ? ( value > best_bound ) : ( value < best_bound );
  if( better )
   best_bound = value;

  // is the new value worse than the current worst?
  const bool worse = f_max ? ( value < worst_bound ) : ( value > worst_bound );

  if( v_best_sol.size() < f_MaxSol ) {
   // there is free space, just throw the new Solution in and re-heap
   v_best_sol.emplace_back( f_Block->get_Solution( nullptr , false ) ,
                            value );
   std::push_heap( v_best_sol.begin() , v_best_sol.end() , sol_cmp );
   if( worse )
    worst_bound = value;
   return;
   }

  if( worse )           // no space, and worse than worst: discard
   return;

  // no space, but better than the worst: find the worst entry and replace
  // it; meanwhile compute the second-worst to refresh worst_bound after
  double second_worst = f_max ? Inf< double >() : - Inf< double >();
  auto bad = v_best_sol.end();
  for( auto it = v_best_sol.begin() ; it != v_best_sol.end() ; ++it )
   if( it->second == worst_bound )
    bad = it;
   else {
    const bool s_worse = f_max ? ( it->second < second_worst )
                               : ( it->second > second_worst );
    if( s_worse )
     second_worst = it->second;
    }

  if( bad == v_best_sol.end() )  // no entry at worst_bound (numerical drift)
   return;                       // play it safe and discard the new value

  delete bad->first;
  *bad = { f_Block->get_Solution( nullptr , false ) , value };
  // re-heap from scratch (the replacement may have broken the order)
  std::make_heap( v_best_sol.begin() , v_best_sol.end() , sol_cmp );

  // the worst is now either the second-worst-before or the new value
  worst_bound = f_max ? std::min( second_worst , value )
                      : std::max( second_worst , value );
  };

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // helper: tell whether the current values of the sub-Block Variables are
 // all 0/1 to within 1e-6 (the binary integrality tolerance). Used to gate
 // the feasibility check in the main loop and in the event handler.

 auto is_integer_solution = [ this ]() -> bool {
  const auto n_sub = f_Block->get_number_nested_Blocks();
  for( Index index = 0 ; index < n_sub ; ++index )
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto si = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
    if( si * ( 1 - si ) > 1e-6 )
     return( false );
    }
  return( true );
  };

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LOG_VERB( 2 )
  *f_log << "PrimalProximalHeur::compute: NumStatVar = " << NumStatVar
         << std::endl;

 best_bound = f_max ? - Inf< double >() : Inf< double >();
 valid_bound = f_max ? Inf< double >() : - Inf< double >();

 std::vector< double > sol( NumStatVar );

 // process pending Modification before solving - - - - - - - - - - - - - - -

 const bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "PrimalProximalHeur::compute: unable to lock the Block" ) );

 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );

 // warm-start: solve the LP relaxation of f_Block with an auxiliary
 // :MILPSolver picked from WarmStartCfg.txt via the Configuration factory,
 // so that the dual variables (and hence the initial Lambda multipliers
 // of the inner Lagrangian Dual) start from a meaningful point. The
 // concrete Solver (CPLEX / Gurobi / SCIP / HiGHS) is selected by the
 // config file, not hard-wired in code.

 LOG_VERB( 2 )
  *f_log << "PrimalProximalHeur::compute: solving MILP relaxation"
         << std::endl;

 auto warmstart = new_aux_solver( "WarmStartCfg.txt" );
 warmstart->set_Block( f_Block );
 warmstart->compute( changedvars );
 warmstart->get_dual_solution();
 warmstart->get_var_solution();

 delete warmstart;

 // re-sync the Lagrangian multipliers with the just-computed LP duals: the
 // Lambda Variables of the Lagrangian Dual Block were initialised from the
 // Constraint duals when set_Block() built it, i.e., before the warm-start
 // ran, so without this the warm-start would not reach the inner Solver
 {
  auto Ls = LagrDual->get_static_variable_v< ColVariable >( "Lambda_s" );
  auto Lsit = Ls->begin();
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el ,
                        [ & Lsit ]( FRowConstraint & con ) {
                         ( Lsit++ )->set_value( con.get_dual() );
                         } ,
                        un_any_type< FRowConstraint >() );

  auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );
  auto Ldit = Ld->begin();
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_dynamic( el ,
                         [ & Ldit ]( FRowConstraint & con ) {
                          ( Ldit++ )->set_value( con.get_dual() );
                          } ,
                         un_any_type< FRowConstraint >() );
  }

 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index iters = 0;
 bool is_the_same = false;
 bool is_integer  = false;
 int  res         = kOK;

 while( true ) {

  LOG_VERB( 2 )
   *f_log << std::endl << "PrimalProximalHeur::compute: iteration "
          << iters << std::endl;

  // read the current sub-Block Variable values: they are the proximal
  // center of this iteration
  {
   const auto n_sub = f_Block->get_number_nested_Blocks();
   Index kvar = 0;
   for( Index index = 0 ; index < n_sub ; ++index )
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
     sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
     ++kvar;
     }
   }

  previous_sol.assign( sol.begin() , sol.end() );

  // apply the proximal term to the inner objective(s)- - - - - - - - - - - -
  // with R == 0 the proximal machinery is inert: the penalty terms would
  // rewrite every objective coefficient to its own value (a storm of
  // Modification per inner iteration via the event handler below, heavy at
  // scale) and the handler could only re-record the same un-penalised
  // points, so both are skipped and PrimalProximalHeur degenerates into a
  // warm-started LagrangianDualSolver (plus the final primal recovery)
  //
  // the first iteration is not penalized either: it solves the Lagrangian
  // Dual of the original objective, which is the only valid bound on the
  // original problem the heuristic can produce, every penalized iteration
  // bounding the penalized objective instead. get_lb() / get_ub() report
  // it and gap_closed() measures the solution found against it. The
  // iteration is not lost, since the proximal center of the next one is
  // the fractional solution the Lagrangian Dual converges to, which is
  // what the penalty is meant to be built on
  const bool penalized = ( R != 0 ) && ( iters > 0 );

  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: adding penalty terms"
          << std::endl;

  if( penalized )
   add_penalty_terms();

  // register an event handler on the inner Solver that checks at every
  // iteration whether the current point is feasible-and-integer, and if so
  // tries to record it in v_best_sol via the record_feasible() helper

  Index index_event = 0;

  if( penalized )
   LagrangianDualSolver::set_event_handler(
    ThinComputeInterface::eEverykIteration ,
    [ this , &sol , &index_event , &record_feasible , &is_integer_solution ]
    () {
     if( index_event != InnerSolver->get_elapsed_iterations() - 1 )
      return( ThinComputeInterface::eContinue );
     ++index_event;

     value_FUNCTION = get_funct_value();
     add_penalty_terms();

     // refresh sol[] with the latest sub-Block Variable values
     const auto n_sub = f_Block->get_number_nested_Blocks();
     Index kvar = 0;
     for( Index index = 0 ; index < n_sub ; ++index )
      for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
       sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
       ++kvar;
       }

     // recompute the linearized penalty contribution
     addterm = 0;
     penalty = 0;
     for( Index ivar = 0 ; ivar < NumStatVar ; ++ivar ) {
      const auto si  = sol[ ivar ];
      const auto psi = previous_sol[ ivar ];
      penalty += R * ( psi - si ) * ( psi - si );
      addterm += R * si * ( 1.0 - 2.0 * psi );
      }

     if( ! is_integer_solution() ) {
      LOG_VERB( 2 )
       *f_log << "  (event) IS_NOT_INTEGER_SOL" << std::endl;
      return( ThinComputeInterface::eContinue );
      }
     LOG_VERB( 2 )
      *f_log << "  (event) IS_INTEGER_SOL, LB = "
             << ( InnerSolver->get_lb() - addterm ) << std::endl;

     if( f_Block->is_feasible() ) {
      LOG_VERB( 2 )
       *f_log << "  (event) IS_FEASIBLE_SOL: " << value_FUNCTION
              << std::endl;
      double rec_cost;
      if( recover_primal( rec_cost ) )
       record_feasible( rec_cost );
      }
     else LOG_VERB( 2 )
      *f_log << "  (event) IS_INFEASIBLE_SOL: " << value_FUNCTION
             << std::endl;

     return( ThinComputeInterface::eContinue );
     } );

  // delegate to the inner Solver - - - - - - - - - - - - - - - - - - - - - -

  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: solving Lagrangian Dual"
          << std::endl;

  res = InnerSolver->compute( changedvars );

  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: Lagrangian Dual solved"
          << std::endl;

  // the first iteration solved the Lagrangian Dual of the original
  // objective: its bound is the valid one on the original problem
  if( ! penalized )
   valid_bound = f_max ? InnerSolver->get_ub() : InnerSolver->get_lb();

  // read back the current sub-Block Variable values - - - - - - - - - - - -
  {
   const auto n_sub = f_Block->get_number_nested_Blocks();
   Index kvar = 0;
   for( Index index = 0 ; index < n_sub ; ++index )
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
     sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
     ++kvar;
     }
   }

  // strip the proximal term from the inner objective(s)- - - - - - - - - - -

  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: removing penalty terms"
          << std::endl;

  if( penalized )
   remove_penalty_terms();

  // recompute the linearized penalty contribution at this iteration- - - - -

  penalty = 0;
  addterm = 0;
  for( Index ivar = 0 ; ivar < NumStatVar ; ++ivar ) {
   const auto si  = sol[ ivar ];
   const auto psi = previous_sol[ ivar ];
   penalty += R * ( psi - si ) * ( psi - si );
   addterm += R * si * ( 1.0 - 2.0 * psi );
   }

  // bound from the inner Solver, corrected by the proximal contribution
  const auto value_bound = f_max ? ( InnerSolver->get_ub() - addterm )
                                 : ( InnerSolver->get_lb() - addterm );

  value_FUNCTION = get_funct_value();

  LOG_VERB( 2 ) {
   if( std::abs( value_bound - value_FUNCTION ) /
       std::max( std::abs( value_bound ) , std::abs( value_FUNCTION ) )
       >= 1e-3 )
    *f_log << "  ERROR: bound mismatch = " << ( value_bound - value_FUNCTION )
           << std::endl;
   *f_log << "  iter = "          << iters                   << std::endl
          << "  InnerSolver LB = "<< InnerSolver->get_lb()   << std::endl
          << "  InnerSolver UB = "<< InnerSolver->get_ub()   << std::endl
          << "  penalty = "       << penalty                 << std::endl
          << "  addterm = "       << addterm                 << std::endl
          << "  corrected bound = " << value_bound           << std::endl
          << "  objective value = " << value_FUNCTION        << std::endl
          << "  best value = "    << ( f_max ? get_lb() : get_ub() )
          << std::endl;
   }

  // stopping criterion: previous and current solution coincide (within
  // 1e-3 in penalty terms, i.e., identical to within sqrt( 1e-3 / R ))
  is_the_same = ( penalty < 1e-3 );

  #ifdef BIN_VARS
   is_integer = is_integer_solution();
   LOG_VERB( 2 )
    *f_log << ( is_integer ? "  IS_INTEGER_SOL" : "  IS_NOT_INTEGER_SOL" )
           << std::endl;
  #endif

  // try to record the current solution if it is feasible (and integer in
  // BIN_VARS mode) and we are past the first iteration
  #ifdef BIN_VARS
   const bool can_record = f_Block->is_feasible() && is_integer && ( iters >= 1 );
  #else
   const bool can_record = f_Block->is_feasible() && ( iters >= 1 );
  #endif
  if( can_record ) {
   LOG_VERB( 2 )
    *f_log << "  IS_FEASIBLE_SOL" << std::endl;
   double rec_cost;
   if( recover_primal( rec_cost ) )
    record_feasible( rec_cost );
   }
  else LOG_VERB( 2 )
   *f_log << "  IS_INFEASIBLE_SOL" << std::endl;

  ++iters;

  // the accuracy required of the heuristic has been reached: the best
  // feasible solution found is within RelAcc of the optimum, since the
  // latter lies between valid_bound and it, so there is nothing left to
  // look for
  if( gap_closed() )
   break;

  if( is_the_same || ( iters >= maxIter ) )
   break;

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if the inner Block is not an R3B-copy, restore the original inner
 // objectives "as if nothing had happened"
 if( ! iBCopy )
  for( auto lbf : v_LBF )
   lbf->cleanup_inner_objective();

 // since the inner Solver solves the dual of the original Block, swap
 // unbounded and infeasible return codes back to the original sense
 if( res == kUnbounded ) {
  res = kInfeasible;
  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: INFEASIBLE, iters = "
          << ( iters - 1 ) << std::endl;
  unlock();
  return( res );
  }
 if( res == kInfeasible ) {
  res = kUnbounded;
  LOG_VERB( 2 )
   *f_log << "PrimalProximalHeur::compute: UNBOUNDED, iters = "
          << ( iters - 1 ) << std::endl;
  unlock();
  return( res );
  }

 // final primal recovery: a Lagrangian point in general violates the coupling
 // constraints of f_Block. Now that the inner objectives have been restored to
 // the original ones, fix the (rounded) proximal binaries and solve the
 // restricted problem on f_Block, whose optimum enforces the coupling
 // constraints and is therefore a genuinely feasible completion with its true
 // primal cost, hence a valid bound. This is done unconditionally: value_FUNCTION
 // is the Lagrangian value of the point (a lower bound), never a valid upper
 // bound. If the restricted problem is infeasible the point is discarded.
 {
  double rec_cost;
  if( recover_primal( rec_cost ) )
   record_feasible( rec_cost );
  }

 LOG_VERB( 2 )
  *f_log << "PrimalProximalHeur::compute: "
         << ( is_the_same ? "converged" : "stopped" ) << " after "
         << ( iters - 1 ) << " iterations, LB = " << get_lb()
         << ", UB = " << get_ub()
         << ", best feasible = " << ( f_max ? get_lb() : get_ub() )
         << std::endl;

 changed_penalties = false;

 unlock();  // unlock the Solver mutex

 // an error of the inner Solver is an error of the heuristic; otherwise
 // what is returned says whether the solution found is as accurate as it
 // was required to be, which is what an heuristic can promise: kOK if the
 // gap between it and the valid bound is within RelAcc, kLowPrecision if a
 // solution was found but no such guarantee comes with it
 if( res >= kError )
  return( res );

 if( gap_closed() )
  return( kOK );

 return( kLowPrecision );

 }  // end( PrimalProximalHeur::compute )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::process_outstanding_Modification( void )
{
 // PrimalProximalHeur does not chain to LagrangianDualSolver::process_-
 // outstanding_Modification(); it scans the locally queued Modification
 // and decides whether the sub-Block dictionaries need to be rebuilt
 // (any change in the inner objectives invalidates the linear / quadratic
 // coefficient cache) and/or whether the stored best Solutions need to
 // be re-checked for feasibility (any change in the inner constraints
 // may invalidate them).

 bool reload            = false;
 bool check_feasibility = false;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire the lock, spin on failure

 for( auto & mod : v_mod ) {
  for( const auto & sbi : f_Block->get_nested_Blocks() ) {
   if( mod->get_Block() != sbi )
    continue;
   if( dynamic_cast< C05FunctionModLin * >( mod.get() ) ) {
    // inner-objective coefficient change: dictionaries are stale
    reload = true;
    }
   else {
    // inner-constraint change: stored Solutions may be infeasible
    check_feasibility = true;
    reload            = true;
    }
   }
  }

 LOG_VERB( 2 )
  *f_log << "PrimalProximalHeur::process_outstanding_Modification: "
         << "reload = " << reload
         << ", check_feasibility = " << check_feasibility << std::endl;

 if( reload && ( ! changed_penalties ) ) {
  guts_of_destructor();
  initialize();
  }

 if( check_feasibility ) {
  // rebuild v_best_sol keeping only the Solutions that are still feasible
  std::vector< sol_value > v_best_sol_new;
  v_best_sol_new.reserve( v_best_sol.size() );

  for( auto & sol : v_best_sol ) {
   sol.first->write( f_Block );
   if( f_Block->is_feasible() )
    v_best_sol_new.push_back( sol );
   else
    delete sol.first;
   }

  v_best_sol = std::move( v_best_sol_new );
  std::make_heap( v_best_sol.begin() , v_best_sol.end() ,
                  [ this ]( const sol_value & a , const sol_value & b ) {
                   return( f_max ? ( a.second < b.second )
                                 : ( a.second > b.second ) );
                   } );
  }

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release the lock

 }  // end( PrimalProximalHeur::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

CDASolver * PrimalProximalHeur::new_aux_solver( const std::string & cfgname )
{
 // resolve the config file next to this source file: __FILE__ is baked in
 // at compile time and points to LagrangianDualSolver/src/<this>.cpp, so
 // the sibling LagrangianDualSolver/<cfgname> is reachable regardless of
 // the process working directory
 const std::string cfgfile =
  ( std::filesystem::path( __FILE__ ).parent_path().parent_path() / cfgname
    ).string();

 auto cfg = Configuration::deserialize( cfgfile );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ( ! bsc ) || bsc->get_SolverNames().empty() ) {
  delete cfg;
  throw( std::runtime_error( "PrimalProximalHeur::new_aux_solver: " +
                             cfgfile + " is not a valid BlockSolverConfig"
                           ) );
  }

 auto slvr = dynamic_cast< CDASolver * >(
                            Solver::new_Solver( bsc->get_SolverName( 0 ) ) );
 if( ! slvr ) {
  delete bsc;
  throw( std::runtime_error( "PrimalProximalHeur::new_aux_solver: the Solver"
                             " in " + cfgfile + " is not a CDASolver" ) );
  }

 if( bsc->num_ComputeConfig() > 0 )
  if( auto cc = bsc->get_SolverConfig( 0 ) )
   slvr->set_ComputeConfig( cc );

 delete bsc;
 return( slvr );

 }  // end( PrimalProximalHeur::new_aux_solver )

/*--------------------------------------------------------------------------*/

bool PrimalProximalHeur::recover_primal( double & cost )
{
 // fix the proximal binaries at their rounded values; eNoMod keeps the
 // fixing invisible to the Solver attached to f_Block, and it is undone
 // below before anyone else can compute()
 for( const auto & sbd : idx_to_var_sbi1 )
  for( const auto & dv : sbd ) {
   const auto pv = dv.second;
   pv->set_value( std::round( pv->get_value() ) );
   pv->is_fixed( true , eNoMod );
   }

 // solve the restricted problem with the recovery Solver, which sees the
 // fixed binaries as bounds and enforces the coupling constraints
 auto solve_restricted = [ & ]( const char * stage ) -> bool {
  auto recovery = new_aux_solver( "RecoveryCfg.txt" );

  Index index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobj = static_cast< Function * >(
                  static_cast< p_FRO >( sbi->get_objective()
                                        )->get_function() );

  if( ! is_linear[ index ] ) {
   auto qf = static_cast< p_DQF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    const auto c2 = idx_to_var_sbi2[ index ][ ivar ].first;
    if( qf->get_num_active_var() > idx_in_obj )
     #ifdef BIN_VARS
      qf->modify_linear_coefficient( idx_in_obj , c1 , mp );
     #else
      qf->modify_term( idx_in_obj , c1 , c2 , mp );
     #endif
    }
   }
  else {
   auto lf = static_cast< p_LF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    if( lf->get_num_active_var() > idx_in_obj )
     lf->modify_coefficient( idx_in_obj , c1 , mp );
    }
   }

  sbi->close_channel( chnl );
  ++index;
  }

  recovery->set_Block( f_Block );

  const auto rc = recovery->compute( true );
  const bool ok = ( rc >= kOK ) && ( rc < kError ) &&
                  recovery->has_var_solution();
  if( ok ) {
   recovery->get_var_solution();  // the completion into the Block Variables

   // evaluate the recovered solution on the pristine copies of the inner
   // objectives: the live ones only had the true costs restored on the
   // binary Variables, so their value is not the true cost of the solution
   double value = 0;
   Index idx = 0;
   for( const auto & sbi : f_Block->get_nested_Blocks() ) {
    if( is_linear[ idx ] ) {
     Funct_sbi[ idx ].compute( true );
     value += Funct_sbi[ idx ].get_value();
     }
    else {
     Funct_sbi_quad[ idx ].compute( true );
     value += Funct_sbi_quad[ idx ].get_value();
     }
    ++idx;
    }
   cost = value;
   }

  LOG_VERB( 2 ) {
   *f_log << "  recover_primal[ " << stage << " ]: ";
   if( ok )
    *f_log << "cost = " << cost << std::endl;
   else
    *f_log << "infeasible" << std::endl;
   }

  delete recovery;
  return( ok );
  };

 bool ok = solve_restricted( "full" );

 if( ! ok ) {
  // the fully-fixed point admits no feasible completion: relax the
  // restriction one-sidedly, keeping only the binaries at 1 fixed, so
  // more can be activated (e.g. more units committed to cover a demand
  // the fixed set cannot); the problem remains a restriction of the
  // original one, so any of its solutions still yields a valid bound.
  // Note that with equality couplings this direction does not always
  // help (an over-active set can be as infeasible as an under-active
  // one), whence the last resort below.
  for( const auto & sbd : idx_to_var_sbi1 )
   for( const auto & dv : sbd )
    if( dv.second->get_value() < 0.5 )
     dv.second->is_fixed( false , eNoMod );

  ok = solve_restricted( "ones" );
  }

 // un-fix the proximal binaries
 for( const auto & sbd : idx_to_var_sbi1 )
  for( const auto & dv : sbd )
   dv.second->is_fixed( false , eNoMod );

 if( ! ok )
  // last resort: nothing fixed, i.e. the original problem within the
  // recovery Solver's own budget; whatever incumbent it finds is still
  // a valid bound, only no longer tied to the proximal point
  ok = solve_restricted( "free" );

 return( ok );

 }  // end( PrimalProximalHeur::recover_primal )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::add_penalty_terms( void )
{
 // add the proximal penalty term to the (quadratic) inner objective Function
 // of every sub-Block, using the dictionaries built in initialize(). The
 // update is bundled inside an Observer channel so all the Modification land
 // atomically as one GroupModification.

 Index pos   = 0;
 Index index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobj = static_cast< Function * >(
                  static_cast< p_FRO >( sbi->get_objective()
                                        )->get_function() );

  if( ! is_linear[ index ] ) {
   // quadratic inner objective: update both linear and (BIN_VARS off only)
   // quadratic coefficients
   auto qf = static_cast< p_DQF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    const auto c2 = idx_to_var_sbi2[ index ][ ivar ].first;
    if( qf->get_num_active_var() > idx_in_obj ) {
     #ifdef BIN_VARS
      qf->modify_linear_coefficient(
                              idx_in_obj ,
                              c1 + R * ( 1.0 - 2.0 * previous_sol[ pos ] ) ,
                              mp );
     #else
      qf->modify_term( idx_in_obj ,
                       c1 - R * 2.0 * previous_sol[ pos ] ,
                       c2 + R ,
                       mp );
     #endif
     }
    else if( R > 0 )
     #ifdef BIN_VARS
      qf->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
                        R * ( 1.0 - 2.0 * previous_sol[ pos ] ) ,
                        0.0 ,
                        mp );
     #else
      qf->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
                        - R * 2.0 * previous_sol[ pos ] ,
                        R ,
                        mp );
     #endif
    ++pos;
    }
   }
  else {
   // linear inner objective: only the linear coefficient changes
   auto lf = static_cast< p_LF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    if( lf->get_num_active_var() > idx_in_obj )
     lf->modify_coefficient(
                       idx_in_obj ,
                       c1 + R * ( 1.0 - 2.0 * previous_sol[ pos ] ) ,
                       mp );
    else if( R > 0 )
     lf->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
                       R * ( 1.0 - 2.0 * previous_sol[ pos ] ) ,
                       mp );
    ++pos;
    }
   }

  sbi->close_channel( chnl );
  ++index;
  }
 }  // end( PrimalProximalHeur::add_penalty_terms )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::remove_penalty_terms( void )
{
 // restore the original linear (and quadratic, when BIN_VARS is undefined)
 // coefficients of every sub-Block's inner objective Function, using the
 // dictionaries built in initialize(). Bundled inside an Observer channel
 // like add_penalty_terms().

 Index index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobj = static_cast< Function * >(
                  static_cast< p_FRO >( sbi->get_objective()
                                        )->get_function() );

  if( ! is_linear[ index ] ) {
   auto qf = static_cast< p_DQF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    const auto c2 = idx_to_var_sbi2[ index ][ ivar ].first;
    if( qf->get_num_active_var() > idx_in_obj )
     #ifdef BIN_VARS
      qf->modify_linear_coefficient( idx_in_obj , c1 , mp );
     #else
      qf->modify_term( idx_in_obj , c1 , c2 , mp );
     #endif
    }
   }
  else {
   auto lf = static_cast< p_LF >( fobj );
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    const auto idx_in_obj = fobj->is_active(
                              idx_to_var_sbi1[ index ][ ivar ].second );
    const auto c1 = idx_to_var_sbi1[ index ][ ivar ].first;
    if( lf->get_num_active_var() > idx_in_obj )
     lf->modify_coefficient( idx_in_obj , c1 , mp );
    }
   }

  sbi->close_channel( chnl );
  ++index;
  }
 }  // end( PrimalProximalHeur::remove_penalty_terms )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::guts_of_destructor( void )
{
 for( auto & el : v_best_sol )
  delete el.first;
 v_best_sol.clear();

 idx_to_var_sbi1.clear();
 idx_to_var_sbi2.clear();
 Funct_sbi.clear();
 Funct_sbi_quad.clear();
 is_linear.clear();
 previous_sol.clear();
 pos_id_sbi.clear();

 }  // end( PrimalProximalHeur::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*------------------- End File PrimalProximalHeur.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
