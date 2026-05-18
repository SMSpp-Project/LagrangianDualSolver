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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PrimalProximalHeur.h"

#include "ColVariable.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using p_FRO = FRealObjective *;
using p_LF  = LinearFunction *;
using p_DQF = DQuadFunction *;

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

 if( f_log && ( logVerb >= 2 ) )
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
  // zero-clamp tiny linear coefficients to keep the proximal term clean
  auto register_binvar = [ & ]( ColVariable * pv ) {
   double c1 = 0;
   double c2 = 0;
   const auto i_in_obj = fobj->is_active( pv );
   if( i_in_obj < fobj->get_num_active_var() ) {
    if( ! fobj->is_linear() ) {
     auto qf = static_cast< p_DQF >( fobj );
     c1 = qf->get_linear_coefficient( i_in_obj );
     c2 = qf->get_quadratic_coefficient( i_in_obj );
     }
    else
     c1 = static_cast< p_LF >( fobj )->get_coefficient( i_in_obj );
    if( std::abs( c1 ) < 1e-6 )
     c1 = 0;
    }
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

 if( NumStatVar == 0 )
  throw( std::invalid_argument(
   "PrimalProximalHeur::initialize: no static binary Variable" ) );

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
  default:
   LagrangianDualSolver::set_par( par , value );
  }
 }  // end( PrimalProximalHeur::set_par( double ) )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int PrimalProximalHeur::compute( bool changedvars )
{
 lock();  // lock the Solver mutex

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // helper: try to insert a new candidate solution (with the given objective
 // value) into the v_best_sol heap, maintaining the heap invariant and
 // updating best_bound / worst_bound. Used both by the inner-Solver event
 // handler (which discovers feasible-integer solutions at a Bundle SS) and
 // by the main loop (which checks at the end of each PPH iteration).

 auto record_feasible = [ this ]( double value ) {
  // is the new value better than the current best?
  const bool better = f_max ? ( value > best_bound ) : ( value < best_bound );
  if( better )
   best_bound = value;

  // is the new value worse than the current worst?
  const bool worse = f_max ? ( value < worst_bound ) : ( value > worst_bound );

  if( v_best_sol.size() < f_MaxSol ) {
   // there is free space, just throw the new Solution in and re-heap
   v_best_sol.emplace_back( f_Block->get_Solution() , value );
   std::push_heap( v_best_sol.begin() , v_best_sol.end() );
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

  *bad = { f_Block->get_Solution() , value };
  // re-heap from scratch (the replacement may have broken the order)
  std::make_heap( v_best_sol.begin() , v_best_sol.end() ,
                  [ this ]( const sol_value & a , const sol_value & b ) {
                   return( f_max ? ( a.second < b.second )
                                 : ( a.second > b.second ) );
                   } );

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

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "PrimalProximalHeur::compute: NumStatVar = " << NumStatVar
         << std::endl;

 best_bound = f_max ? - Inf< double >() : Inf< double >();

 std::vector< double > sol( NumStatVar );

 // process pending Modification before solving - - - - - - - - - - - - - - -

 const bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "PrimalProximalHeur::compute: unable to lock the Block" ) );

 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );

 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index iters = 0;
 bool is_the_same = false;
 bool is_integer  = false;
 int  res         = kOK;

 while( true ) {

  if( f_log && ( logVerb >= 2 ) )
   *f_log << std::endl << "PrimalProximalHeur::compute: iteration "
          << iters << std::endl;

  // at iteration 0 seed the "previous solution" with a random 0/1 vector
  if( iters == 0 )
   for( Index kvar = 0 ; kvar < NumStatVar ; ++kvar )
    sol[ kvar ] = rand() % 2;
/*
  if( iters == 0 ){
          Index index = 0;
	        Index kvar = 0;
          for( const auto & sbi : f_Block->get_nested_Blocks() ) {
            for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ){
	             sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
              auto si = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
              kvar++;
            }
            index++;
          }
        }
*/

  previous_sol.insert( previous_sol.begin() , sol.begin() , sol.end() );

  // apply the proximal term to the inner objective(s)- - - - - - - - - - - -

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: adding penalty terms"
          << std::endl;

  add_penalty_terms();

  // register an event handler on the inner Solver that checks at every
  // iteration whether the current point is feasible-and-integer, and if so
  // tries to record it in v_best_sol via the record_feasible() helper

  Index index_event = 0;

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
      if( f_log && ( logVerb >= 2 ) )
       *f_log << "  (event) IS_NOT_INTEGER_SOL" << std::endl;
      return( ThinComputeInterface::eContinue );
      }
     if( f_log && ( logVerb >= 2 ) )
      *f_log << "  (event) IS_INTEGER_SOL, LB = "
             << ( InnerSolver->get_lb() - addterm ) << std::endl;

     if( f_Block->is_feasible() ) {
      if( f_log && ( logVerb >= 2 ) )
       *f_log << "  (event) IS_FEASIBLE_SOL: " << value_FUNCTION
              << std::endl;
      record_feasible( value_FUNCTION );
      }
     else if( f_log && ( logVerb >= 2 ) )
      *f_log << "  (event) IS_INFEASIBLE_SOL: " << value_FUNCTION
             << std::endl;

     return( ThinComputeInterface::eContinue );
     } );

  // delegate to the inner Solver - - - - - - - - - - - - - - - - - - - - - -

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: solving Lagrangian Dual"
          << std::endl;

  res = InnerSolver->compute( changedvars );

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: Lagrangian Dual solved"
          << std::endl;

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

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: removing penalty terms"
          << std::endl;

  remove_penalty_terms();

  // recompute the linearized penalty contribution at this iteration- - - - -

  penalty = 0;
  addterm = 0;
  if( iters >= 1 )
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

  if( f_log && ( logVerb >= 2 ) ) {
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
  // 1e-3 in penalty terms, equivalent too identical to within sqrt( 1e-3 / R ))
  if( iters >= 1 )
   is_the_same = ( penalty < 1e-3 );

  #ifdef BIN_VARS
   is_integer = is_integer_solution();
   if( f_log && ( logVerb >= 2 ) )
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
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "  IS_FEASIBLE_SOL" << std::endl;
   record_feasible( value_FUNCTION );
   }
  else if( f_log && ( logVerb >= 2 ) )
   *f_log << "  IS_INFEASIBLE_SOL" << std::endl;

  ++iters;

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
  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: INFEASIBLE, iters = "
          << ( iters - 1 ) << std::endl;
  unlock();
  return( res );
  }
 if( res == kInfeasible ) {
  res = kUnbounded;
  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: UNBOUNDED, iters = "
          << ( iters - 1 ) << std::endl;
  unlock();
  return( res );
  }

 // if we stopped because of solution convergence, also record the current
 // (integer) value as a candidate Solution so it shows up in v_best_sol
 if( is_the_same ) {
  double bound;
  if( ! is_integer )
   bound = f_max ? get_lb() : get_ub();
  else {
   bound = f_max ? std::max( get_lb() , value_FUNCTION )
                 : std::min( get_ub() , value_FUNCTION );
   v_best_sol.emplace_back( f_Block->get_Solution() , value_FUNCTION );
   }
  best_bound = bound;
  if( f_log && ( logVerb >= 2 ) )
   *f_log << "PrimalProximalHeur::compute: converged after " << ( iters - 1 )
          << " iterations, LB = " << InnerSolver->get_lb()
          << ", UB = " << InnerSolver->get_ub()
          << ", best feasible = " << bound << std::endl;
  }

 changed_penalties = false;

 unlock();  // unlock the Solver mutex

 return( res );

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

 if( f_log && ( logVerb >= 2 ) )
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
  const auto mp   = Observer::make_par( eModBlck , chnl );
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
  const auto mp   = Observer::make_par( eModBlck , chnl );
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
