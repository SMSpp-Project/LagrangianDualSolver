/*--------------------------------------------------------------------------*/
/*--------------------- File PrimalProximalHeur.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework.
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

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PrimalProximalHeur.h"

#include "FRealObjective.h"

#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using p_FRO = FRealObjective *;
using p_LF = LinearFunction *;
using p_DQF = DQuadFunction *;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

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
 if( f_Block )           // was attached to a Block
  guts_of_destructor();  // cleanup
 LagrangianDualSolver::set_Block( block );  // call method of base class
 if( f_Block ) {
  initialize();
  best_bound = f_max ? - Inf< double >() : Inf< double >();
  worst_bound = - best_bound;
  addterm = 0.0;
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::initialize( void )
{
 // count and check the static ColVariable - - - - - - - - - - - - - - - - - 
 // meanwhile construct the static dictionaries for the linear and quadratic
 // terms of the (quadratic) objective functions of the sub-blocks

 NumStatVar = 0;

 pos_id_sbi.resize( f_Block->get_number_nested_Blocks() );
 idx_to_var_sbi1.resize( f_Block->get_number_nested_Blocks() );
 idx_to_var_sbi2.resize( f_Block->get_number_nested_Blocks() );
 Funct_sbi.resize( f_Block->get_number_nested_Blocks());
 
 Index index = 0;
 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
   Funct_sbi[index] = static_cast< Function * >(static_cast< p_FRO >( sbi->get_objective())->get_function());
   pos_id = 0;
   double addval1;
   double addval2;
   for( const auto & el : sbi->get_static_variables() ) {
   // Singles
   if( un_any_thing_0( ColVariable , el ,
		  {
      #ifdef BIN_VARS 
      if( var.is_integer() && var.get_lb() == 0.0 && var.get_ub() == 1.0 ){
      #endif
        if(static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(&var) ){ 
          auto indexz = static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function())->is_active(&var);
          auto fobj_sbi = static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function());
          if( ! fobj_sbi->is_linear() ){
            addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
            addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
          } else {
            addval1 = static_cast< p_LF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_coefficient(indexz);
          }
          idx_to_var_sbi1[index].push_back(double_var( addval1 , &var ));
          idx_to_var_sbi2[index].push_back(double_var( addval2 , &var ));
        } else {
          idx_to_var_sbi1[index].push_back(double_var( 0.0 , &var ));
          idx_to_var_sbi2[index].push_back(double_var( 0.0 , &var ));
        }
        pos_id++;
			  NumStatVar++;
        #ifdef BIN_VARS 
        } 
        #endif
        } ) )
     continue;
   // Vectors
   if( un_any_thing_1( ColVariable , el ,
		  {
			for (Index j = 0 ; j < var.size() ; ++j){
        #ifdef BIN_VARS 
        if( (var.data()+j)->is_integer() && (var.data()+j)->get_lb() == 0.0 && (var.data()+j)->get_ub() == 1.0 ){
        #endif
          if(static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(var.data()+j)){
          auto indexz = static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()
              )->get_function())->is_active(var.data()+j);
          auto fobj_sbi = static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()
              )->get_function());
          if( ! fobj_sbi->is_linear() ){
            addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
            addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
          } else {
            addval1 = static_cast< p_LF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_coefficient(indexz);
          }
          idx_to_var_sbi1[index].push_back(double_var( addval1 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( addval2 , var.data()+j ));
        } else {
          idx_to_var_sbi1[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( 0.0 , var.data()+j ));
			  }
        pos_id++;
        NumStatVar++;
        #ifdef BIN_VARS 
        } 
        #endif
        } } ) )
    continue;
   // Multiarrays
   if( un_any_thing_K( ColVariable , el ,
		  {
      for (Index j = 0 ; j < var.num_elements() ; ++j){
        #ifdef BIN_VARS 
			  if( (var.data()+j)->is_integer() && (var.data()+j)->get_lb() == 0.0 && (var.data()+j)->get_ub() == 1.0 ){
        #endif
          if(static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(var.data()+j)){ 
          auto indexz = static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective()
          )->get_function())->is_active(var.data()+j);
          if( ! static_cast< Function * >( static_cast< p_FRO >( sbi->get_objective())->get_function())->is_linear() ){
            addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
            addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
          } else {
            addval1 = static_cast< p_LF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_coefficient(indexz);
          }
          idx_to_var_sbi1[index].push_back(double_var( addval1 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( addval2 , var.data()+j ));
        } else {
          idx_to_var_sbi1[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( 0.0 , var.data()+j ));
			  }
			  pos_id++;
        NumStatVar++;
			  #ifdef BIN_VARS 
        } 
        #endif
        } } ) )
     continue;
   //throw( std::invalid_argument(
   //             "PrimalProximalHeur: static variable not a static binary ColVariable" ) );
   }
    pos_id_sbi[index] = pos_id;
    index++;
 }

 if( NumStatVar == 0 )
  throw( std::invalid_argument(
		      "PrimalProximalHeur::setBlock: no static binary variable" ) );

 }  // end( PrimalProximalHeur::initialize )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIterPPH ): maxIter = value; break;
  case( intMaxSol ):  f_MaxSol = value; break;
  case( intLogVerb ):
   logVerb = value & 3;
   LagrangianDualSolver::set_par( par , std::max( 0 , value >> 2 ) );
   break;
  default: LagrangianDualSolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dbl_penaltyFactor ): R = value; break;
  default:                   LagrangianDualSolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int PrimalProximalHeur::compute( bool changedvars )
{
 Index iters = 0;
 bool is_the_same = false;
 int res;

 best_bound = f_max ? -Inf<double>() : Inf<double>();

 lock();  // lock the mutex

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "\nNumStatVar: " << NumStatVar << "\n";

 std::vector< double > sol( NumStatVar );

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "PrimalProximalHeur: unable to lock the Block" ) );
                  
 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );

 while( ( ! is_the_same ) && ( iters < maxIter ) ) {
  // main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "\niteration = " << iters << "\n";
/*
  if( iters == 0 )
   for( Index kvar = 0 ; kvar < NumStatVar ; ++kvar )
    sol[ kvar ] = rand() % 2;
*/
  if( iters >= 1 )
    previous_sol.insert( previous_sol.begin() , & sol[ 0 ] ,
		       & sol[ NumStatVar ] );

  if( iters >= 1 ) {
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "\nADDING PENALTY TERMS\n\n";

   add_penalty_terms();
   //changed_penalties = true;
   }

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "COMPUTE SOLUTION" << std::endl;

  if( iters >= 1 ) { 
   LagrangianDualSolver::set_event_handler(
     ThinComputeInterface::eEverykIteration ,
     [ this ] () { 
        if( f_Block->is_feasible() ) {
          Index index = 0;
          for( const auto & sbi : f_Block->get_nested_Blocks() ) {
            for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ){
              auto si = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
              if( si * ( 1 - si ) > 1e-3 ){
                *f_log << "IS_NOT_INTEGER_SOL" << std::endl;
                return( ThinComputeInterface::eContinue );
              }
            }
            index++;
          }
          value_FUNCTION = get_funct_value() ;
          add_penalty_terms();
          // if new solution is feasible, add to v_best_sol 
          // and possibly update the best bound 
          if( f_log && ( logVerb >= 2 ) )
            *f_log << "IS_FEASIBLE_SOL: " << value_FUNCTION << std::endl;

          // better than the best
          bool better = f_max ? ( value_FUNCTION > best_bound ) :
                                ( value_FUNCTION < best_bound );
          if( better )
            best_bound = value_FUNCTION;

          // worse than the worst
          bool worse = f_max ? ( value_FUNCTION < worst_bound ) :
                                ( value_FUNCTION > worst_bound );

          if( v_best_sol.size() < f_MaxSol ) {
            // there is free space, just throw the new solution in
            v_best_sol.push_back( std::pair( f_Block->get_Solution() ,
                  value_FUNCTION ) );
            std:: push_heap( v_best_sol.begin() , v_best_sol.end() );
            if( worse )
            worst_bound = value_FUNCTION;
            }
          else
            if( ! worse ) {
            // there is no space, so throw it in only if it's better than the worst
            // find the position of the element with the worst value, and find the
            // second-worst value to update worst_bound;
            double second_worst = f_max ? Inf< double >() : - Inf< double >();
            std::vector< sol_value >::iterator bad;
            for( auto it = v_best_sol.begin() ; it != v_best_sol.end() ;
            ++it )
              if( it->second == worst_bound )
              bad = it;
              else {
                bool s_worse = f_max ? ( it->second < second_worst ) :
                                        ( it->second > second_worst );
                if( s_worse )
                  second_worst = it->second;
              }
            // replace the worst element with the new one
            *bad = std::pair( f_Block->get_Solution() , value_FUNCTION );
            // re-make the heap, since we may have just invalidated it
            std::make_heap( v_best_sol.begin() , v_best_sol.end() ,
                [ this ]( const sol_value & a , const sol_value & b ) {
                  return( f_max ? ( a.second < b.second ) :
                            ( a.second > b.second ) );
                  } );
            // the worst solution is now either the second worst of the new one
            if( f_max )
              worst_bound = std::min( second_worst , value_FUNCTION );
            else
              worst_bound = std::max( second_worst , value_FUNCTION );
            }
          }
          else
          if( f_log && ( logVerb >= 2 ) )
            *f_log << "IS_INFEASIBLE_SOL: " << value_FUNCTION << std::endl;
        return( ThinComputeInterface::eContinue );
      }  // end of lambda
		);
  }

  res = InnerSolver->compute( changedvars );
  if( f_log && ( logVerb >= 2 ) )
   *f_log << "SOLUTION COMPUTED" << std::endl;

  if( iters >= 0 ) {
   Index kvar = 0;
   Index index = 0;

   for( const auto & sbi : f_Block->get_nested_Blocks() ) {
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ){
     sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
     //std::cout << "sol[ " << kvar << " ] = " << sol[ kvar ] << std::endl;
     kvar++;
    }
    index++;
    }
   }

/*
  if( f_log && ( logVerb >= 2 ) ){
   auto solu = f_Block->get_Solution();
   solu->print( *f_log );
  }
*/

  if( iters >= 0 ) {
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "\nREMOVING PENALTY TERMS\n\n";

   remove_penalty_terms();
   //!changed_penalties = true;
   //!process_outstanding_Modification();
   }

  penalty = 0.0;
  addterm = 0.0;

  if(iters >= 1 ){
    for( int ivar = 0 ; ivar < NumStatVar ; ++ivar ) {
      auto si = sol[ ivar ];
      auto psi = previous_sol[ ivar ];
      penalty += R * ( psi - si ) * ( psi - si );
      addterm += R * si * ( 1.0 - 2.0 * psi );
    }
  }

   auto value_bound = f_max ? InnerSolver->get_ub() - addterm :
    InnerSolver->get_lb() - addterm;

   value_FUNCTION = get_funct_value();
 
   if( f_log && ( logVerb >= 2 ) )    
    if( std::abs( value_bound - value_FUNCTION ) / 
      std::max( std::abs( value_bound ) , std::abs( value_FUNCTION ) ) >= 1e-3 )
      *f_log << "ERROR: " << value_bound - value_FUNCTION << std::endl;

  if( f_log && ( logVerb >= 2 ) ) {
   *f_log << "ITERS: " << iters << std::endl;
   *f_log << "InnerSolver LB: " << InnerSolver->get_lb() << std::endl;
   *f_log << "InnerSolver UB: " << InnerSolver->get_ub() << std::endl;
   *f_log << "penalty: " << penalty << std::endl;
   *f_log << "addterm: " << addterm << std::endl;
   *f_log << "SOL1: " << value_bound << std::endl;
   *f_log << "SOL2: " << value_FUNCTION << std::endl;
   auto bound = f_max ? get_lb() : get_ub();
   *f_log << "Best Value: " << bound << std::endl;
  } 

  if( iters >= 1 ) {  // stop criterion : check if the current and
                      // previous solutions are the same
   // TODO: put some parameter instead 1e-6
  
   /*
   is_the_same = true;
   for( int ivar = 0 ; ivar < NumStatVar ; ++ivar )
    if( std::abs( sol[ ivar ] - previous_sol[ ivar ] ) > 1e-6 ) { 
     is_the_same = false;
     break;
     }
  */

  is_the_same = penalty < 1e-3 ? true : false;
  }

  Index index = 0;
  bool is_integer = true;
  for( const auto & sbi : f_Block->get_nested_Blocks() ) {
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ){
      auto si = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
      if( si * ( 1 - si ) > 1e-3 ){
        is_integer = false;
        *f_log << "IS_NOT_INTEGER_SOL: " << si << std::endl;
        break;
      }
      if( ! is_integer )
        break;
    }
    index++;
   }

   /// if( f_Block->is_feasible() && iters >= 1 ) {
   if( f_Block->is_feasible() && is_integer ) {
   // if new solution is feasible, add to v_best_sol 
   // and possibly update the best bound 
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "IS_FEASIBLE_SOL" << std::endl;

   // better than the best
   bool better = f_max ? ( value_FUNCTION > best_bound ) :
                         ( value_FUNCTION < best_bound );
   if( better )
    best_bound = value_FUNCTION;

   // worse than the worst
   bool worse = f_max ? ( value_FUNCTION < worst_bound ) :
                        ( value_FUNCTION > worst_bound );

   if( v_best_sol.size() < f_MaxSol ) {
    // there is free space, just throw the new solution in
    v_best_sol.push_back( std::pair( f_Block->get_Solution() ,
					 value_FUNCTION ) );
    std:: push_heap( v_best_sol.begin() , v_best_sol.end() );
    if( worse )
     worst_bound = value_FUNCTION;
    }
   else
    if( ! worse ) {
     // there is no space, so throw it in only if it's better than the worst
     // find the position of the element with the worst value, and find the
     // second-worst value to update worst_bound;
     double second_worst = f_max ? Inf< double >() : - Inf< double >();
     std::vector< sol_value >::iterator bad;
     for( auto it = v_best_sol.begin() ; it != v_best_sol.end() ;
	  ++it )
      if( it->second == worst_bound )
       bad = it;
      else {
       bool s_worse = f_max ? ( it->second < second_worst ) :
                              ( it->second > second_worst );
       if( s_worse )
	      second_worst = it->second;
       }
     // replace the worst element with the new one
     *bad = std::pair( f_Block->get_Solution() , value_FUNCTION );
     // re-make the heap, since we may have just invalidated it
     std::make_heap( v_best_sol.begin() , v_best_sol.end() ,
		     [ this ]( const sol_value & a , const sol_value & b ) {
		      return( f_max ? ( a.second < b.second ) :
			              ( a.second > b.second ) );
		      } );
     // the worst solution is now either the second worst of the new one
     if( f_max )
      worst_bound = std::min( second_worst , value_FUNCTION );
     else
      worst_bound = std::max( second_worst , value_FUNCTION );
     }
   }
  else
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "IS_INFEASIBLE_SOL" << std::endl;

  iters++;

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if iBCopy == false, bring back the inner Block to its original objective
 // "like if nothing ever happened",
 if( ! iBCopy ) {
  for( auto lbf : v_LBF )
   lbf->cleanup_inner_objective();
  }

 // because the inner Solver is solving the dual of the original Block,
 // the unbounded an unfeasible return states have to be exchanged
 if( res == kUnbounded ) {
  res = kInfeasible;
  if( f_log && ( logVerb >= 2 ) ) {
   *f_log << "\nINFEASIBLE!\n";
   *f_log << "NUMBER ITERS: " << iters-1 << "\n";
   *f_log << "EXIT CODE: " << res << "\n\n";
   }
  return( res );
  }
  else {
  if( res == kInfeasible ){
   res = kUnbounded;
   if( f_log && ( logVerb >= 2 ) ){
    *f_log << "\nUNBOUNDED!\n";
    *f_log << "NUMBER ITERS: " << iters-1 << "\n";
    *f_log << "EXIT CODE: " << res << "\n\n";
    }
   return( res );
   }
  }

 if( is_the_same )
  if( f_log && ( logVerb >= 2 ) ) {
   *f_log << "R = " << R << std::endl;
   *f_log << "maxIterPPH = " << maxIter << std::endl;
   *f_log << "IS_THE_SAME = TRUE\n"; 
   *f_log << "NUMBER ITERS: " << iters-1 << "\n";
   *f_log << "LB: " << InnerSolver->get_lb() << "\n";
   *f_log << "UB: " << InnerSolver->get_ub() << "\n";
   auto bound = f_max ? get_lb() : get_ub();
   *f_log << "Best Feasible solution: " << bound << std::endl;
   }

 changed_penalties = false;

 unlock();  // unlock the mutex


 return( res );

 }  // end( PrimalProximalHeur::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::process_outstanding_Modification( void )
{

 //!LagrangianDualSolver::process_outstanding_Modification();    

 bool reload = false;
 bool check_feasibility = false;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod ){
   for( const auto & sbi : f_Block->get_nested_Blocks() ) {
    if( mod->get_Block() == sbi ) {
     if( auto tmod = dynamic_cast< C05FunctionModLin * >( mod.get() ) ) {
      // modifications in the objective of sub-blocks: possible reloading
      // of the static dictionaries for objective (linear) terms
      reload = true;
      //! changed_penalties = false;
     } else {
      // modifications in the constraints of the objective of the sub-block:
      // possibly some of the (best) solutions become infeasible
      // (check feasibility) 
      check_feasibility = true;
      reload = true;
     }
    }
   }
  }

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "check_feasibility = " << check_feasibility << std::endl;

 if( reload and !changed_penalties ) {
  if( f_log && ( logVerb >= 2 ) )
    *f_log << "reload..." << std::endl;
  guts_of_destructor();
  initialize();
  }

 if( check_feasibility ) {
  // check feasibility of the saved solutions, remove those that have
  // become unfeasible
  std::vector< sol_value > v_best_sol_new;

  for( auto & sol : v_best_sol ) {
   //std::cout << v_best_sol.size() << std::endl;
   sol.first->write( f_Block );
   if( f_Block->is_feasible() )
    v_best_sol_new.push_back( sol );
   else
    delete sol.first;
   }

  v_best_sol = v_best_sol_new;
  std::make_heap( v_best_sol.begin() , v_best_sol.end() ,
		  [ this ]( const sol_value & a , const sol_value & b ) {
		   return( f_max ? ( a.second < b.second ) :
			           ( a.second > b.second ) );
		   } );
  }

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( v_mod.empty() )  // no Modification coming directly from f_Block
  return;             // all done

 }  // end( PrimalProximalHeur::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::add_penalty_terms( void )
{
 // adding the penality terms to the (quadratic) objective functions of the 
 // sub-blocks using the static dictionaries 

 Index pos = 0;
 Index index = 0;
 index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobji = static_cast< Function * >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );

  if( ! fobji->is_linear() ) {       
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    auto indexz = fobji->is_active( idx_to_var_sbi1[ index ][ ivar ].second );
    auto addval1 = idx_to_var_sbi1[ index ][ ivar ].first;
    auto addval2 = idx_to_var_sbi2[ index ][ ivar ].first;
    auto fobji1 = static_cast< p_DQF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );
    if( fobji1->get_num_active_var() > indexz ){
      #ifdef BIN_VARS 
        fobji1->modify_linear_coefficient( indexz , addval1 + R * 
            ( 1.0 - 2.0 * previous_sol[ pos ] ) , mp );
      #else
        fobji1->modify_term( indexz , addval1 - R * 2.0 * previous_sol[ pos ] , 
            addval2 + R , mp );
      #endif
    } else
      if( R > 0.0 )
        #ifdef BIN_VARS 
          fobji1->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
            R * ( 1.0 - 2.0 * previous_sol[ pos ] ) , 0.0 , mp );   
        #else
          fobji1->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
            - R * 2.0 * previous_sol[ pos ] , R , mp );   
        #endif
    pos++;
    }
  } else {
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    auto indexz = fobji->is_active( idx_to_var_sbi1[ index ][ ivar ].second );
    auto addval1 = idx_to_var_sbi1[ index ][ ivar ].first;
    auto fobji1 = static_cast< p_LF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );
    if( fobji1->get_num_active_var() > indexz )
      fobji1->modify_coefficient( indexz , addval1 + R * ( 1.0 - 2.0 * previous_sol[ pos ] ) , mp );
    else
      if( R > 0.0 )
        fobji1->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
            R * ( 1.0 - 2.0 * previous_sol[ pos ] ) , mp );       
    pos++;
    }
  }

  sbi->close_channel( chnl );
  index++;
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::remove_penalty_terms( void )
{
 // removing the penality terms to the (quadratic) objective functions of the 
 // sub-blocks using the static dictionaries 

 Index pos = 0;
 Index index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobji = static_cast< Function * >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );

  if( ! fobji->is_linear() ) {       
    for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    auto indexz = fobji->is_active( idx_to_var_sbi1[ index ][ ivar ].second );
    auto addval1 = idx_to_var_sbi1[ index ][ ivar ].first;
    auto addval2 = idx_to_var_sbi2[ index ][ ivar ].first;
    auto fobji1 = static_cast< p_DQF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );
    if( fobji1->get_num_active_var() > indexz ){
        #ifdef BIN_VARS 
        fobji1->modify_linear_coefficient( indexz , addval1 , mp );
        #else
        fobji1->modify_term( indexz , addval1 , addval2 , mp ); 
        #endif
    }
    pos++;
    }
  } else {
   for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
    auto indexz = fobji->is_active( idx_to_var_sbi1[ index ][ ivar ].second );
    auto addval1 = idx_to_var_sbi1[ index ][ ivar ].first;
    auto fobji1 = static_cast< p_LF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );
    if( fobji1->get_num_active_var() > indexz )
      fobji1->modify_coefficient( indexz , addval1 , mp );
    }
  }

  sbi->close_channel(chnl);
  index += 1;
  }
}

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::guts_of_destructor( void )
{
 for( auto & el : v_best_sol )
  delete el.first;
 v_best_sol.clear();
 idx_to_var_sbi1.clear();
 Funct_sbi.clear();
 previous_sol.clear();
 pos_id_sbi.clear();

 }  // end( PrimalProximalHeur::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*------------------- End File PrimalProximalHeur.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
