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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PrimalProximalHeur.h"

#include "FRealObjective.h"

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
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::initialize( void )
{
 // count and check the static ColVariable - - - - - - - - - - - - - - - - - 
 // meanwhile construct the static dictionaries for the linear and quadratic
 // terms of the (quadratic) objective functions of the sub-blocks

 NumStatVar = 0;

{
 pos_id_sbi.resize( f_Block->get_number_nested_Blocks() );
 idx_to_var_sbi1.resize( f_Block->get_number_nested_Blocks() );
 idx_to_var_sbi2.resize( f_Block->get_number_nested_Blocks() );
 Funct_sbi.resize( f_Block->get_number_nested_Blocks()) ;
 
 Index t = -1;
 Index index = 0;
 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  Funct_sbi[index] = static_cast< p_DQF >(static_cast< p_FRO >( sbi->get_objective())->get_function());
  pos_id = 0;
  for( const auto & el : sbi->get_static_variables() ) {
   // Singles
    t++;
   if( un_any_thing_0( ColVariable , el ,
		  {
      if(static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(&var)){ 
        auto indexz = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
        )->get_function())->is_active(&var);
        auto addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
        auto addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
			  idx_to_var_sbi1[index].push_back(double_var( addval1 , &var ));
        idx_to_var_sbi2[index].push_back(double_var( addval2 , &var ));
        idx_to_var1.push_back(double_var( addval1 , & var ));
        idx_to_var2.push_back(double_var( addval2 , & var ));
      } else {
        idx_to_var_sbi1[index].push_back(double_var( 0.0 , &var ));
        idx_to_var_sbi2[index].push_back(double_var( 0.0 , &var ));
        idx_to_var1.push_back(double_var( 0.0 , & var ));
        idx_to_var2.push_back(double_var( 0.0 , & var ));
      }
      pos_id++;
			NumStatVar++;
          } ) )
     continue;
   // Vectors
   if( un_any_thing_1( ColVariable , el ,
		  {
			for (Index j = 0 ; j < var.size() ; ++j){
        if(static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(var.data()+j)){ 
          auto indexz = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
          )->get_function())->is_active(var.data()+j);
          auto addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
          auto addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
          idx_to_var_sbi1[index].push_back(double_var( addval1 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( addval2 , var.data()+j ));
          idx_to_var1.push_back(double_var( addval1 , var.data()+j ));
          idx_to_var2.push_back(double_var( addval2 , var.data()+j ));
        } else {
          idx_to_var_sbi1[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var1.push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var2.push_back(double_var( 0.0 , var.data()+j ));
			  }
        pos_id++;
      }
			NumStatVar += var.size();
          } ) )
    continue;
   // Multiarrays
   if( un_any_thing_K( ColVariable , el ,
		  {
			NumStatVar += var.num_elements();
      for (Index j = 0 ; j < var.num_elements() ; ++j){
			  if(static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()                                                                                                         
            )->get_function())->is_active(var.data()+j)){ 
          auto indexz = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
          )->get_function())->is_active(var.data()+j);
          auto addval1 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_linear_coefficient(indexz);
          auto addval2 = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_quadratic_coefficient(indexz);
          idx_to_var_sbi1[index].push_back(double_var( addval1 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( addval2 , var.data()+j ));
          idx_to_var1.push_back(double_var( addval1 , var.data()+j ));
          idx_to_var2.push_back(double_var( addval2 , var.data()+j ));
        } else {
          idx_to_var_sbi1[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var_sbi2[index].push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var1.push_back(double_var( 0.0 , var.data()+j ));
          idx_to_var2.push_back(double_var( 0.0 , var.data()+j ));
			  }
			  pos_id++;
			}
          } ) )
     continue;
   throw( std::invalid_argument(
                "PrimalProximalHeur: static variable not a ColVariable" ) );
   }
     pos_id_sbi[index] = pos_id;
    index += 1;
  }
 }

 if( idx_to_var_sbi1.empty() && idx_to_var_sbi2.empty())
  throw( std::invalid_argument(
		      "PrimalProximalHeur::setBlock: no static variable" ) );

 std::sort( var_to_idx.begin() , var_to_idx.end() );



 }  // end( PrimalProximalHeur::initialize )

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): maxIter = value; break;
  case( intMaxSol ):  f_MaxSol = value; break;
  case( intLogVerb ):
   logVerb = value & 3;
   LagrangianDualSolver::set_par( par , std::max( 0 , value >> 2 ) );
   break;
  case( intMaxIterLD ): LagrangianDualSolver::set_par( intMaxIter , value );
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
 double integer_viol;
 double integer_viol_sum;
 double funct_old = 0.0;
 value_FUNCTION = Inf< double >();

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "\nNumStatVar: " << NumStatVar << "\n";

 std::vector< double > sol( NumStatVar );

 while( ( ! is_the_same ) && ( iters < maxIter ) ) {
  // main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "\niteration = " << iters << "\n";

  changed_penalties = false;

  bool owned = f_Block->is_owned_by( f_id );
  if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
   throw( std::runtime_error(
                       "PrimalProximalHeur: unable to lock the Block" ) );

  process_outstanding_Modification();

  if( ! owned )
   f_Block->unlock( f_id );

  if( iters == 0 )
   for( Index kvar = 0 ; kvar < NumStatVar ; ++kvar )
    sol[ kvar ] = rand() % 2;

  previous_sol.insert( previous_sol.begin() , & sol[ 0 ] ,
		       & sol[ NumStatVar ] );

  if( iters >= 1 ) {
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "\nADDING PENALTY TERMS\n\n";

   add_penalty_terms();
   changed_penalties = true;
   process_outstanding_Modification();
   }

  if( f_log && ( logVerb >= 2 ) )
   *f_log << "COMPUTE SOLUTION" << std::endl;

  res = InnerSolver->compute( changedvars );
  if( f_log && ( logVerb >= 2 ) )
   *f_log << "SOLUTION COMPUTED" << std::endl;

  if( iters >= 0 ) {
   Index kvar = 0;
   double auxl = 0.0;
   for( Index kvar = 0 ; kvar < NumStatVar ; ++kvar )
    sol[ kvar ] = idx_to_var1[ kvar ].second->get_value();
   }

  if( iters >= 1 ) {
   if( f_log && ( logVerb >= 2 ) )
    *f_log << "\nREMOVING PENALTY TERMS\n\n";

   remove_penalty_terms();
   changed_penalties = true;
   process_outstanding_Modification();
   }

  if( f_log && ( logVerb >= 2 ) ) {
   *f_log << "ITERS: " << iters << std::endl;
   *f_log << "SOL LB: " << InnerSolver->get_lb() << std::endl;
   *f_log << "SOL UB: " << InnerSolver->get_ub() << std::endl;
   *f_log << "FUNCT Value: " << get_funct_value() << std::endl;
   auto bound = f_max ? get_lb() : get_ub();
   *f_log << "Best Value: " << bound << std::endl;
   }

  penalty = 0.0;
  double addterm = 0.0;
  for( int ivar = 0 ; ivar < NumStatVar ; ++ivar ) {
   auto si = sol[ ivar ];
   auto psi = previous_sol[ ivar ];
   penalty += R * ( psi - si ) * ( psi - si );
   addterm += R * si * ( si - 2.0 * psi );
   }

  if( f_log && ( logVerb >= 2 ) ) {
   *f_log << "penalty: " << addterm << std::endl;
   *f_log << "SOL: " << InnerSolver->get_lb() - addterm << std::endl;
   }

  if( f_Block->is_feasible() ) {
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
    
  if( iters >= 0 ) {  // stop criterion : check if the current and
                      // previous solutions are the same
   // TODO: put some parameter instead 1e-6
   is_the_same = true;
   for( int ivar = 0 ; ivar < NumStatVar ; ++ivar )
    if( std::abs( sol[ ivar ] - previous_sol[ ivar ] ) > 1e-6 ) { 
     is_the_same = false;
     break;
     }
   }

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
   *f_log << "IS_THE_SAME = TRUE\n"; 
   *f_log << "NUMBER ITERS: " << iters-1 << "\n";
   *f_log << "LB: " << InnerSolver->get_lb() << "\n";
   *f_log << "UB: " << InnerSolver->get_ub() << "\n";
   auto bound = f_max ? get_lb() : get_ub();
   *f_log << "Best Feasible solution: " << bound << std::endl;
   }

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
 bool reload = false;
 bool check_feasibility = false;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod ){
  if( mod->get_Block() == f_Block ) {
    //LagrangianDualSolver::process_outstanding_Modification();
   }
  else { 
   for( const auto & sbi : f_Block->get_nested_Blocks() ) {
    if( mod->get_Block() == sbi && ! changed_penalties ) {
     if( auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod.get() ) ) 
      // modifications in the objective of sub-blocks: possible reloading
      // of the static dictionaries for objective (linear and quadratic) terms
      reload = true;
     else
      // modifications in the constraints of the objective of the sub-block:
      // possibly some of the (best) solutions become infeasible
      // (check feasibility) 
      check_feasibility = true;
      }
    }
   }
  }

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "check_feasibility = " << check_feasibility << std::endl;

 if( reload ) {
  if( f_log && ( logVerb >= 2 ) )
    *f_log << "reload..." << std::endl;
  initialize();
  }

 if( check_feasibility ) {
  // check feasibility of the saved solutions, remove those that have
  // become unfeasible
  std::vector< sol_value > v_best_sol_new;
  for( auto & sol : v_best_sol ) {
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
 penalty = 0.0;
 index = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobji = static_cast< p_DQF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );

  for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar ) {
   auto indexz = fobji->is_active( idx_to_var_sbi1[ index ][ ivar ].second );
   if( fobji->get_num_active_var() > indexz ) {
    auto addval1 = idx_to_var_sbi1[ index ][ ivar ].first;
    auto addval2 = idx_to_var_sbi2[ index ][ ivar ].first;  
    fobji->modify_term( indexz , addval1 - R * 2.0 * previous_sol[ pos ] ,
			addval2 + R , mp );
    }
   else
    if( R > 0.0 )
     fobji->add_variable( idx_to_var_sbi1[ index ][ ivar ].second ,
			  - R * 2.0 * previous_sol[ pos ] , R , mp );           
   pos++;
   }

  sbi->close_channel( chnl );
  index += 1;
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::remove_penalty_terms( void )
{
 // removing the penality terms to the (quadratic) objective functions of the 
 // sub-blocks using the static dictionaries 

 Index pos = 0;
 Index index = 0;
 double value_FUNCTION1 = 0.0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
  auto fobji = static_cast< p_DQF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );

  for( Index ivar = 0 ; ivar < pos_id_sbi[index] ; ++ivar ) {
   if( fobji->get_num_active_var() >
       fobji->is_active(idx_to_var_sbi1[index][ivar].second)) {
    auto indexz = fobji->is_active(idx_to_var_sbi1[index][ivar].second);
    auto addval1 = idx_to_var_sbi1[index][ivar].first; 
    auto addval2 = idx_to_var_sbi2[index][ivar].first; 
    fobji->modify_term( indexz , addval1 , addval2, mp );
    }
   else {
    fobji->add_variable( idx_to_var_sbi1[index][ivar].second ,
			 idx_to_var_sbi1[index][ivar].first ,
			 idx_to_var_sbi2[index][ivar].first , mp );
    }
   pos++;
   }
  
  // compute the corresponding the objective function value for the current 
  // solution (use the original objective functions without penalty terms)

  Funct_sbi[index]->compute(true);
  value_FUNCTION1 += Funct_sbi[index]->get_value();  
  sbi->close_channel(chnl);
  index += 1;
  }

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "value_FUNCTION1 = " << value_FUNCTION1 << std::endl;

 value_FUNCTION = 0.0;
 index = 0;

 // compute the corresponding the objective function value for the current 
 // solution (use the original objective functions without penalty terms)
 // double check: value_FUNCTION == value_FUNCTION1

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  auto fobji = static_cast< p_DQF >( static_cast< p_FRO >(
				  sbi->get_objective() )->get_function() );
  fobji->compute( true );
  value_FUNCTION += fobji->get_value();
  }

 if( f_log && ( logVerb >= 2 ) )
  *f_log << "value_FUNCTION = " << value_FUNCTION << std::endl;
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::guts_of_destructor( void )
{
 for( auto & el : v_best_sol )
  delete el.first;
 v_best_sol.clear();
 idx_to_var_sbi2.clear();
 idx_to_var_sbi1.clear();
 Funct_sbi.clear();
 idx_to_var2.clear();
 idx_to_var1.clear();
 var_to_idx.clear();
 previous_sol.clear();
 pos_id_sbi.clear();

 }  // end( PrimalProximalHeur::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*------------------- End File PrimalProximalHeur.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
