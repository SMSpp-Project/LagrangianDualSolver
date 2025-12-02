/*--------------------------------------------------------------------------*/
/*--------------------- File PrimalProximalHeur.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework for a "generic"
 * Lagrangian-based Solver.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BlockSolverConfig.h"

#include "PrimalProximalHeur.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "RBlockConfig.h"

#include "ColVariableSolution.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Bitwise-coded macro that activates costly checks that should never be done
  * in production, but can be useful during debugging. Currently supported
  * checks are:
  *
  * - bit 0 (+ 1): is_correct() is called on the Lagrangian Dual Block to
  *   verify that all Variable and Constraint are properly linked. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using VarValue = Function::FunctionValue;
using c_VarValue = Function::c_FunctionValue;
using Vec_FunctionValue = Function::Vec_FunctionValue;

using Vec_VarValue = Function::Vec_FunctionValue;
using c_Vec_VarValue = Function::c_Vec_FunctionValue;

using coeff_pair = LinearFunction::coeff_pair;
using v_coeff_pair = LinearFunction::v_coeff_pair;
using v_c_coeff_pair = LinearFunction::v_c_coeff_pair;

using LinearCombination = C05Function::LinearCombination;
using c_LinearCombination = C05Function::c_LinearCombination;

using dual_pair = LagBFunction::dual_pair;
using v_dual_pair = std::vector< dual_pair >;
using v_c_dual_pair = const v_dual_pair;

using p_AB = AbstractBlock *;
using p_FRC = FRowConstraint *;
using p_FRO = FRealObjective *;
using p_LF = LinearFunction *;
using p_LBF = LagBFunction *;
using p_DQF = DQuadFunction *;

using p_BC = BlockConfig *;
using p_BSC = BlockSolverConfig *;

using p_SConf_p_p = SimpleConfiguration< std::pair< Configuration * ,
						    Configuration * > > *;

/*--------------------------------------------------------------------------*/
/* PrimalProximalHeur always need to have an "inner CDASolver" set, but at
 * the beginning it would have none: it will then have a FakeCDASolver one,
 * which does nothing. Clearly, a "real" one will have to be set of
 * PrimalProximalHeur::compute() is to work. */
/*
class FakeCDASolver : public CDASolver {
 public:
 FakeCDASolver( void ) : CDASolver() {}
 virtual ~FakeCDASolver() {}
 int compute( bool changedvars = true ) override final { return( kError ); }
 bool has_var_solution( void ) override final { return( false ); }
 void get_var_solution( Configuration *solc = nullptr ) override final {}
 bool has_dual_solution( void ) override final { return( false ); }
 void get_dual_solution( Configuration *solc = nullptr ) override final {}
 void add_Modification( sp_Mod & mod ) override final {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_0( FakeCDASolver );
*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr VarValue NaNshift
                              = std::numeric_limits< VarValue >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

static constexpr VarValue INFshift = Inf< VarValue >();
 ///< convenience constexpr for "Infty"

static constexpr Index InINF = SMSpp_di_unipi_it::Inf< Index >();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register PrimalProximalHeur to the Solver factory

SMSpp_insert_in_factory_cpp_0( PrimalProximalHeur );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF PrimalProximalHeur ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::initialize(){

  // count and check the binary ColVariable - - - - - - - - - - - - - - - - - 
 // meanwhile construct the static dictionaries

 NumBinStatVar = 0;

 int svn = 0;

 for( const auto & sbi : f_Block->get_nested_Blocks() )
  svn += sbi->get_static_variables().size();

 //std::cout << "svn: " << svn << "\n";
 //var_to_idx.resize( svn );
 //idx_to_var.resize( svn );
{
 pos_id_sbi.resize(f_Block->get_number_nested_Blocks());
 idx_to_var_sbi1.resize(f_Block->get_number_nested_Blocks());
 idx_to_var_sbi2.resize(f_Block->get_number_nested_Blocks());
 Funct_sbi.resize(f_Block->get_number_nested_Blocks());
 
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
      ints.push_back(t);
			int_vars.push_back(&var);
      pos_id++;
			NumBinStatVar++;
          } ) )
     continue;
   // Vectors
   if( un_any_thing_1( ColVariable , el ,
		  {
			for (Index j = 0 ; j < var.size() ; ++j){
			  int_vars.push_back(var.data()+j);
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
			NumBinStatVar += var.size();
			ints.push_back(t);
          } ) )
    continue;
   // Multiarrays
   if( un_any_thing_K( ColVariable , el ,
		  {
			NumBinStatVar += var.num_elements();
			ints.push_back(t);
      for (Index j = 0 ; j < var.num_elements() ; ++j){
			  int_vars.push_back(var.data()+j);
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

}

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , double value )
{
  switch( par ) {
    case( dbl_penaltyFactor ): 
      R = value;
      break;
  default:
    LagrangianDualSolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intLogVerb ):
    InnerSolver->set_par( int_par_lds( par ) , std::max( 0 , value - 2 ));
    logVerb = value;
  break;
  default:
   LagrangianDualSolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int PrimalProximalHeur::compute( bool changedvars )
{

 initialize();

 Index iters = 0;
 bool is_the_same = false;
 bool is_integer_sol = false; 
 Solution *solution;
 int res;
 double integer_viol;
 double integer_viol_sum;
 double funct_old = 0.0;
 value_FUNCTION = Inf<double>();

 if( logVerb - get_int_par( intLogVerb ) >= 2 ){
  std::cout << "\nNumBinStatVar: " << NumBinStatVar << "\n";
 }

 double sol[NumBinStatVar];

 while ( !is_the_same and iters < maxIter ){

 if( logVerb - get_int_par( intLogVerb ) >= 2 ){
   std::cout << "\niteration = " << iters << "\n";
 }

 changed_penalties = false;

 lock();  // lock the mutex

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "PrimalProximalHeur: unable to lock the Block" ) );

 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );

 /* This is no longer needed, since these Modification happen when
    f_play_dumb == true in the inner LagBFunction

 // if iBCopy == false, inhibit all Modification from the UpdateSolver; these
 // would reach the (disconnected) original Block, but there is no reason for
 // this because these are all "local" changes that will be undone at the end
 // of the solution process, so the Solver attached to the father (and other
 // ancestors) have no reason to act upon them
 if( ! iBCopy )
  for( auto us : v_US )
   us->inhibit_Modification( true );
*/

  //int ivar = 0;

 if(iters==0){
  for( Index kvar = 0 ; kvar < NumBinStatVar ; ++kvar ){
    sol[kvar] = rand() % 2;
  }
 }
 
 previous_sol.insert(previous_sol.begin(), &sol[0], &sol[NumBinStatVar]);

  if( iters >= 1 ){
    if( logVerb - get_int_par( intLogVerb ) >= 2 ){
      std::cout << "\nADDING PENALTY TERMS\n\n";
    }
    add_penalty_terms();
    changed_penalties = true;
    process_outstanding_Modification();
  }

  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "COMPUTE SOLUTION" << std::endl;
  }
  auto res = InnerSolver->compute( changedvars );
  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "SOLUTION COMPUTED" << std::endl;
  }

if ( iters >= 0 ){
    InnerSolver->get_dual_solution();
    solution = f_Block->get_Solution();

    auto sols = new ColVariableSolution;
    sols->read( f_Block );
    sols->write( f_Block );

    Index kvar = 0;
    double auxl = 0.0;
    for( Index kvar = 0 ; kvar < NumBinStatVar ; ++kvar ){
      sol[kvar] = idx_to_var1[kvar].second->get_value();
    }
  }

  if( iters >= 1 ){
    if( logVerb - get_int_par( intLogVerb ) >= 2 ){
      std::cout << "\nREMOVING PENALTY TERMS\n\n";
    }
    remove_penalty_terms();
    changed_penalties = true;
    process_outstanding_Modification();
  }

  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "ITERS: " << iters << std::endl;
    std::cout << "SOL LB: " << InnerSolver->get_lb() << std::endl;
    std::cout << "SOL UB: " << InnerSolver->get_ub() << std::endl;
    std::cout << "FUNCT Value: " << get_funct_value() << std::endl;
    std::cout << "Best Value: " << get_best_bound() << std::endl;
  }

  penalty = 0.0;
  double addterm = 0.0;
  for( int ivar = 0 ; ivar < NumBinStatVar ; ++ivar ){
    penalty += R* pow((previous_sol[ivar] - sol[ivar]), 2.0);
    addterm += R* pow(sol[ivar], 2.0) - 2.0 * R* (previous_sol[ivar] * sol[ivar]);
  }

  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "penalty: " << addterm << std::endl;
    std::cout << "SOL: " << InnerSolver->get_lb() - addterm << std::endl;
  }

if(f_Block->is_feasible()){
  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "IS_FEASIBLE_SOL" << std::endl;
  }
        if( f_max ){
              if(value_FUNCTION > best_bound){
                best_bound = value_FUNCTION;
              } 
              best_solutions.push_back(std::pair( f_Block->get_Solution() , value_FUNCTION ));
              std::sort(best_solutions.begin(), best_solutions.end(), [](auto &left, auto &right) {
                  return left.second < right.second;
              });
          } else {
              if(value_FUNCTION < best_bound){
                  best_bound = value_FUNCTION;
              }
              best_solutions.push_back(std::pair( f_Block->get_Solution() , value_FUNCTION ));
              std::sort(best_solutions.begin(), best_solutions.end(), [](auto &left, auto &right) {
                return left.second < right.second;
              });
          }

  if(size(best_solutions) > get_int_par( intMaxSol ))
    best_solutions.erase(best_solutions.begin());
       
 } else {
  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
    std::cout << "IS_INFEASIBLE_SOL" << std::endl;
  }
 } 
    
  if ( iters >= 0 ){
    is_the_same = true;
    for( int ivar = 0 ; ivar < NumBinStatVar ; ++ivar ){
     if (!(sol[ivar] < previous_sol[ivar]+1e-6 && sol[ivar] > previous_sol[ivar]-1e-6)) { //1e-12
            is_the_same = false;
            break;
      }
    }

    //if(f_Block->is_feasible())
      //is_the_same = true;
  }

  if ( iters >= 0 ){
    is_integer_sol = true;
    integer_viol = -Inf<double>();
    integer_viol_sum = 0.0;
    for( int ivar = 0 ; ivar < NumBinStatVar ; ++ivar ){
      integer_viol_sum += std::min(std::abs(sol[ivar]),std::abs(1-sol[ivar]));
      if(std::min(std::abs(sol[ivar]),std::abs(1-sol[ivar])) > integer_viol)
        integer_viol = std::min(std::abs(sol[ivar]),std::abs(1-sol[ivar]));
    }
    if(integer_viol >= 1e-3)
      is_integer_sol = false;
  }

  iters++;

 // if iBCopy == false, bring back the inner Block to its original objective
 // "like if nothing ever happened",
 if( ! iBCopy ) {
  for( auto lbf : v_LBF )
   lbf->cleanup_inner_objective();
   /* and re-enable the Modification from the UpdateSolver - no longer needed
   for( auto us : v_US )
    us->inhibit_Modification( false );
    */
  }

   unlock();  // unlock the mutex

   // because the inner Solver is solving the dual of the original Block,
  // the unbounded an unfeasible return states have to be exchanged
  if( res == kUnbounded ){
    res = kInfeasible;
    if( logVerb - get_int_par( intLogVerb ) >= 2 ){
      std::cout << "\nINFEASIBLE!\n";
      std::cout << "NUMBER ITERS: " << iters-1 << "\n";
      std::cout << "EXIT CODE: " << res << "\n\n";
    }
    return( res );
  } else {
    if( res == kInfeasible ){
      res = kUnbounded;
     if( logVerb - get_int_par( intLogVerb ) >= 2 ){
        std::cout << "\nUNBOUNDED!\n";
        std::cout << "NUMBER ITERS: " << iters-1 << "\n";
        std::cout << "EXIT CODE: " << res << "\n\n";
      }
      return( res );
    }
  }

  if( is_integer_sol ){
    if( logVerb - get_int_par( intLogVerb ) >= 2 ){
      std::cout << "R = " << R << std::endl;
      std::cout << "IS_INTEGER_SOL = TRUE\n"; 
      std::cout << "NUMBER ITERS: " << iters-1 << "\n";
    }
  }

  if( is_the_same ){
  if( logVerb - get_int_par( intLogVerb ) >= 2 ){
      std::cout << "R = " << R << std::endl;
      std::cout << "IS_THE_SAME = TRUE\n"; 
      if( is_integer_sol ){
        std::cout << "IS_INTEGER_SOL = TRUE\n"; 
      }
      std::cout << "NUMBER ITERS: " << iters-1 << "\n";
      std::cout << "LB: " << InnerSolver->get_lb() << "\n";
      std::cout << "UB: " << InnerSolver->get_ub() << "\n";
      std::cout << "Best Feasible solution: " << get_best_bound() << std::endl;
  }

    return( res );
  }
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

void PrimalProximalHeur::add_penalty_terms()
{

 Index pos = 0;
 Index index = 0;
 penalty = 0.0;
 index = 0;

for( const auto & sbi : f_Block->get_nested_Blocks() ) {

  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );
 
 for( Index ivar = 0 ; ivar < pos_id_sbi[index] ; ++ivar ){  
  	 auto indexz = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
        )->get_function())->is_active(idx_to_var_sbi1[index][ivar].second);  

       if(static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->                                                                                                       
            get_num_active_var() > indexz){                                                                            
         auto addval1 = idx_to_var_sbi1[index][ivar].first;
         auto addval2 = idx_to_var_sbi2[index][ivar].first;  
         static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function()
             )->modify_term( indexz , addval1 - R * 2.0 * previous_sol[pos] , addval2 + R, mp );                                                                                                                                                                                              
       } 

       else {

      if(R > 0.0){                                                                                                                                                                                        
	   	    static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function()                                                                                                                                             
					    )->add_variable( idx_to_var_sbi1[index][ivar].second, -R * 2.0 * previous_sol[pos] , R , mp );           
       }}

    pos++;
 }

sbi->close_channel(chnl);
index += 1;
}
}

/*--------------------------------------------------------------------------*/

void PrimalProximalHeur::remove_penalty_terms()
{

 Index pos = 0;
 Index index = 0;
 double value_FUNCTION1 = 0.0;

 for( const auto & sbi : f_Block->get_nested_Blocks() ) {

  const auto chnl = sbi->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );

   for( Index ivar = 0 ; ivar < pos_id_sbi[index] ; ++ivar ){   
        if(static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->
              get_num_active_var() > static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
                )->get_function())->is_active(idx_to_var_sbi1[index][ivar].second)){
        
          auto indexz = static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()
              )->get_function())->is_active(idx_to_var_sbi1[index][ivar].second);                                                                                                                      
          auto addval1 = idx_to_var_sbi1[index][ivar].first; 
          auto addval2 = idx_to_var_sbi2[index][ivar].first; 
          static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function()                                                                                                                                             
              )->modify_term( indexz , addval1 , addval2, mp );
        } else {
          static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective()                                                                                                                      
                                             )->get_function()                                                                                                                                             
					    )->add_variable( idx_to_var_sbi1[index][ivar].second, idx_to_var_sbi1[index][ivar].first , idx_to_var_sbi2[index][ivar].first , mp );
        }
    pos++;
  }   

Funct_sbi[index]->compute(true);
value_FUNCTION1 += Funct_sbi[index]->get_value();  
sbi->close_channel(chnl);
index += 1;
}

if( logVerb - get_int_par( intLogVerb ) >= 2 ){
  std::cout << "value_FUNCTION1 = " << value_FUNCTION1 << std::endl;
}

value_FUNCTION = 0.0;
index = 0;

for( const auto & sbi : f_Block->get_nested_Blocks() ) {
  static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->compute(true);
  value_FUNCTION += static_cast< p_DQF >( static_cast< p_FRO >( sbi->get_objective())->get_function())->get_value();
}

if( logVerb - get_int_par( intLogVerb ) >= 2 ){
  std::cout << "value_FUNCTION = " << value_FUNCTION << std::endl;
}

}

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
  } else { 
    for( const auto & sbi : f_Block->get_nested_Blocks() ) {
      if( mod->get_Block() == sbi && ! changed_penalties ) {
        if( auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod.get() ) ) 
          reload = true;
        else
          check_feasibility = true;
      }
    }
  }
 }

 if( logVerb - get_int_par( intLogVerb ) >= 2 )
  std::cout << "check_feasibility = " << check_feasibility << std::endl;

 if(reload){
  if( logVerb - get_int_par( intLogVerb ) >= 2 )
    std::cout << "reload..." << std::endl;
  initialize();
 }
 int indexSol = 0;

 if(check_feasibility)
  for( auto sol : best_solutions ){
    sol.first->write( f_Block );
    if(! f_Block->is_feasible()){
      best_solutions.erase(best_solutions.begin() + indexSol);
    }
    indexSol++;
 }

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( v_mod.empty() )  // no Modification coming directly from f_Block
  return;                 // all done
  
 }  // end( PrimalProximalHeur::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------- End File PrimalProximalHeur.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
