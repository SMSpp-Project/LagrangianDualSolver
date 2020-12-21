/*--------------------------------------------------------------------------*/
/*--------------------- File LagrangianDualSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the LagrangianDualSolver class, which implements the
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
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2020 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "LagrangianDualSolver.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LagBFunction.h"

#include "LinearFunction.h"

#include "RBlockConfig.h"

#include "UpdateSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LdsLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

#define LdsLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using p_LF = LinearFunction *;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using VarValue = Function::FunctionValue;
using c_VarValue = Function::c_FunctionValue;
using Vec_FunctionValue = Function::Vec_FunctionValue;

using Vec_VarValue = Function::Vec_FunctionValue;
using c_Vec_VarValue = Function::c_Vec_FunctionValue;

using v_coeff_pair = LinearFunction::v_coeff_pair;
using v_c_coeff_pair = LinearFunction::v_c_coeff_pair;

using LinearCombination = C05Function::LinearCombination;
using c_LinearCombination = C05Function::c_LinearCombination;

using dual_pair = LagBFunction::dual_pair;
using v_dual_pair = std::vector< dual_pair >;
using v_c_dual_pair = const v_dual_pair;

using p_AB = AbstractBlock *;
using p_LF = LinearFunction *;
using p_LBF = LagBFunction *;

using SConf_p_p = SimpleConfiguration< std::pair< Configuration * ,
						  Configuration * > >;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr VarValue NaNshift
                              = std::numeric_limits< VarValue >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

static constexpr VarValue INFshift
                               = std::numeric_limits< VarValue >::infinity();
 ///< convenience constexpr for "Infty"

static constexpr cIndex InINF = SMSpp_di_unipi_it::Inf<Index>();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register LagrangianDualSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( LagrangianDualSolver );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names

const std::vector< std::string > LagrangianDualSolver::int_pars_str = {
 "int_LDSlv_iBCopy" ,
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of double parameters names

const std::vector< std::string > LagrangianDualSolver::dbl_pars_str = {
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of string parameters names

const std::vector< std::string > LagrangianDualSolver::str_pars_str = {
 "str_LDSlv_ISName"
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for int parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::int_pars_map = {
 { "int_LDSlv_iBCopy" , LagrangianDualSolver::int_LDSlv_iBCopy  } ,
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for double parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::dbl_pars_map = {
 { "str_LDSlv_ISName" , LagrangianDualSolver::str_LDSlv_ISName  } ,
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for string parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::str_pars_map = {
 };


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default int parameters

const std::vector< int > LagrangianDualSolver::dflt_int_par = {
 0  // int_LDSlv_iBCopy
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default double parameters

const std::vector< double > LagrangianDualSolver::dflt_dbl_par = {
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default double parameters

const std::vector<double> LagrangianDualSolver::dflt_str_par = {
 "BundleSolver"  // str_LDSlv_ISName
 };

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF LagrangianDualSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_Block( Block * block )
{
 if( f_Block ) {  // changing from a previous Block- - - - - - - - - - - - - -
                  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 // lock the Block - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );


 // check conditions on the Block- - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // the Block must not contain any variable- - - - - - - - - - - - - - - - -
    
 if( ! f_Block->get_static_variables().empty() )
  throw( std::invalid_argument(
		    "LagrangianDualSolver: static Variable not allowed" ) );
    
 if( ! f_Block->get_dynamic_variables().empty() )
  throw( std::invalid_argument(
		   "LagrangianDualSolver: dynamic Variable not allowed" ) );
    
 // there must be no Objective- - - - - - - - - - - - - - - - - - - - - - - -

 if( f_Block->get_objective() )
  throw( std::invalid_argument(
			   "LagrangianDualSolver: Objective not allowed" ) );

 // children are required to exist - - - - - - - - - - - - - - - - - - - - -

 const auto & sb = f_Block->get_nested_Blocks();
 f_nsb = sb.size();
 if( ! f_nsb )
  throw( std::logic_error( "LagrangianDualSolver: no sub-Block" ) );

 // children must have a FRealObjective with a LinearFunction inside, all
 // of them must have the same "verse"- - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < f_nsb ; ++i ) {
  // generate the Objective (if not there already)
  sb[ i ]->generate_objective();
  auto osbi = dynamic_cast< cost FRealObjective * >(
						  sb[ i ]->get_objective() );
  if( ! osbi )
   throw( std::invalid_argument(
		       "LagrangianDualSolver: wrong sub-Block Objective" ) );
  if( ! i )
   f_convex = ( osbi->get_sense() == Objective::eMax );
  else
   if( f_convex != ( osbi->get_sense() == Objective::eMax ) )
    throw( std::invalid_argument(
	  "LagrangianDualSolver: sub-Block Objective with mixed min/max" ) );
  }
 
 // count and check the FRowConstraint in the Block - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = 0;

 // count and check the static FRowConstraint - - - - - - - - - - - - - - - -
 // meanwhile construct the static dictionaries

 // number of groups of static constraints
 auto scn = f_Block->get_static_constraints().size();

 // resize the static constraints<-->Lagrangian-variables dictionaries
 scon_to_idx.resize( scn );
 idx_to_scon.resize( scn );

 {
  Index pos = 0;
  for( const auto & el : f_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] = std::make_pair( & var , NumVar );
			idx_to_scon[ pos++ ] =
			 std::make_pair( NumVar++ , & var );
		        } ) )
    continue;
   // Vectors
   if( un_any_thing_1( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] =
			 std::make_pair( var.data() , NumVar );
			idx_to_scon[ pos++ ] =
			 std::make_pair( NumVar , var.data() );
			NumVar += var.size();
		        } ) )
    continue;
   // Multiarrays
   if( un_any_thing_K( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] =
			 std::make_pair( var.data() , NumVar );
			idx_to_scon[ pos++ ] =
			 std::make_pair( NumVar , var.data() );
			NumVar += var.num_elements();
		        } ) )
    continue;
   throw( std::invalid_argument(
    "LagrangianDualSolver: some static constraint not a FRowConstraint" ) );
   }
  }

 static_cons = NumVar;

 // sort the static constraints-->Lagrangian-variables dictionary
 // this is not necessary for the Lagrangian-variables-->static constraints
 // one since it's surely sorted already
 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 for( const auto & el : f_Block->get_dynamic_constraints() ) {
  // Singles lists
  if( un_any_thing_0( FRowConstraint , el , { NumVar += var.size() } ) )
   continue;
  // Vectors of lists
  if( un_any_thing_1( FRowConstraint , el ,
                      {
		       for( auto & el: var )
			NumVar += el.size();
		       } ) )
    continue;
  // Multiarrays of lists
  if( un_any_thing_K( FRowConstraint , el ,
		      {
		       auto it = var.data();
		       for( auto i = var.num_elements() ; i-- ; ++it )
			NumVar += it->size();
		       } ) )
   continue;
  throw( std::invalid_argument(
   "LagrangianDualSolver: some dynamic constraint not a FRowConstraint" ) );
  }

 // create the Lagrangian Dual and its sub-Block- - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LagrDual = new AbstractBlock;  // create the AbstractBlock

 // resize the sub-Block dictionary
 blck_to_idx.resize( f_nsb );

 // first loop: create the sub-Block and their LagBFunction - - - - - - - - -
 // here comes the crucial decision: copy the sub-Block or "evict" them

 for( Index i = 0 ; i < f_nsb ; ++i ) {
  auto sbi = new AbstractBlock;
  std::get< 0 >( blck_to_idx[ i ] ) = sbi;
  std::get< 1 >( blck_to_idx[ i ] ) = i;

  Block * csbi;
  if( LPar1 ) {  // copying the sub-Block
   csbi = sb[ i ]->get_R3_Block( nullptr );  // the copy R3B
   // immediately register an UpdateSolver to the original sub-Block so
   // that any Modification to the original sub-Block is immediately
   // forwarded to the copy
   sb[ i ]->register_Solver( new UpdateSolver( csbi ) );
   }
  else {         // evicting the sub-Block
   csbi = sb[ i ];  // use the original sub-Block; note that its father
                    // will be changed when used in LagBFunction constructor
   csbi->register_Solver( new UpdateSolver( f_Block ) );
   // immediately register an UpdateSolver to the original sub-Block so
   // that any Modification to the original sub-Block is immediately
   // forwarded to the former father as it it were still its son
   }

  auto lbfi = new LagBFunction( csbi );
  std::get< 2 >( blck_to_idx[ i ] ) = lbfi;
  auto osbi = new FRealObjective( sbi , lbfi );
  osbi->set_sense( f_convex ? Objective::eMin : Objective::eMax , eNoMod );
  sbi->set_objective( osbi );

  LagrDual->add_nested_Block( sbi );  // add the sub-Block
  }

 // sort the sub-Block dictionary by Block address
 std::sort( blck_to_idx.begin() , blck_to_idx.end() );

 // create the static and dynamic Lagrangian variables- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto Ls = new std::vector< ColVariable >( static_cons );
 auto Ld = new std::list< ColVariable >( NumVar - static_cons );

 // pass the Lagrangian variables to the Lagrangian Dual
 LagrDual->add_static_variable( *Ls , "Lambda_s" );
 LagrDual->add_dynamic_variable( *Ld , "Lambda_d" );

 // create the Objective of the Lagrangian Dual - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto lf = new LinearFunction;
 auto obj = new FRealObjective( LagrDual , lf );
 obj->set_sense( f_convex ? Objective::eMin : Objective::eMax , eNoMod );
 LagrDual->set_objective( obj , eNoMod );
 v_coeff_pair objcf( NumVar );
 
 // scan all FRowConstraints- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // meanwhile construct the linear objective function

 auto objit = objcf.begin();

 // construct the auxiliary data structure to hold the Lagrangian terms;
 // LagTerms[ i ][ h ] contains the v_coeff_pair corresponding to the
 // Lagrangian term of sub-Block h for the i-th variable

 std::vector< std::vector< v_coeff_pair > > LagTerms( NumVar );
 auto LTit = LagTerms.begin();

 // scan all static FRowConstraints - - - - - - - - - - - - - - - - - - - - -
 {
  auto Lit = Ls->begin();

  // define a lambda that does the job
  auto scan = [ & ]( FRowConstraint & con ) -> void {
   // check the LHS/RHS
   auto con_lhs = con.get_lhs();
   auto con_rhs = con.get_rhs();

   if( ( ( con_lhs == -Inf< double >() ) && ( con_rhs == Inf< double >() ) )
       || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( *(Lit++) , 0 );
    ++LTit;
    return;
    }

   if( ( con_lhs > -Inf< double >() ) && ( con_rhs < Inf< double >() ) &&
       ( con_lhs != con_rhs ) )
    throw( std::invalid_argument(
     "LagrangianDualSolver: ranged static constraints not supported yet" ) );

   // define the sign constraints on the multiplier (if any)
   if( f_convex ) {  // for a max problem
    if( con_lhs == -Inf< double >() )     // a <= constraint 
     Lit->is_positive( true , eNoMod );   // ==> a >= multiplier
    else
     if( con_rhs == Inf< double >() )     // a >= constraint 
      Lit->is_negative( true , eNoMod );  // ==> a <= multiplier     
    }
   else {            // for a min problem
    if( con_lhs == -Inf< double >() )     // a <= constraint 
     Lit->is_negative( true , eNoMod );   // ==> a <= multiplier
    else
     if( con_rhs == Inf< double >() )     // a >= constraint 
      Lit->is_positive( true , eNoMod );  // ==> a >= multiplier
    }

   // write the coefficient in the objective
   *(objit++) = std::make_pair( *(Lit++) , con_rhs == Inf< double >()
				           ? con_lhs : con_rhs );

   // split the linear constraint among the sub-Block
   split_constraint( con , *(LTit++) );
   };

  // finally apply the lambda to all static constraints
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el , scan , un_any_type< FRowConstraint >() );
  }

 // scan all dynamic FRowConstraints- - - - - - - - - - - - - - - - - - - - -
 // meanwhile construct the dynamic dictionaries

 // resize the dynamic constraints<-->Lagrangian-variables dictionaries
 dcon_to_idx.resize( NumVar - static_cons );
 idx_to_dcon.resize( NumVar - static_cons );

 {
  Index i = static_cons;
  auto Lit = Ld->begin();
  auto dc2iit = dcon_to_idx.begin();
  auto i2dcit = idx_to_dcon.begin();

  // define a lambda that does the job
  auto scan = [ & ]( FRowConstraint & con ) -> void {
   // first write the dictonaries
   *(dc2iit++) = std::make_pair( & con , i++ );
   *(i2dcit++) = & con;

   // then check the LHS/RHS
   auto con_lhs = con.get_lhs();
   auto con_rhs = con.get_rhs();

   if( ( ( con_lhs == -Inf< double >() ) && ( con_rhs == Inf< double >() ) )
       || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( *(Lit++) , 0 );
    ++LTit;
    return;
    }

   if( ( con_lhs > -Inf< double >() ) && ( con_rhs < Inf< double >() ) &&
       ( con_lhs != con_rhs ) )
    throw( std::invalid_argument(
    "LagrangianDualSolver: ranged dynamic constraints not supported yet" ) );

   // define the sign constraints on the multiplier (if any)
   if( f_convex ) {  // for a max problem
    if( con_lhs == -Inf< double >() )     // a <= constraint 
     Lit->is_positive( true , eNoMod );   // ==> a >= multiplier
    else
     if( con_rhs == Inf< double >() )     // a >= constraint 
      Lit->is_negative( true , eNoMod );  // ==> a <= multiplier     
    }
   else {            // for a min problem
    if( con_lhs == -Inf< double >() )     // a <= constraint 
     Lit->is_negative( true , eNoMod );   // ==> a <= multiplier
    else
     if( con_rhs == Inf< double >() )     // a >= constraint 
      Lit->is_positive( true , eNoMod );  // ==> a >= multiplier
    }

   // write the coefficient in the objective
   *(objit++) = std::make_pair( *(Lit++) , con_rhs == Inf< double >()
				           ? con_lhs : con_rhs );

   // split the linear constraint among the sub-Block
   split_constraint( con , *(LTit++) );
   };

  // finally apply the lambda to all dynamic constraints
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_dynamic( el , scan , un_any_type< FRowConstraint >() );
  }

 // sort the dynamic constraints-->Lagrangian-variables dictionary
 // this is not necessary for the Lagrangian-variables-->dynamic constraints
 // one since it's surely sorted already
 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

 // pass the Lagrangian terms to the LagBFunction - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index h = 0 ; h < f_nsb ; ++h ) {
  v_dual_pair dp( NumVar );  // construct the dual pairs

  Index i = 0;
  for( ; i < static_cons ; ++i ) {
   dp[ i ].first = & (*Ls)[ i ];
   dp[ i ].second = new LinearFunction( std::move( LagTerms[ i ][ h ] ) );
   }

  auto Lit = Ld->begin();
  for( ; i < NumVar ; ++i ) {
   dp[ i ].first = & (*Lit++);
   dp[ i ].second = new LinearFunction( std::move( LagTerms[ i ][ h ] ) );
   }

  auto SBi = static_cast< p_AB >( LagrDual->get_nested_Block( h ) );
  auto LBF = static_cast< p_LBF >(
		  SBi->get_objective< FRealObjective >()->get_function() );
  LBF->set_dual_pairs( std::move( lp ) );
  }

 // configure the LagrDual Solver- - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



 // finally, release the Block - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! owned )
  f_Block->unlock( f_id );
 
 }  // end( LagrangianDualSolver::set_Block )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( int_LDSlv_iBCopy ):
   if( value < 0 )
    throw( std::invalid_argument( "LPar1 must be >= 0" ) );
   iBCopy = bool( value );
   break;
  default:
   InnerSolver->set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , double value )
{
 InnerSolver->set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , std::string && value )
{
 switch( par ) {
  case( str_LDSlv_ISName ):
   value = SMSpp_classname_normalise( std::move( value ) );
   if( ISName != value ) {
    ISName = value;
    unregister_inner_Solver();
    InnerSolver = new_Solver( ISName );
    register_inner_Solver();
    }
   break;
  default:
   InnerSolver->set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_ComputeConfig( ComputeConfig * scfg )
{
 ThinComputeInterface::set_ComputeConfig( scfg );



 

 }  // end( LagrangianDualSolver::set_ComputeConfig )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::compute( bool changedvars )
{
 // !ToDO: to be completed
 return( InnerSolver->compute() );

 }  // end( LagrangianDualSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_var_solution( Configuration *solc )
{
 // !TODO: to be implemented

 }  // end( LagrangianDualSolver::get_var_solution() )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_dual_solution( Configuration *solc )
{
 // !TODO: to be implemented

 }  // end( LagrangianDualSolver::get_dual_solution() )

/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( int_LDSlv_iBCopy ): return( iBCopy );
  }

 return( InnerSolver->get_int_par( int_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualSolver::get_dbl_par( idx_type par ) const
{
 return( InnerSolver->get_dbl_par( dbl_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

const std::string &  LagrangianDualSolver::get_str_par( idx_type par ) const
{
 switch( par ) {
  case( str_LDSlv_ISName ): return( ISName );
  }

 return( InnerSolver->get_str_par( str_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::register_inner_Solver( void )
{
 if( ! LagrDual )
  return;

 LagrDual->register_Solver( InnerSolver );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::unregister_inner_Solver( void )
{
 if( LagrDual )
  LagrDual->unregister_Solver( InnerSolver , true );
 else
  delete InnerSolver;

 InnerSolver = nullptr;
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_default_inner_BlockSolverConfig( void )
{
 if( auto inner_block = get_inner_block() ) {
  auto solver_config = new RBlockSolverConfig( inner_block );
  solver_config->clear();
  solver_config->apply( inner_block );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::configure_LagrangianDualBlock( void )
{
 if( ! f_LDBConfig )
  return;

 
}

/*--------------------------------------------------------------------------*/

Index LagrangianDualSolver::index_of_static_constraint(
						 const FRowConstraint * con )
{
 if( scon_to_idx.empty() )
  return( Inf< Index >() );

 assert( std::is_sorted( scon_to_idx.begin() , scon_to_idx.end() ) );
 auto it = upper_bound( scon_to_idx.begin() , scon_to_idx.end() ,
                        std::make_tuple( con , 0 , 0 ),
                        []( auto & p1 , auto & p2 ) {
                         return( std::get< 0 >( p1 ) < std::get< 0 >( p2 ) );
                         } );

 // it now refers to the first (group of) element(s) greater than i
 if( it == scon_to_idx.begin() )  // all elements are greater
  return( Inf< Index >() );       // it is not there

 --it;  // the previous group is the one it belongs to

 // first element of the constraint group
 //const auto first = std::get< 0 >( *it );
 //auto distance = std::distance( first , con );
 auto dist = std::distance( first , std::get< 0 >( *it ) );

 if( ( dist >= 0 ) && ( Index( dist ) < std::get< 2 >( *it ) ) )
  return( std::get< 1 >( *it ) + Index( dist ) );  // it belongs to this group

 return( Inf< Index >() );  // it doesn't exist

 }  // end( LagrangianDualSolver::index_of_static_constraint )

/*--------------------------------------------------------------------------*/

Index LagrangianDualSolver::index_of_dynamic_constraint(
						  const FRowConstraint * con )
{
 assert( std::is_sorted( dcon_to_idx.begin() , dcon_to_idx.end() ) );
 auto it = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
                        std::make_pair( con , 0 ) ,
                        []( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                         } );

 if( ( it != dcon_to_idx.end() ) && ( it->first == con ) )
  return( it->second );

 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

FRowConstraint * LagrangianDualSolver::static_constraint_with_index( Index i )
{
 #ifdef NDEBUG
  if( ( i >= static_cons ) || idx_to_scon.empty() )
   throw( std::invalid_argument(
	       "LagrangianDualSolver::invalid index of static constraint" ) );

  assert( std::is_sorted( idx_to_scon.begin() , idx_to_scon.end() ) );
 #endif

 auto it = upper_bound( idx_to_scon.begin() , idx_to_scon.end() ,
                        std::make_pair( i , nullptr ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                         } );

 // it now refers to the first (group of) element(s) greater than i
 #ifdef NDEBUG
  if( it == idx_to_scon.begin() )  // all elements are greater
   throw( std::invalid_argument(
			 "LagrangianDualSolver::inconsistent idx_to_scon" ) );
 #endif

 --it;  // the previous group is the one it belongs to

 return( it->second + ( i - it->first ) );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::split_constraint( const FRowConstraint & con ,
			 std::vector< LinearFunction::v_coeff_pair > & split )
{
 auto lf = dynamic_cast< const LinearFunction * >( con.get_function() );
 if( ! lf )
  throw( std::invalid_argument(
			"LagrangianDualSolver: FRowConstraint not linear" ) );

 auto & vc = lf->get_v_var();

 split.resize( f_nsb );

 if( f_nsb == 1 ) {    // easy case: only one sub-Block, nothing to split
  split.front() = vc;
  return;
  }

 for( auto & el : split )
  el.clear();

 if( vc.empty )   // easy case: empty constraint, nothing to split
  return;

 std::vector< Index > blckidx( vc.size() );  // Block to which the var belongs
 std::vector< Index > cntr( f_nsb , 0 );

 // first pass: count the size of each split[ h ]; meanwhile save the
 // Variable-to-sub-Block-index information in blckidx to avoid computing
 // it twice;
 for( Index i = 0 ; i < vc.size() ; ) {
  auto bi = Block2Index( vc[ i ].first->get_Block() );
  blckidx[ i++ ] = bi;
  ++cntr[ bi ];
  }

 // properly size all split[ h ]; meanwhile, reset the counter
 for( Index h = 0 ; h < f_nsb ; ) {
  split[ h ].resize( cntr[ h ] );
  cntr[ h++ ] = 0;
  }
 
 // second pass: construct all split[ h ]
 for( Index i = 0 ; i < vc.size() ; ++i )
  split[ blckidx[ i ] ][ cntr[ blckidx[ i ]++ ] ] = vc[ i ];

 }  // end( LagrangianDualSolver::split_constraint )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::guts_of_destructor( void )
{
 unregister_inner_Solver();

 //!!LamVcblr.clear();

 // !TODO: to be completed

 }  // end( LagrangianDualSolver:guts_of_destructor )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::process_outstanding_Modification( void )
{
 // !TODO: to be done


}  // end( LagrangianDualSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------- End File LagrangianDualSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
