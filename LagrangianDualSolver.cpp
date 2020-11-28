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

using LinearCombination = C05Function::LinearCombination;
using c_LinearCombination = C05Function::c_LinearCombination;

using dual_pair = LagBFunction::dual_pair;
using v_dual_pair = std::vector< dual_pair >;
using v_c_dual_pair = const v_dual_pair;

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
 "intLPar1" ,
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of double parameters names

const std::vector< std::string > LagrangianDualSolver::dbl_pars_str = {
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for int parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::int_pars_map = {
 { "intLPar1" , LagrangianDualSolver::intLPar1  } ,
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for double parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::dbl_pars_map = {
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default int parameters

const std::vector< int > LagrangianDualSolver::dflt_int_par = {
  0  // intLPar1
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default double parameters

const std::vector<double> LagrangianDualSolver::dflt_dbl_par = {
 };

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF LagrangianDualSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_Block( Block * block )
{
    
 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done
  
 // the block does not contain any variable
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
 if( f_Block->get_static_variables().size() )
  throw( std::logic_error( "static Variable are not allowed" ) );
    
 if( f_Block->get_dynamic_variables().size() )
  throw( std::logic_error( "dynamic Variable are not allowed" ) );
    
 // children are required to exist   - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const auto & sb = f_Block->get_nested_Blocks();
 if( sb.empty() )
  throw( std::logic_error( "children are required to exist" ) );

 // the objective function of the block must be a LinearFunction- - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( !f_Block->get_objective() )  {

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::logic_error( "the objective is not a real function" ) );
  }

 // create a father block as an abstract one (it is a copy of
 // f_Block), if LDSPar1 is true make a copy of
 // the r3 block type
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 AbstractBlock *NB; // the copy of the father block
 LagBvect.resize( sb.size() ); // LagBFunction

 if( LPar1 ) { // using the r3 blocks

  AbstractBlock *r3b = dynamic_cast< AbstractBlock * >(
	   f_Block->get_R3_Block( nullptr , new AbstractBlock() ) );

  for( Index i = 0 ; i < sb.size() ; ++i ) {
   sb[ i ]->register_Solver( new UpdateSolver( f_Block ) );
   LagBvect[ i ]->set_inner_block( r3b->get_nested_Blocks()[ i ] );
   }

  NB = dynamic_cast< AbstractBlock * >( r3b );
  }
 else {
  NB = dynamic_cast< AbstractBlock * >( f_Block );
  for( Index i = 0 ; i < sb.size() ; ++i )
   LagBvect[ i ]->set_inner_block( sb[ i ] );
  }

 // construct the Lagrangian variables: lambda
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Vec_any sc = f_Block->get_static_constraints();
 auto LagVars = new std::vector< ColVariable >( sc.size() );

 // set the variables in the father abstract block
 NB->add_static_variable( *LagVars , "" , true );

 // the relaxed constraints are in the form of
 // sum_h A_h x_h = b, so the objective function
 // of the father block has to be the linear function
 // lambda * b;

 // construct the objective function of the abstract block
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Function::Vec_FunctionValue bVector( sc.size() );

 double bterm;  // read the b vector from
 for( Index i = 0 ; i < sc.size() ; ++i ) {
  if( un_any_const_static( sc[ i ] ,
	 [&bterm]( LinearConstraint & cnst ) {
	  bterm = (cnst.get_linear_function())->get_constant_term(); } ,
	  un_any_type< LinearConstraint >() ) );
  bVector[i] = bterm;
  }

 LinearFunction::v_coeff_pair vars( sc.size() );
 for( Index i = 0; i < sc.size() ; ++i )
  vars[ i ] = std::make_pair( &(*LagVars)[i] , bVector[i] );

 NB->set_objective( new FRealObjective( NB ,
			new LinearFunction( std::move( vars ) ) ) , eNoMod );

 // create the children of NB as many as the children
 // of f_Block
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index h = 0 ; h < sb.size() ; ++h ) {

  (LagBvect[ h ]->get_inner_block())->set_f_Block( NB );

  // LagPairs are the pairs to be relaxed in the
  // h-th child

  v_dual_pair LagPairs( 0 );
  for( Index i = 0; i < sc.size() ; ++i ) {

   LinearFunction * relaxed_function;
   // relaxed_function represents the linear function
   // a_i^T x_i, where x_i contains the variables
   // of the component only

   LinearFunction *gi; // the i-th constraint in the father block
   un_any_static( sc[i] , [ & ]( LinearFunction & linfun ) { gi = & linfun; } ,
		  un_any_type<LinearFunction>() );

   // get the pairs if the constraint gi
   const auto & rp = gi->get_v_var();

   // get the objective of the child
   const auto frobj = LagBvect[ h ]->get_inner_block()->get_objective< FRealObjective >();
   auto obj = dynamic_cast< p_LF >( frobj->get_function() );

   for( Index h = 0 ; h < rp.size() ; ++h ) {

	// if the variable of the constraint is active
	// in the child, this has to be added to
	// the relaxed constraint

    auto j = obj->is_active( rp[ h ].first );

    if( j < obj->get_num_active_var() ) {
     relaxed_function->add_variable( rp[ h ].first ,
    		 rp[ h ].second );
     }
    } // end scanning of variables

   // add the pair to be relaxed in the h-th child
   if( relaxed_function->get_num_active_var() )
    LagPairs.push_back( std::make_pair( &(*LagVars)[i] ,
       	  relaxed_function ) );
   }
  } // end children creation
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 // the set of "active" Variable
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = sc.size();
 LamVcblr.resize( NumVar );
 for( Index i = 0 ; i < NumVar ; ++i )
  LamVcblr[ i ] = &(*LagVars)[i];

 // allocate memory  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda.resize( NumVar );    // the point

 // !TODO: to be completed

 // the
 // bundle solver registration - - - - - - - - - - - - - - - - - - - - - - - -

 NB->register_Solver(BndSlv, true);

 }  // end( LagrangianDualSolver::set_Block )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( const idx_type par , const int value )
{
 switch( par ) {
  case( intLPar1 ):
   if( value < 0 )
    throw( std::invalid_argument( "LPar1 must be >= 0" ) );
   LPar1 = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( LagrangianDualSolver::set_par( int ) )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( const idx_type par , const double value )
{
 switch( par ) {
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( LagrangianDualSolver::set_par( double ) )

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
 #ifndef NDEBUG
  if( ! InnrSlv )
   throw( std::logic_error( "inner CDASolver not initialised yet" ) );
 #endif

 // !ToDO: to be completed
 return( InnrSlv->compute() );



 }  // end( LagrangianDualSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_var_solution( Configuration *solc )
{
 #ifndef NDEBUG
  if( ! InnrSlv )
   throw( std::logic_error( "inner CDASolver not initialised yet" ) );
 #endif
 // !TODO: to be implemented

 }  // end( LagrangianDualSolver::get_var_solution() )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_dual_solution( Configuration *solc )
{
 #ifndef NDEBUG
  if( ! InnrSlv )
   throw( std::logic_error( "inner CDASolver not initialised yet" ) );
 #endif

 // !TODO: to be implemented

 }  // end( LagrangianDualSolver::get_dual_solution() )

/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::get_int_par( const idx_type par ) const
{
 switch( par ) {
  case( intLPar1 ):
   return( LPar1 );
   break;
  default:
   return( CDASolver::get_dflt_int_par( par ) );
  }
 }  // end( LagrangianDualSolver::get_int_par )

/*--------------------------------------------------------------------------*/

double LagrangianDualSolver::get_dbl_par( const idx_type par ) const
{
 switch( par ) {
  default:
   return( CDASolver::get_dflt_dbl_par( par ) );
  }
 }  // end( LagrangianDualSolver::get_dbl_par )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::Log1( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 // !TODO: to be completed
 } // end( LagrangianDualSolver::Log1 )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::Log2( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 // !TODO: to be completed

 } // end( LagrangianDualSolver::Log2 )

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
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::guts_of_destructor( void )
{
 
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
