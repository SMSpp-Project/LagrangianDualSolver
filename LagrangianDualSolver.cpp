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

using p_BC = BlockConfig *;
using p_BC = BlockSolverConfig *;

using p_SConf_p_p = SimpleConfiguration< std::pair< Configuration * ,
						    Configuration * > > *;

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
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for string parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::str_pars_map = {
 { "str_LDSlv_ISName" , LagrangianDualSolver::str_LDSlv_ISName } ,
 { "str_LDBlck_BCfg" , LagrangianDualSolver::str_LDBlck_BCfg } ,
 { "str_LDBlck_BSlvCfg" , LagrangianDualSolver::str_LDBlck_BSlvCfg }
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
 "BundleSolver" ,  // str_LDSlv_ISName
 "" ,              // str_LDBlck_BCfg
 ""                // str_LDBlck_BSlvCfg
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
  cleanup_LagrDual();
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
  throw( std::invalid_argument( "LagrangianDualSolver: no sub-Block" ) );

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
	      "LagrangianDualSolver: mixed min/max sub-Block Objective" ) );
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
       "LagrangianDualSolver: static constraint not a FRowConstraint" ) );
   }
  }

 static_cons = NumVar;

 // sort the static constraints-->Lagrangian-variables dictionary
 // this is not necessary for the Lagrangian-variables-->static constraints
 // one since it's surely sorted already
 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 // count and check the dynamic FRowConstraint- - - - - - - - - - - - - - - -
 for( const auto & el : f_Block->get_dynamic_constraints() ) {
  // Singles lists
  if( un_any_thing_0( FRowConstraint , el , { NumVar += var.size() } ) )
   continue;
  // Vectors of lists
  if( un_any_thing_1( FRowConstraint , el ,
                      {
		       for( auto & lel: var )
			NumVar += lel.size();
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
         "LagrangianDualSolver: dynamic constraint not a FRowConstraint" ) );
  }

 // create the Lagrangian Dual and its sub-Block- - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LagrDual = new AbstractBlock;  // create the AbstractBlock

 // resize the sub-Block dictionary and the pointers to the LagBFunction
 blck_to_idx.resize( f_nsb );
 f_LBF.resize( f_nsb );

 // first loop: create the sub-Block and their LagBFunction - - - - - - - - -
 // here comes the crucial decision: copy the sub-Block or "evict" them
 f_US.resize( f_nsb );  // meanwhile handle the UpdateSolver

 for( Index i = 0 ; i < f_nsb ; ++i ) {
  auto sbi = new AbstractBlock;
  blck_to_idx[ i ].first = sbi;
  blck_to_idx[ i ].second = i;

  Block * csbi;
  if( iBCopy ) {  // copying the sub-Block
   csbi = sb[ i ]->get_R3_Block( nullptr );  // the copy R3B
   // immediately register an UpdateSolver to the original sub-Block so
   // that any Modification to the original sub-Block is immediately
   // forwarded to the copy
   f_US[ i ] = new UpdateSolver( csbi );
   sb[ i ]->register_Solver( f_US[ i ] );
   }
  else {         // evicting the sub-Block
   csbi = sb[ i ];  // use the original sub-Block; note that its father
                    // will be changed when used in LagBFunction constructor
   f_US[ i ] = new UpdateSolver( f_Block , 2 );
   csbi->register_Solver( f_US[ i ] );
   // immediately register an UpdateSolver to the original sub-Block so
   // that any Modification to the original sub-Block is immediately
   // forwarded to the former father as it it were still its son; note that
   // the UpdateSolver forwards, as opposed to map_*, the Modification
   }

  auto lbfi = new LagBFunction( csbi );
  f_LBF[ i ] = lbfi;
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
   auto lhs = con.get_lhs();
   auto rhs = con.get_rhs();

   if( ( ( lhs == -INFshift ) && ( rhs == INFshift ) ) || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( *(Lit++) , 0 );
    ++LTit;
    return;
    }

   if( ( lhs > -INFshift ) && ( rhs < INFshift ) && ( lhs != rhs ) )
    throw( std::invalid_argument(
     "LagrangianDualSolver: ranged static constraints not supported yet" ) );

   // define the sign constraints on the multiplier (if any)
   if( f_convex ) {  // for a max problem
    if( lhs == -INFshift )                // a <= constraint 
     Lit->is_positive( true , eNoMod );   // ==> a >= multiplier
    else
     if( rhs == INFshift )                // a >= constraint 
      Lit->is_negative( true , eNoMod );  // ==> a <= multiplier     
    }
   else {            // for a min problem
    if( lhs == -INFshift )                // a <= constraint 
     Lit->is_negative( true , eNoMod );   // ==> a <= multiplier
    else
     if( con_rhs == INFshift )            // a >= constraint 
      Lit->is_positive( true , eNoMod );  // ==> a >= multiplier
    }

   // write the coefficient in the objective
   *(objit++) = std::make_pair( *(Lit++) , rhs == INFshift ? lhs : rhs );

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
   auto lhs = con.get_lhs();
   auto rhs = con.get_rhs();

   if( ( ( lhs == -INFshift ) && ( rhs == INFshift ) ) || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( *(Lit++) , 0 );
    ++LTit;
    return;
    }

   if( ( lhs > -INFshift ) && ( rhs < INFshift ) && ( lhs != rhs ) )
    throw( std::invalid_argument(
    "LagrangianDualSolver: ranged dynamic constraints not supported yet" ) );

   // define the sign constraints on the multiplier (if any)
   if( f_convex ) {  // for a max problem
    if( lhs == -INFshift )                // a <= constraint 
     Lit->is_positive( true , eNoMod );   // ==> a >= multiplier
    else
     if( rhs == INFshift )                // a >= constraint 
      Lit->is_negative( true , eNoMod );  // ==> a <= multiplier     
    }
   else {            // for a min problem
    if( lhs == -INFshift )                // a <= constraint 
     Lit->is_negative( true , eNoMod );   // ==> a <= multiplier
    else
     if( rhs == INFshift )                // a >= constraint 
      Lit->is_positive( true , eNoMod );  // ==> a >= multiplier
    }

   // write the coefficient in the objective
   *(objit++) = std::make_pair( *(Lit++) , rhs == INFshift ? lhs : rhs );

   // split the linear constraint among the sub-Block
   split_constraint( con , *(LTit++) );
   };

  // finally apply the lambda to all dynamic constraints
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_dynamic( el , scan , un_any_type< FRowConstraint >() );
  }

 // sort the dynamic constraints-->Lagrangian-variables dictionary
 // this must not be done for the Lagrangian-variables-->dynamic constraints
 // one since the mapping is positional
 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

 // release the Block- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // all required data has been read, now we only need to finish constructing
 // the Lagrangian Dual
 
 if( ! owned )
  f_Block->unlock( f_id );

 // pass the Lagrangian terms to the corresponding LagBFunction - - - - - - -
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

 // if a BlockConfig is not present but a name is, load it from file
 if( ( ! f_BCfg ) && ( ! f_BCfg_name.empty() ) ) {
  auto BC = Configuration::deserialize( f_BCfg_name );
  f_BCfg = dynamic_cast< BlockConfig * >( BC );
  if( ! f_BCfg ) {
   delete BC;
   throw( std::invalid_argument(
		   "LagrangianDualSolver: invalid BlockConfig from file" ) );

   }
  }

 // if a BlockConfig is present, apply() it
 if( f_BCfg )
  f_BCfg->apply( LagrDual );

 // if a BlockSolverConfig is not present but a name is, load it from file
 if( ( ! f_BSlvCfg ) && ( ! f_BSlvCfg_name.empty() ) ) {
  auto BSC = Configuration::deserialize( f_BSlvCfg_name );
  f_BSlvCfg = dynamic_cast< BlockSolverConfig * >( BSC );
  if( ! f_BSlvCfg ) {
   delete BCS;
   throw( std::invalid_argument(
	    "LagrangianDualSolver: invalid BlockSolverConfig from file" ) );

   }
  }

 // if a BlockConfig is present, apply() it
 if( f_BSlvCfg )
  f_BSlvCfg->apply( LagrDual );

 // register the inner Solver to the Lagrangian Dual Block
 register_inner_Solver();

 }  // end( LagrangianDualSolver::set_Block )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( int_LDSlv_iBCopy ):
   iBCopy = bool( value );
   break;
  default:
   InnerSolver->set_par(  int_par_lds( par ) , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , double value )
{
 InnerSolver->set_par( dbl_par_lds( par ) , value );
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
  case( str_LDBlck_BCfg ):
   f_BCfg_name = value;
   break;
  case( str_LDBlck_BSlvCfg ):
   f_BSlvCfg_name = value;
   break;
  default:
   InnerSolver->set_par( str_par_lds( par ) , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_ComputeConfig( ComputeConfig * scfg )
{
 if( ! scfg ) {  // factory reset
  delete f_BCfg;
  f_BCfg = nullptr;
  delete f_BSlvCfg;
  f_BSlvCfg = nullptr;
  ThinComputeInterface::set_ComputeConfig();
  return;
  }

 // first of all check if the inner Solver is changing and act upon it
 for( const auto & pair : scfg->str_pars )
  if( pair.first == "str_LDSlv_ISName" )
   set_par( str_LDSlv_ISName , pair.second );

 // now call the base ThinComputeInterface to do the bulk of work; note
 // that str_LDSlv_ISName is called twice (if ever), but set_par( string )
 // checks if the class remains the same and does nothing
 ThinComputeInterface::set_ComputeConfig( scfg );

 // finally, take care of extra_Configuration (if any)
 if( ! scfg->f_extra_Configuration )
  return;
 
 if( auto scpp = dynamic_cast< p:SConf_p_p >(
			                   scfg->f_extra_Configuration ) ) {
  if( scpp->f_value.first ) {
   clear_inner_BlockSolverConfig();
   f_BSlvCfg = dynamic_cast< p_BSC >( scpp->f_value.first );
   if( ! f_BSlvCfg )
    throw( std::invalid_argument(
       "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config.fist"
				 ) );
   }

  if( scpp->f_value.second ) {
   clear_inner_BlockConfig();
   f_BCfg = dynamic_cast< p_BC >( scpp->f_value.second );
   if( ! f_BCfg )
    throw( std::invalid_argument(
       "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config.second"
				 ) );
   }
  return;
  }

 if( auto BSC = dynamic_cast< p_BSC >( scfg->f_extra_Configuration ) ) {
  clear_inner_BlockSolverConfig();
  f_BSlvCfg = BSC;
  return;
  }
  
 if( auto BC = dynamic_cast< p_BC >( scfg->f_extra_Configuration ) ) {
  clear_inner_BlockConfig();
  f_BCfg = BC;
  return;
  }

 throw( std::invalid_argument(
         "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config" ) );
   
 }  // end( LagrangianDualSolver::set_ComputeConfig )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::compute( bool changedvars )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );
 
 return( InnerSolver->compute() );

 }  // end( LagrangianDualSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_var_solution( Configuration * solc )
{
 if( ! LagrDual )
  throw( std::logic_error(
    "LagrangianDualSolver::get_var_solution: Lagrangian Dual not formed" ) );

 // first ensure that the optimal convex multipliers are written as the
 // important_linearization_coefficients() of the LagBFunction
 InnerSolver->get_dual_solution();

 // define a lambda that does the solution (computation and) retrieval
 // for a specific sub-Block
 auto [ this ] getsoli( Index i ) -> void {
  auto szi = v_LBF[ i ]->get_dflt_int_par( C05Function::intGPMaxSz );
  if( ! szi )
   throw( std::invalid_argument(
           "LagrangianDualSolver::get_var_solution: no Solution stored" ) );
  auto & lc = v_LBF[ i ]->get_important_linearization_coefficients();
  if( lc.empty() )
   throw( std::invalid_argument(
       "LagrangianDualSolver::get_var_solution: no coefficients stored" ) );

  Index pos;
  if( lc.size() == 1 )        // the solution is already computed
   pos = lc.front().first();  // this is its position
  else {
   // the solution need be computed: find a free spot where to put it
   for( pos = 0 ; pos < szi ; ++pos )
    if( ! v_LBF[ i ]->is_linearization_there( pos ) )
     break;

   // if no free spot can be found, put it anywhere
   if( pos == szi )
    pos = 0;

   // compute the solution and store it in the global pool
   v_LBF[ i ]->store_combination_of_linearizations( lc , pos );

   // now the important linearization is in the pool, recall this
   LinearCombination nlc( { pos , 1 } );
   v_LBF[ i ]->set_important_linearization( nlc );
   }

  // if necessary move the solution back from the global pool to the Block
  v_LBF[ i ]->global_pool_to_block( pos );

  // if sub-Block is a copy, map_back the solution to the original
  if( iBCopy )
    f_Block->get_nested_Block( i )->map_back_solution(
			    v_LBF[ i ]->get_nested_Block( 0 ) , nullptr );
  };

 auto SCvi = dynamic_cast< SimpleConfiguration< std::vector< int > >
			   >( solc );
 if( SCvi ) {
  assert( std::is_ordered( SCvi->value.begin() , SCvi->value.end() ) );
  if( ( SCvi->value.front() < 0 ) || ( SCvi->value.back() >= f_nsb ) )
   throw( std::invalid_argument(
  "LagrangianDualSolver::get_var_solution: wrong indices in solc->value" ) );

  for( auto el : SCvi->value )
   getsoli( el );
  }
 else
  for( Index i = 0 ; i < f_nsb ; ++i )
   getsoli( i );
 
 }  // end( LagrangianDualSolver::get_var_solution() )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_dual_solution( Configuration * solc )
{
 if( ! LagrDual )
  throw( std::logic_error(
    "LagrangianDualSolver::get_var_solution: Lagrangian Dual not formed" ) );

 bool get_duals = true;
 if( solc ) {
  auto SCvp = dynamic_cast< SimpleConfiguration< std::vector< Configuration *
							      > > >( solc );
  if( ! SCvp )
   throw( std::logic_error(
      "LagrangianDualSolver::get_var_solution: wrong Configuration type" ) );

  if( SCvp->value.size() < f_nsb )
   throw( std::logic_error(
         "LagrangianDualSolver::get_var_solution: solc->value too short" ) );

  if( SCvp->value.size() == f_nsb )
   get_duals = false;

  // get the dual solution of the constraints inside the sub-Block
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   auto LSBi = v_LBF[ i ]->get_nested_Block( 0 );
   if( LSBi->get_registered_solvers().empty() )
    continue;

   auto SBSi = dynamic_cast< CDASolver * >(
                                    LSBi->get_registered_solvers().front() );
   if( ! SBSi )
    continue;
   SBSi->get_dual_solution( SCvp->value[ i ] );
   if( iBCopy )  // the sub-Block is a copy
    f_Block->get_nested_Block( i )->map_back_solution( LSBi , nullptr ,
						       SCvp->value[ i ] );
   }
  }

 if( get_duals ) {  // get the dual solution of the relaxed constraints
  // get the static part
  auto Ls = LagrDual->get_static_variable_v< ColVariable >( "Lambda_s" );
  auto Lsit = Ls->begin();
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el , [ & Lsit ]( FRowConstraint & con ) -> void {
                              con.set_dual( (Lsit++)->get_value() );
                              } , un_any_type< FRowConstraint >() );
  // get the dynamic part
  auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );
  auto Ldit = Ld->begin();
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_static( el , [ & Lsit ]( FRowConstraint & con ) -> void {
                              con.set_dual( (Ldit++)->get_value() );
                              } , un_any_type< FRowConstraint >() );
  } 
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
  case( str_LDSlv_ISName ):   return( ISName );
  case( str_LDBlck_BCfg ):    return( f_BCfg_name );
  case( str_LDBlck_BSlvCfg ): return( f_BSlvCfg_name );
  }

 return( InnerSolver->get_str_par( str_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
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
  return( InInf );                // it is not there

 --it;  // the previous group is the one it belongs to

 // first element of the constraint group
 //const auto first = std::get< 0 >( *it );
 //auto distance = std::distance( first , con );
 auto dist = std::distance( first , std::get< 0 >( *it ) );

 if( ( dist >= 0 ) && ( Index( dist ) < std::get< 2 >( *it ) ) )
  return( std::get< 1 >( *it ) + Index( dist ) );  // it belongs to this group

 return( InINF );  // it doesn't exist

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

 return( InInf );
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
 auto lf = dynamic_cast< const p_LP >( con.get_function() );
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

void cleanup_LagrDual( void )
{
 if( ! LagrDual )  // nothing to be cleaned up
  return;          // all done

 // first detach the inner Solver
 unregister_inner_Solver();

 // unregister and delete the UpdateSolver
 const auto & lsb = LagrDual->get_nested_Blocks();
 for( Index i = 0 ; i < f_nsb ; ++i )
  lsb[ i ]->unregister_Solver( f_US[ i ] , true );

 f_US.clear();

 // cleanup the BlockConfig (but keep it, so use a clear()-ed clone())
 if( f_BCfg ) {
  auto BC = f_BCfg->clone();
  BC->clear();
  BC->apply( LagrDual );
  delete BC;
  }

 // cleanup the BlockSolverConfig (but keep it, so use a clear()-ed clone())
 if( f_BSlvCfg ) {
  auto BSC = f_BSlvCfg->clone();
  BSC->clear();
  BSC->apply( LagrDual );
  delete BSC;
  }

 // if necessary put back the sub_Block
 if( ! iBCopy ) {
  bool owned = f_Block->is_owned_by( f_id );
  if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
   throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

  // remove the sub-Block from the LagBFunction (but do not delete it)
  for( auto el : f_LBF ) {
   el->set_inner_block( nullptr , false );

   // re-attach the sub-Block to their original father
   const auto & sb = f_Block->get_nested_Blocks();
   for( Index i = 0 ; i < f_nsb ; ++i )
    sb[ i ]->set_f_Block( f_Block );
   }

  if( ! owned )
   f_Block->unlock( f_id );
  }

 delete LagrDual;
 LagrDual = nullptr;

 }  // end( cleanup_LagrDual )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::guts_of_destructor( void )
{
 unregister_inner_Solver();
 delete InnerSolver;
 InnerSolver = nullptr;
 clear_inner_BlockConfig();
 clear_inner_BlockSolverConfig();
 cleanup_LagrDual();
		 
 }  // end( LagrangianDualSolver:guts_of_destructor )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::process_outstanding_Modification( void )
{
 // !TODO: to be done


}  // end( LagrangianDualSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------- End File LagrangianDualSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
