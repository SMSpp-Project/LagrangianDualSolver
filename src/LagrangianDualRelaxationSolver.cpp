/*--------------------------------------------------------------------------*/
/*--------------- File LagrangianDualRelaxationSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the classes LagrangianDualRelaxationSolver and
 * LagrangianChange.
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Filippo Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cmath>

#include <filesystem>

#include <fstream>

#include <iomanip>

#include "LagrangianDualRelaxationSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;

using OFValue = Solver::OFValue;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

// how far from the closest integer a value has to be to be fractional
static constexpr double FracEps = 1e-6;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register LagrangianDualRelaxationSolver to the Solver factory, and
// LagrangianChange to the Change one

SMSpp_insert_in_factory_cpp_0( LagrangianDualRelaxationSolver );

SMSpp_insert_in_factory_cpp_0( LagrangianChange );

/*--------------------------------------------------------------------------*/
/*---------------- METHODS OF LagrangianDualRelaxationSolver ---------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intApplyStrategy ):  f_apply_strategy = value;  break;
  case( intBranchStrategy ): f_branch_strategy = value; break;
  case( intStrongCands ):    f_strong_cands = value;    break;
  default: PrimalProximalHeur::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intApplyStrategy ):  return( eMaster );
  case( intBranchStrategy ): return( eMostFractional );
  case( intStrongCands ):    return( 10 );
  default: return( PrimalProximalHeur::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intApplyStrategy ):  return( f_apply_strategy );
  case( intBranchStrategy ): return( f_branch_strategy );
  case( intStrongCands ):    return( f_strong_cands );
  default: return( PrimalProximalHeur::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolver::int_par_str2idx(
					       const std::string & name ) const
{
 if( name == "intApplyStrategy" )
  return( intApplyStrategy );
 if( name == "intBranchStrategy" )
  return( intBranchStrategy );
 if( name == "intStrongCands" )
  return( intStrongCands );
 return( PrimalProximalHeur::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::int_par_idx2str(
						       idx_type par ) const
{
 static const std::string apply = "intApplyStrategy";
 static const std::string branch = "intBranchStrategy";
 static const std::string cands = "intStrongCands";
 switch( par ) {
  case( intApplyStrategy ):  return( apply );
  case( intBranchStrategy ): return( branch );
  case( intStrongCands ):    return( cands );
  default: return( PrimalProximalHeur::int_par_idx2str( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_par( idx_type par ,
					      std::string && value )
{
 if( par == strStrongLog )
  f_strong_log = std::move( value );
 else
  PrimalProximalHeur::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::get_dflt_str_par(
						     idx_type par ) const
{
 static const std::string empty;
 if( par == strStrongLog )
  return( empty );
 return( PrimalProximalHeur::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::get_str_par(
						     idx_type par ) const
{
 if( par == strStrongLog )
  return( f_strong_log );
 return( PrimalProximalHeur::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolver::str_par_str2idx(
					      const std::string & name ) const
{
 if( name == "strStrongLog" )
  return( strStrongLog );
 return( PrimalProximalHeur::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::str_par_idx2str(
						     idx_type par ) const
{
 static const std::string log = "strStrongLog";
 if( par == strStrongLog )
  return( log );
 return( PrimalProximalHeur::str_par_idx2str( par ) );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_global_information(
						    GlobalInformation * gi )
{
 ChangeSolver::set_global_information( gi );
 map_varToSol = nullptr;
 if( ! gi )
  return;

 if( ! gi->exists( str_VarToSol ) )
  gi->add_to_Universe< PurgedColumn >( str_VarToSol );
 map_varToSol = gi->get_from_Universe< PurgedColumn >( str_VarToSol );
 if( map_varToSol && ( ! map_varToSol->contains( str_PurgedColumns ) ) )
  map_varToSol->write( str_PurgedColumns , PurgedColumn{} );

 }  // end( LagrangianDualRelaxationSolver::set_global_information )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolver::compute( bool changedvars )
{
 const int status = PrimalProximalHeur::compute( changedvars );

 const bool done = ( status >= kOK ) &&
                   ( ( status == kLowPrecision ) || ( status < kError ) );
 const auto & sol = get_Lagrangian_convexified_solution();
 if( ( ! done ) || ( sol.size() != NumStatVar ) )
  return( status );

 // the Lagrangian Dual not solved (stopped by a budget, or kLowPrecision)
 // gives a valid bound, and a convexified solution that need not satisfy
 // the relaxed constraints: the node can still be branched on a variable
 // fractional in it, but if there is none it can be neither branched nor
 // fenced, and what the inner Solver returned is said
 const int ld_status = get_Lagrangian_initial_status();
 if( ld_status != kOK )
  if( std::none_of( sol.begin() , sol.end() , []( double v ) {
       return( std::abs( v - std::round( v ) ) > FracEps ); } ) )
   return( ld_status );

 return( kOK );

 }  // end( LagrangianDualRelaxationSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solution * LagrangianDualRelaxationSolver::get_Solution( Configuration * solc )
{
 f_Block->lock( this );
 LagrangianDualSolver::get_var_solution( solc );
 auto solution = f_Block->get_Solution( solc , false );  // loaded
 f_Block->unlock( this );
 return( solution );
 }

/*--------------------------------------------------------------------------*/

Solution * LagrangianDualRelaxationSolver::get_true_solution(
						      Configuration * solc )
{
 f_Block->lock( this );
 get_true_var_solution( solc );
 auto solution = f_Block->get_Solution( solc , false );  // loaded
 f_Block->unlock( this );
 return( solution );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

std::vector< Change * > LagrangianDualRelaxationSolver::branch( void )
{
 // the convexified Lagrangian solution, one value per binary variable of
 // the sub-Blocks: that of the sub-Blocks alone is integral whenever they
 // are, and it says nothing about where the relaxation is fractional
 const auto & sol = get_Lagrangian_convexified_solution();
 if( sol.empty() )
  throw( std::runtime_error( "LagrangianDualRelaxationSolver::branch: no "
			     "Lagrangian solution available" ) );
 if( sol.size() != NumStatVar )
  throw( std::runtime_error( "LagrangianDualRelaxationSolver::branch: the "
			     "Lagrangian solution has the wrong size" ) );

 // the candidates: the fractional variables, the most fractional first
 struct Cand {
  ColVariable * var;  // the variable
  double value;       // its value in the convexified solution
  double frac;        // its fractionality
  Index k;            // its index among the binary static variables
  Index sb;           // its sub-Block
  Index pos;          // its position among those of the sub-Block
  double cost;        // its cost in the objective of the sub-Block
  };
 std::vector< Cand > cands;
 {
  Index k = 0;
  for( Index sb = 0 ; sb < idx_to_var_sbi1.size() ; ++sb )
   for( Index pos = 0 ; pos < idx_to_var_sbi1[ sb ].size() ; ++pos , ++k ) {
    const auto & dv = idx_to_var_sbi1[ sb ][ pos ];
    const double v = sol[ k ];
    const double frac = std::abs( v - std::round( v ) );
    if( frac > FracEps )
     cands.push_back( { dv.second , v , frac , k , sb , pos , dv.first } );
    }
  }
 ++f_n_branch;

 // no fractional variable: an integral Lagrangian solution of a
 // well-terminated Lagrangian Dual satisfies the relaxed constraints, i.e.,
 // it solves the node, which is then fenced by bound and never branched;
 // were it to be branched, the two children would be the same node, hence
 // this is an error, of the termination of the inner Solver (say, dblNZEps
 // too large)
 if( cands.empty() )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::branch: the "
			   "Lagrangian solution is integral, the relaxed "
			   "constraints are violated beyond the termination "
			   "tolerance of the Lagrangian Dual (dblNZEps)" ) );

 std::stable_sort( cands.begin() , cands.end() ,
		   []( const Cand & a , const Cand & b ) {
		    return( a.frac > b.frac ); } );

 switch( f_branch_strategy ) {
  case( eMostFractional ):
   return( branchings( cands.front().var , cands.front().value ) );

  case( eStrongBranching ): {
   if( f_strong_cands <= 0 )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::branch: "
				  "intStrongCands must be positive" ) );
   if( cands.size() > Index( f_strong_cands ) )
    cands.resize( f_strong_cands );

   // the improvement of a child bound over that of the node, on the side
   // of the relaxation; at least eps, so that a candidate with one child
   // not improving still ranks by the other one
   const double node = valid_bound;
   const double eps = 1e-6 * std::max( std::abs( node ) , 1.0 );
   auto gain = [ & ]( double child ) {
    const double g = f_max ? node - child : child - node;
    return( std::max( g , eps ) );
    };

   Index best_c = 0;
   double best_score = -1;
   std::vector< std::array< double , 3 > > evals( cands.size() );
   for( Index c = 0 ; c < cands.size() ; ++c ) {
    auto changes = branchings( cands[ c ].var , cands[ c ].value );
    const double down = child_bound( changes[ 0 ] );
    const double up = child_bound( changes[ 1 ] );
    for( auto ch : changes )
     delete ch;

    const double score = gain( down ) * gain( up );
    evals[ c ] = { down , up , score };
    if( f_log && ( logVerb >= 1 ) )
     *f_log << "LagrangianDualRelaxationSolver::branch: candidate " << c
	    << " value " << cands[ c ].value << " down " << down << " up "
	    << up << " score " << score << std::endl;
    if( score > best_score ) {
     best_score = score;
     best_c = c;
     }
    }

   // the data of the candidates, for a branching rule to learn from them
   // [see strStrongLog]
   if( ! f_strong_log.empty() ) {
    const bool is_new = ! std::filesystem::exists( f_strong_log );
    std::ofstream out( f_strong_log , std::ios::app );
    if( ! out )
     throw( std::runtime_error( "LagrangianDualRelaxationSolver::branch: "
				"cannot open " + f_strong_log ) );
    if( is_new )
     out << "call,node,fixed,rank,sb,pos,value,lagr,frac,cost,down,up,"
	 << "score,chosen" << std::endl;
    Index n_fixed = 0;
    for( const auto & sbd : idx_to_var_sbi1 )
     for( const auto & dv : sbd )
      if( dv.second->is_fixed() )
       ++n_fixed;
    const auto & lagr = get_Lagrangian_initial_solution();
    out << std::setprecision( 10 );
    for( Index c = 0 ; c < cands.size() ; ++c ) {
     const auto & cd = cands[ c ];
     const auto n_sb = idx_to_var_sbi1[ cd.sb ].size();
     out << f_n_branch << "," << node << ","
	 << double( n_fixed ) / NumStatVar << "," << c << "," << cd.sb << ","
	 << ( n_sb > 1 ? double( cd.pos ) / ( n_sb - 1 ) : 0.0 ) << ","
	 << cd.value << "," << ( cd.k < lagr.size() ? lagr[ cd.k ] : cd.value )
	 << "," << cd.frac << "," << cd.cost << "," << evals[ c ][ 0 ] << ","
	 << evals[ c ][ 1 ] << "," << evals[ c ][ 2 ] << ","
	 << ( c == best_c ? 1 : 0 ) << std::endl;
     }
    }

   return( branchings( cands[ best_c ].var , cands[ best_c ].value ) );
   }

  default:
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::branch: "
				 "unknown intBranchStrategy" ) );
  }

 }  // end( LagrangianDualRelaxationSolver::branch )

/*--------------------------------------------------------------------------*/

std::vector< Change * > LagrangianDualRelaxationSolver::branchings(
				      ColVariable * var , double value ) const
{
 const std::vector< AbstractPath > path{ AbstractPath( var , path_base() ) };
 const double lo = std::floor( value );
 const double up = std::ceil( value );
 std::vector< Change * > changes;
 switch( f_apply_strategy ) {
  case( eMaster ):
   changes.push_back( new LagrangianChange( AbstractChange::eChgUB ,
					    { lo } , path ) );
   changes.push_back( new LagrangianChange( AbstractChange::eChgLB ,
					    { up } , path ) );
   break;
  case( eSubproblem ):
   changes.push_back( new LagrangianChange( AbstractChange::eFixX ,
					    { lo } , path ) );
   changes.push_back( new LagrangianChange( AbstractChange::eFixX ,
					    { up } , path ) );
   break;
  default:
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::"
				 "branchings: unknown intApplyStrategy" ) );
  }

 return( changes );

 }  // end( LagrangianDualRelaxationSolver::branchings )

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::child_bound( Change * change )
{
 auto undo = apply( change , true );
 const int res = LagrangianDualSolver::compute( false );
 double bound;
 if( res == kInfeasible )
  bound = f_max ? - Inf< double >() : Inf< double >();
 else
  bound = f_max ? LagrangianDualSolver::get_ub()
                : LagrangianDualSolver::get_lb();
 apply( undo , false );
 delete undo;
 return( bound );

 }  // end( LagrangianDualRelaxationSolver::child_bound )

/*--------------------------------------------------------------------------*/

Change * LagrangianDualRelaxationSolver::apply( Change * change ,
						bool doUndo )
{
 auto c = dynamic_cast< LagrangianChange * >( change );
 if( ! c )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the Change is not a LagrangianChange" ) );

 const auto type = c->get_type();
 if( ( type != AbstractChange::eChgLB ) &&
     ( type != AbstractChange::eChgUB ) &&
     ( type != AbstractChange::eFixX ) &&
     ( type != AbstractChange::eUnfixX ) )
  return( c->apply( f_Block , doUndo ) );  // not a branching: to the Block

 if( c->get_paths().empty() )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the LagrangianChange has no path" ) );
 auto pv = c->get_paths()[ 0 ].get_element< ColVariable >( path_base() );
 if( ! pv )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the variable is not in the Block" ) );
 auto lbf = LagBF_of( pv );

 const auto undo = [ & ]( int utype , std::vector< double > data )
  -> Change * {
  if( ! doUndo )
   return( nullptr );
  return( new LagrangianChange( utype , std::move( data ) ,
				c->get_paths() ) );
  };

 // the bounds, in the Lagrangian Dual: a new dual pair of the LagBFunction,
 // or the constant term of the one already there - - - - - - - - - - - - -

 if( ( type == AbstractChange::eChgLB ) ||
     ( type == AbstractChange::eChgUB ) ) {
  if( f_apply_strategy != eMaster )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: a "
				 "bound Change needs intApplyStrategy == "
				 "eMaster" ) );
  if( c->get_data().empty() )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				 "the bound Change has no data" ) );
  const double value = c->get_data()[ 0 ];
  const bool lb = ( type == AbstractChange::eChgLB );

  // the pair of the bound: g( x ) = value - x <= 0 for a lower bound,
  // g( x ) = x - value <= 0 for an upper one
  auto & pair = map_varToLF[ pv ];
  auto & g = lb ? pair.second : pair.first;
  double old;
  if( g ) {
   old = lb ? g->get_constant_term() : - g->get_constant_term();
   g->set_constant_term( lb ? value : - value );
   }
  else {
   old = lb ? pv->get_lb() : pv->get_ub();
   g = new LinearFunction( LinearFunction::v_coeff_pair{
					   { pv , lb ? -1.0 : 1.0 } } ,
			   lb ? value : - value );
   auto y = new ColVariable();
   y->is_positive( true );
   lbf->add_dual_pairs( LagBFunction::v_dual_pair{ { y , g } } );
   }
  return( undo( type , { old } ) );
  }

 // the fixings, in the sub-Block - - - - - - - - - - - - - - - - - - - - -

 if( f_apply_strategy != eSubproblem )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: a "
				"fixing Change needs intApplyStrategy == "
				"eSubproblem" ) );
 if( ! f_global_information )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::apply: no "
			   "GlobalInformation" ) );
 if( ! map_varToSol )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::apply: no "
			   "Collection of the purged columns" ) );

 if( type == AbstractChange::eFixX ) {
  if( c->get_data().empty() )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				 "the fixing Change has no data" ) );
  Change * ret = pv->is_fixed() ? undo( AbstractChange::eFixX ,
					{ pv->get_value() } )
	                        : undo( AbstractChange::eUnfixX , {} );

  // the columns that the fixing purges from the global pool are kept, by
  // the variable, in the GlobalInformation: the handler is in place while
  // the fixing is done, as the LagBFunction purges them right then
  auto mvts = map_varToSol;
  const auto id = lbf->set_event_handler( LagBFunction::eColumnPurged ,
					  [ lbf , pv , mvts ]( void ) -> int {
   auto el = lbf->release_current_purged_solution();
   if( el.sol )
    mvts->write_with( str_PurgedColumns , [ & ]( PurgedColumn & pc ) {
     pc[ pv ].push_back( std::move( el ) );
     } );
   return( eContinue );
   } );
  pv->set_value( c->get_data()[ 0 ] );
  pv->is_fixed( true );
  lbf->reset_event_handler( LagBFunction::eColumnPurged , id );
  return( ret );
  }

 // eUnfixX: the columns purged by the fixing of the variable go back
 Change * ret = undo( AbstractChange::eFixX , { pv->get_value() } );
 pv->is_fixed( false );
 std::vector< LagBFunction::gpool_el > purged;
 map_varToSol->write_with( str_PurgedColumns , [ & ]( PurgedColumn & pc ) {
  if( auto it = pc.find( pv ) ; it != pc.end() ) {
   purged = std::move( it->second );
   pc.erase( it );
   }
  } );
 lbf->restore_purged_solutions( std::move( purged ) );
 return( ret );

 }  // end( LagrangianDualRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

LagBFunction * LagrangianDualRelaxationSolver::LagBF_of( Variable * var )
 const
{
 // the LagBFunction whose inner Block is the Block of the variable or one
 // of its ancestors
 for( Block * b = var->get_Block() ; b ; b = b->get_f_Block() )
  for( auto lbf : v_LBF )
   if( lbf && ( lbf->get_inner_block() == b ) )
    return( lbf );

 throw( std::invalid_argument( "LagrangianDualRelaxationSolver::LagBF_of: "
			       "the variable is in no sub-Block" ) );

 }  // end( LagrangianDualRelaxationSolver::LagBF_of )

/*--------------------------------------------------------------------------*/
/*------------- End File LagrangianDualRelaxationSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/
