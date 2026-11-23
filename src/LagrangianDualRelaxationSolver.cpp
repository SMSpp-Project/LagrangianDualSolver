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

#include <cctype>

#include <unordered_map>

#include <unordered_set>

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
  case( intBranchEstimate ):
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intBranchEstimate is neither 0 nor 1" ) );
   f_branch_estimate = value;
   break;
  case( intHeurFreq ):
   if( value < 0 )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intHeurFreq is negative" ) );
   f_heur_freq = value;
   break;
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
  case( intBranchEstimate ): return( 0 );
  case( intHeurFreq ):       return( 1 );
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
  case( intBranchEstimate ): return( f_branch_estimate );
  case( intHeurFreq ):       return( f_heur_freq );
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
 if( name == "intBranchEstimate" )
  return( intBranchEstimate );
 if( name == "intHeurFreq" )
  return( intHeurFreq );
 return( PrimalProximalHeur::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::int_par_idx2str(
						       idx_type par ) const
{
 static const std::string apply = "intApplyStrategy";
 static const std::string branch = "intBranchStrategy";
 static const std::string cands = "intStrongCands";
 static const std::string estimate = "intBranchEstimate";
 static const std::string heur = "intHeurFreq";
 switch( par ) {
  case( intApplyStrategy ):  return( apply );
  case( intBranchStrategy ): return( branch );
  case( intStrongCands ):    return( cands );
  case( intBranchEstimate ): return( estimate );
  case( intHeurFreq ):       return( heur );
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
  if( par == strBranchModel ) {
   f_branch_model = std::move( value );
   f_model = LearnedBranching();  // read again at the next use
   }
 else
  if( par == strEstimateModel ) {
   f_estimate_model = std::move( value );
   f_pc_model = LearnedPseudoCosts();  // read again at the next use
   }
 else
  PrimalProximalHeur::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::get_dflt_str_par(
						     idx_type par ) const
{
 static const std::string empty;
 if( ( par == strStrongLog ) || ( par == strBranchModel ) ||
     ( par == strEstimateModel ) )
  return( empty );
 return( PrimalProximalHeur::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::get_str_par(
						     idx_type par ) const
{
 if( par == strStrongLog )
  return( f_strong_log );
 if( par == strBranchModel )
  return( f_branch_model );
 if( par == strEstimateModel )
  return( f_estimate_model );
 return( PrimalProximalHeur::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolver::str_par_str2idx(
					      const std::string & name ) const
{
 if( name == "strStrongLog" )
  return( strStrongLog );
 if( name == "strBranchModel" )
  return( strBranchModel );
 if( name == "strEstimateModel" )
  return( strEstimateModel );
 return( PrimalProximalHeur::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::str_par_idx2str(
						     idx_type par ) const
{
 static const std::string log = "strStrongLog";
 static const std::string model = "strBranchModel";
 static const std::string emodel = "strEstimateModel";
 if( par == strStrongLog )
  return( log );
 if( par == strBranchModel )
  return( model );
 if( par == strEstimateModel )
  return( emodel );
 return( PrimalProximalHeur::str_par_idx2str( par ) );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_par( idx_type par ,
				      std::vector< std::string > && value )
{
 if( par == vstrBranchGroups )
  f_branch_groups = std::move( value );
 else  // PrimalProximalHeur has none of its own
  LagrangianDualSolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
LagrangianDualRelaxationSolver::get_dflt_vstr_par( idx_type par ) const
{
 static const std::vector< std::string > empty;
 if( par == vstrBranchGroups )
  return( empty );
 return( PrimalProximalHeur::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
LagrangianDualRelaxationSolver::get_vstr_par( idx_type par ) const
{
 if( par == vstrBranchGroups )
  return( f_branch_groups );
 return( PrimalProximalHeur::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolver::vstr_par_str2idx(
					      const std::string & name ) const
{
 if( name == "vstrBranchGroups" )
  return( vstrBranchGroups );
 return( PrimalProximalHeur::vstr_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::vstr_par_idx2str(
						     idx_type par ) const
{
 static const std::string groups = "vstrBranchGroups";
 if( par == vstrBranchGroups )
  return( groups );
 return( PrimalProximalHeur::vstr_par_idx2str( par ) );
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
 // the heuristic runs at the first call and at one every f_heur_freq
 f_NoHeur = ( f_ncompute > 0 ) &&
            ( ( f_heur_freq == 0 ) || ( f_ncompute % f_heur_freq != 0 ) );
 ++f_ncompute;
 f_solved = false;
 const int status = PrimalProximalHeur::compute( changedvars );

 // the bound of a child just after its parent branched is an observation
 // of the pseudocost of the variable branched upon [see intBranchEstimate];
 // an infeasible child (an infinite bound) is none
 if( f_pending ) {
  const auto & b = *f_pending;
  f_pc.record( b.k , b.dir , b.moved , f_max ? b.node - valid_bound
		                             : valid_bound - b.node );
  f_pending.reset();
  }

 // beyond the cutoff [see dblUpCutOff]: the node is fenced, which is what
 // whoever set the cutoff needs to know
 if( status == kCutOff )
  return( status );

 // the restart that restore_center() asked for is over
 if( f_rst_alg >= 0 ) {
  InnerSolver->set_par( InnerSolver->int_par_str2idx( "intRstAlg" ) ,
			f_rst_alg );
  f_rst_alg = -1;
  }

 // the multipliers are now those of this node [see apply()]
 f_fresh_center = nullptr;

 const bool done = ( status >= kOK ) &&
                   ( ( status == kLowPrecision ) || ( status < kError ) );
 if( ! done )
  return( status );

 // without the convexified solution (say, the master problem of the inner
 // Solver failed, leaving no linearization to combine) the node cannot be
 // branched: what the inner Solver returned is said, or kError if it said
 // kOK
 const int ld_status = get_Lagrangian_initial_status();
 const auto & sol = get_Lagrangian_convexified_solution();
 if( sol.size() != NumStatVar )
  return( ld_status == kOK ? int( kError ) : ld_status );

 // the node is branched on a variable fractional in the convexified
 // solution, among those branched upon, as candidates() does
 const auto rank = branchable();
 bool frac = false;
 for( Index k = 0 ; ( ! frac ) && ( k < sol.size() ) ; ++k )
  frac = ( rank[ k ] < Inf< Index >() ) &&
         ( std::abs( sol[ k ] - std::round( sol[ k ] ) ) > FracEps );
 if( frac )
  return( kOK );

 // with none, the node is solved if the Lagrangian Dual is and a feasible
 // solution is within the accuracy of its bound (the convexified solution
 // itself, if integral and feasible [see PrimalProximalHeur::compute()]),
 // and branch() gives no Change; otherwise it can be neither branched nor
 // fenced, and the status is that of the Lagrangian Dual (stopped by a
 // budget, say), or kLowPrecision if it is solved but a variable not
 // branched upon is fractional, or the integral solution is not feasible
 if( ld_status != kOK )
  return( ld_status );
 if( ! gap_closed() )
  return( kLowPrecision );
 f_solved = true;
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
 if( ! solc ) {
  if( v_best_sol.empty() )
   return( nullptr );
  return( v_best_sol.front().first->clone() );  // the caller owns it
  }

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
 if( f_solved )
  return( std::vector< Change * >() );
 auto cands = candidates();
 ++f_n_branch;
 f_estimates.clear();
 const auto all = f_branch_estimate ? cands : std::vector< Cand >();
 const Index c = choose( cands );   // it may shorten cands
 auto changes = branchings( cands[ c ].var , cands[ c ].value );
 if( f_branch_estimate ) {
  // the children record what their bound will tell of the pseudocosts,
  // the "down" one first [see branchings()]
  const Cand & cb = cands[ c ];
  static_cast< LagrangianChange * >( changes[ 0 ] )->set_branching(
		 { cb.k , 0 , cb.value - std::floor( cb.value ) , valid_bound } );
  static_cast< LagrangianChange * >( changes[ 1 ] )->set_branching(
		 { cb.k , 1 , std::ceil( cb.value ) - cb.value , valid_bound } );
  f_estimates = estimates( all , cb.k );
  }
 return( changes );

 }  // end( LagrangianDualRelaxationSolver::branch )

/*--------------------------------------------------------------------------*/

std::vector< double > LagrangianDualRelaxationSolver::estimates(
				    const std::vector< Cand > & cands , Index k )
{
 auto bc = describe( cands );
 std::size_t b = bc.size();
 for( std::size_t i = 0 ; i < bc.size() ; ++i ) {
  bc[ i ].id = cands[ i ].k;    // the same variable at all the nodes
  if( cands[ i ].k == k )
   b = i;
  }
 if( b == bc.size() )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::estimates: the "
			   "variable branched upon is not a candidate" ) );

 PseudoCosts::Prior prior;
 std::vector< Features > x;
 if( ! f_estimate_model.empty() ) {
  if( f_pc_model.empty() )
   f_pc_model.load( f_estimate_model );
  x = features( cands );
  prior = [ & ]( std::size_t i ) {
   return( f_pc_model.prior( x[ i ] , valid_bound ) );
   };
  }
 const auto e = f_pc.estimates( bc , b , valid_bound , f_max , prior );
 return( std::vector< double >{ e[ 0 ] , e[ 1 ] } );

 }  // end( LagrangianDualRelaxationSolver::estimates )

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::fixed_fraction( void ) const
{
 Index n_fixed = 0;
 for( const auto & sbd : idx_to_var_sbi1 )
  for( const auto & dv : sbd )
   if( dv.second->is_fixed() )
    ++n_fixed;
 return( NumStatVar ? double( n_fixed ) / NumStatVar : 0.0 );
 }

/*--------------------------------------------------------------------------*/

std::vector< BranchCandidate >
LagrangianDualRelaxationSolver::describe( const std::vector< Cand > & cands )
 const
{
 const auto & lagr = get_Lagrangian_initial_solution();
 std::vector< BranchCandidate > bc;
 bc.reserve( cands.size() );
 for( const auto & cd : cands ) {
  const auto n_sb = idx_to_var_sbi1[ cd.sb ].size();
  bc.push_back( { std::size_t( &cd - cands.data() ) , cd.value ,
		  cd.k < lagr.size() ? lagr[ cd.k ] : cd.value , cd.cost ,
		  n_sb > 1 ? double( cd.pos ) / ( n_sb - 1 ) : 0.0 ,
		  static_cast< unsigned int >( cd.sb ) } );
  }
 return( bc );

 }  // end( LagrangianDualRelaxationSolver::describe )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::Cand >
LagrangianDualRelaxationSolver::candidates( void )
{
 // the convexified Lagrangian solution, one value per binary variable of
 // the sub-Blocks: that of the sub-Blocks alone is integral whenever they
 // are, and it says nothing about where the relaxation is fractional
 const auto & sol = get_Lagrangian_convexified_solution();
 if( sol.empty() )
  throw( std::runtime_error( "LagrangianDualRelaxationSolver::candidates: "
			     "no Lagrangian solution available" ) );
 if( sol.size() != NumStatVar )
  throw( std::runtime_error( "LagrangianDualRelaxationSolver::candidates: "
			     "the Lagrangian solution has the wrong size" ) );

 // the candidates: the fractional variables among those that may be
 // branched upon [see vstrBranchGroups], the most fractional first
 std::vector< Cand > cands;
 {
  const auto rank = branchable();
  Index k = 0;
  for( Index sb = 0 ; sb < idx_to_var_sbi1.size() ; ++sb )
   for( Index pos = 0 ; pos < idx_to_var_sbi1[ sb ].size() ; ++pos , ++k ) {
    const auto & dv = idx_to_var_sbi1[ sb ][ pos ];
    const double v = sol[ k ];
    const double frac = std::abs( v - std::round( v ) );
    if( ( rank[ k ] < Inf< Index >() ) && ( frac > FracEps ) )
     cands.push_back( { dv.second , v , frac , k , sb , pos , dv.first ,
			rank[ k ] } );
    }
  }

 // no fractional variable: compute() says kOK for such a node only if it
 // has solved it, and then branch() asks for no candidate, so being asked
 // for them is an error of the caller, which has branched a node whose
 // compute() said it cannot be [see compute()]
 if( cands.empty() )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::candidates: no "
			   "variable branched upon is fractional, and the "
			   "node cannot be branched" ) );

 // those of the first entry of vstrBranchGroups that has any, the most
 // fractional first [see vstrBranchGroups]
 std::stable_sort( cands.begin() , cands.end() ,
		   []( const Cand & a , const Cand & b ) {
		    return( a.rank != b.rank ? a.rank < b.rank
			                     : a.frac > b.frac ); } );
 const auto first = cands.front().rank;
 cands.erase( std::find_if( cands.begin() , cands.end() ,
			    [ first ]( const Cand & c ) {
			     return( c.rank != first ); } ) , cands.end() );
 return( cands );

 }  // end( LagrangianDualRelaxationSolver::candidates )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::Features >
LagrangianDualRelaxationSolver::features( const std::vector< Cand > & cands )
 const
{
 return( BranchingFeatures::features( describe( cands ) , fixed_fraction() ) );

 }  // end( LagrangianDualRelaxationSolver::features )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::StrongEval >
LagrangianDualRelaxationSolver::strong_branching(
					       const std::vector< Cand > & cands )
{
 // the children of a candidate: its two LagrangianChange, each applied,
 // solved and undone
 auto children = [ & ]( const BranchCandidate & bc ) {
  const auto & cd = cands[ bc.id ];
  auto changes = branchings( cd.var , cd.value );
  std::array< double , 2 > bounds;
  try {
   bounds = { child_bound( changes[ 0 ] ) , child_bound( changes[ 1 ] ) };
   }
  catch( ... ) {
   for( auto ch : changes )
    delete ch;
   throw;
   }
  for( auto ch : changes )
   delete ch;
  return( bounds );
  };

 auto evals = StrongBranching::evaluate( describe( cands ) , children ,
					 valid_bound , f_max ,
					 ( f_log && ( logVerb >= 1 ) ) ?
					 f_log : nullptr );

 // the bounds of the children are observations of the pseudocosts too
 if( f_branch_estimate )
  for( std::size_t i = 0 ; i < cands.size() ; ++i ) {
   const double v = cands[ i ].value;
   for( int d = 0 ; d < 2 ; ++d )
    f_pc.record( cands[ i ].k , d , d ? std::ceil( v ) - v
		                       : v - std::floor( v ) ,
		 f_max ? valid_bound - evals[ i ][ d ]
		       : evals[ i ][ d ] - valid_bound );
   }
 return( evals );

 }  // end( LagrangianDualRelaxationSolver::strong_branching )

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::write_strong_log(
					   const std::vector< Cand > & cands ,
					   const std::vector< StrongEval > & evals ,
					   Index chosen )
{
 if( f_strong_log.empty() )
  return;

 StrongBranching::write_log( f_strong_log , f_n_branch , valid_bound ,
			     describe( cands ) , features( cands ) , evals ,
			     chosen );

 }  // end( LagrangianDualRelaxationSolver::write_strong_log )

/*--------------------------------------------------------------------------*/

LagrangianDualRelaxationSolver::Index
LagrangianDualRelaxationSolver::choose( std::vector< Cand > & cands )
{
 if( f_branch_strategy == eMostFractional )
  return( 0 );

 if( ( f_branch_strategy != eStrongBranching ) &&
     ( f_branch_strategy != eLearned ) ) {
  if( f_branch_strategy == eOnline )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::choose: "
				 "intBranchStrategy eOnline needs "
				 "LagrangianDualRelaxationSolverML" ) );
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::choose: "
				"unknown intBranchStrategy" ) );
  }

 if( f_strong_cands <= 0 )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::choose: "
				"intStrongCands must be positive" ) );
 if( cands.size() > Index( f_strong_cands ) )
  cands.resize( f_strong_cands );

 if( f_branch_strategy == eStrongBranching ) {
  const auto evals = strong_branching( cands );
  const Index best_c = StrongBranching::best( evals );
  write_strong_log( cands , evals , best_c );
  return( best_c );
  }

 // eLearned: the candidate the model of strBranchModel scores best
 if( f_model.empty() ) {
  if( f_branch_model.empty() )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::choose: "
				 "eLearned needs strBranchModel" ) );
  f_model.load( f_branch_model );
  }
 return( f_model.best( features( cands ) ) );

 }  // end( LagrangianDualRelaxationSolver::choose )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::Index >
LagrangianDualRelaxationSolver::branchable( void ) const
{
 std::vector< Index > rank( NumStatVar ,
			    f_branch_groups.empty() ? 0 : Inf< Index >() );
 if( f_branch_groups.empty() )
  return( rank );

 // each entry is "name", "name[ i ]" or "name[ i : j ]", blanks being
 // allowed anywhere [see vstrBranchGroups]
 struct Entry {
  std::string name;
  Index lo , hi;
  };
 std::vector< Entry > entries;
 for( const auto & e : f_branch_groups ) {
  std::string s;
  for( auto c : e )
   if( ! std::isspace( static_cast< unsigned char >( c ) ) )
    s.push_back( c );
  const auto br = s.find( '[' );
  if( br == std::string::npos ) {
   entries.push_back( { s , 0 , Inf< Index >() } );
   continue;
   }
  if( s.back() != ']' )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::"
				 "branchable: invalid entry " + e +
				 " of vstrBranchGroups" ) );
  const auto in = s.substr( br + 1 , s.size() - br - 2 );
  const auto col = in.find( ':' );
  const Index lo = std::stoul( in.substr( 0 , col ) );
  const Index hi = ( col == std::string::npos ) ? lo + 1
                   : Index( std::stoul( in.substr( col + 1 ) ) );
  entries.push_back( { s.substr( 0 , br ) , lo , hi } );
  }

 // the first entry that has each variable of the groups of the sub-Blocks,
 // then the ranks in the order in which PrimalProximalHeur has registered
 // the variables
 std::unordered_map< const ColVariable * , Index > first;
 std::vector< bool > named( entries.size() , false );
 for( const auto & sbi : f_Block->get_nested_Blocks() )
  for( const auto & group : sbi->get_static_variable_groups() ) {
   if( ! group )
    continue;
   Index pos = 0;
   for( Index h = 0 ; h < entries.size() ; ++h )
    named[ h ] = named[ h ] || ( entries[ h ].name == group->get_name() );
   group->for_each_as< ColVariable >( [ & ]( ColVariable & var ) {
     for( Index h = 0 ; h < entries.size() ; ++h )
      if( ( entries[ h ].name == group->get_name() ) &&
	  ( pos >= entries[ h ].lo ) && ( pos < entries[ h ].hi ) ) {
       first.try_emplace( & var , h );
       break;
       }
     ++pos;
     } );
   }

 // an entry that names no group branches upon nothing, which is a mistake
 // of the configuration rather than a node with nothing fractional
 for( Index h = 0 ; h < entries.size() ; ++h )
  if( ! named[ h ] )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::"
				 "branchable: entry " + f_branch_groups[ h ] +
				 " of vstrBranchGroups names no group of the "
				 "sub-Blocks" ) );

 Index k = 0;
 for( const auto & sbv : idx_to_var_sbi1 )
  for( const auto & dv : sbv ) {
   if( const auto it = first.find( dv.second ) ; it != first.end() )
    rank[ k ] = it->second;
   ++k;
   }

 return( rank );

 }  // end( LagrangianDualRelaxationSolver::branchable )

/*--------------------------------------------------------------------------*/

std::vector< Change * > LagrangianDualRelaxationSolver::branchings(
				      ColVariable * var , double value ) const
{
 // the variable is one of the sub-Blocks of f_Block [see the class notes],
 // while the Lagrangian Dual solves their copies
 if( iBCopy )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::branchings: "
				"the branching on a variable needs the "
				"sub-Blocks not to be copied, see "
				"int_LDSlv_iBCopy" ) );
 const std::vector< AbstractPath > path{ AbstractPath( var , f_Block ) };
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

 // the children start from the multipliers of this node [see apply()]
 const auto cntr = center();
 for( auto ch : changes )
  static_cast< LagrangianChange * >( ch )->set_center( cntr );
 f_fresh_center = cntr.get();

 return( changes );

 }  // end( LagrangianDualRelaxationSolver::branchings )

/*--------------------------------------------------------------------------*/

std::shared_ptr< const LagrangianChange::Center >
                       LagrangianDualRelaxationSolver::center( void ) const
{
 auto cntr = std::make_shared< LagrangianChange::Center >( v_LBF.size() );
 for( Index h = 0 ; h < v_LBF.size() ; ++h )
  if( auto lbf = v_LBF[ h ] ) {
   auto & values = ( *cntr )[ h ];
   values.resize( lbf->get_num_active_var() );
   for( Index i = 0 ; i < values.size() ; ++i )
    if( auto y = static_cast< ColVariable * >( lbf->get_active_var( i ) ) )
     values[ i ] = y->get_value();
   }
 return( cntr );

 }  // end( LagrangianDualRelaxationSolver::center )

/*--------------------------------------------------------------------------*/

bool LagrangianDualRelaxationSolver::restore_center(
				     const LagrangianChange::Center & center )
{
 const auto par = InnerSolver->int_par_str2idx( "intRstAlg" );
 if( par == Inf< idx_type >() )
  return( false );

 for( Index h = 0 ; h < v_LBF.size() ; ++h )
  if( auto lbf = v_LBF[ h ] ) {
   const auto values = h < center.size() ? &center[ h ] : nullptr;
   for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i )
    if( auto y = static_cast< ColVariable * >( lbf->get_active_var( i ) ) )
     y->set_value( ( values && ( i < values->size() ) ) ? ( *values )[ i ]
		                                            : 0 );
   }

 // the bundle is emptied and the algorithm restarts from the Variable, at
 // the next compute() only
 if( f_rst_alg < 0 ) {
  f_rst_alg = InnerSolver->get_int_par( par );
  InnerSolver->set_par( par , f_rst_alg | 2 | 4 );
  }
 return( true );

 }  // end( LagrangianDualRelaxationSolver::restore_center )

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::child_bound( Change * change )
{
 auto undo = apply( change , true );
 const int res = LagrangianDualSolver::compute( false );
 f_fresh_center = nullptr;  // the multipliers are those of the child now
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

 // the next compute() solves the child of this branching, if any, as the
 // Changes of the path to a node are applied from the root down
 f_pending.reset();
 if( f_branch_estimate )
  f_pending = c->get_branching();

 // a child starts from the multipliers its father had, not from those of
 // whatever node was solved last, which may be anywhere in the tree and may
 // have been infeasible, its multipliers grown without bound: going down to
 // it, each LagrangianChange writes those of its own father, so that the
 // last one, the father of the node, is what is left. The first child
 // solved right after its father branched finds them there already
 if( c->get_center() && ( c->get_center() != f_fresh_center ) )
  restore_center( *c->get_center() );

 const auto type = c->get_type();
 if( ( type != AbstractChange::eChgLB ) &&
     ( type != AbstractChange::eChgUB ) &&
     ( type != AbstractChange::eFixX ) &&
     ( type != AbstractChange::eUnfixX ) )
  return( c->apply( f_Block , doUndo ) );  // not a branching: to the Block

 if( c->get_paths().empty() )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the LagrangianChange has no path" ) );
 auto pv = c->get_paths()[ 0 ].get_element< ColVariable >( f_Block );
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
  // g( x ) = x - value <= 0 for an upper one, the opposite for a
  // maximization problem (as the relaxed constraints are reversed [see
  // to_be_reversed()])
  const bool rev = ( lb != f_max );
  auto & pair = map_varToLF[ pv ];
  auto & g = lb ? pair.second : pair.first;
  double old;
  if( g ) {
   old = rev ? g->get_constant_term() : - g->get_constant_term();
   g->set_constant_term( rev ? value : - value );
   }
  else {
   old = lb ? pv->get_lb() : pv->get_ub();
   g = new LinearFunction( LinearFunction::v_coeff_pair{
					   { pv , rev ? -1.0 : 1.0 } } ,
			   rev ? value : - value );
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
