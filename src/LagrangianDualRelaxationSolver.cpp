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

#include <limits>

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

namespace
{

// Run cleanup on every exit, including stack unwinding. The actions below
// restore an existing parameter, queue a shared snapshot, or undo a branch.
template< class Action >
class ScopeExit
{
 public:
 explicit ScopeExit( Action action ) : f_action( std::move( action ) ) { }
 ScopeExit( const ScopeExit & ) = delete;
 ScopeExit & operator=( const ScopeExit & ) = delete;
 ~ScopeExit() noexcept( noexcept( f_action() ) ) { f_action(); }

 private:
 Action f_action;
 };

// A restart is scoped to exactly one solve, including its error exits.
// PPH dispatches its inner iterations through compute_inner_solver(); the
// no-binary path and strong branching solve through LagrangianDualSolver.
class RestartParameter
{
 public:
 RestartParameter( CDASolver * solver , bool & pending ) : f_solver( solver ) {
  if( ! pending )
   return;
  const auto par = solver->int_par_str2idx( "intRstAlg" );
  if( par != Inf< ThinComputeInterface::idx_type >() ) {
   const int saved = solver->get_int_par( par );
   solver->set_par( par , saved | 2 | 4 );
   f_par = par;
   f_saved = saved;
   }
  pending = false;
  }
 RestartParameter( const RestartParameter & ) = delete;
 RestartParameter & operator=( const RestartParameter & ) = delete;
 ~RestartParameter() {
  if( f_par != Inf< ThinComputeInterface::idx_type >() )
   f_solver->set_par( f_par , f_saved );
  }

 private:
 CDASolver * f_solver;
 ThinComputeInterface::idx_type f_par =
   Inf< ThinComputeInterface::idx_type >();
 int f_saved = 0;
 };

// The cutoff of the inner Solver follows the incumbent of the cell during one
// solve, which the events eEverykIteration read every k iterations; an inner
// Solver that has no events is left as it is [see intCutOffReread].
class CutOffReread
{
 public:
 CutOffReread( CDASolver * solver , int k , std::atomic< double > * cell ,
	       bool max ) : f_solver( solver ) {
  if( ( k <= 0 ) || ( ! cell ) )
   return;
  const auto par = max ? Solver::dblLwCutOff : Solver::dblUpCutOff;
  try {
   f_id = solver->set_event_handler( ThinComputeInterface::eEverykIteration ,
				     [ solver , cell , max , par ]() {
    const double v = cell->load();
    if( std::isfinite( v ) &&
	( max ? ( v > solver->get_dbl_par( par ) )
	      : ( v < solver->get_dbl_par( par ) ) ) )
     solver->set_par( par , v );
    return( ThinComputeInterface::eContinue );
    } );
   }
  catch( std::logic_error & ) {
   return;
   }
  f_on = true;
  f_saved = solver->get_int_par( ThinComputeInterface::intEverykIt );
  solver->set_par( ThinComputeInterface::intEverykIt , k );
  }
 CutOffReread( const CutOffReread & ) = delete;
 CutOffReread & operator=( const CutOffReread & ) = delete;
 ~CutOffReread() {
  if( ! f_on )
   return;
  f_solver->set_par( ThinComputeInterface::intEverykIt , f_saved );
  f_solver->reset_event_handler( ThinComputeInterface::eEverykIteration ,
				 f_id );
  }

 private:
 CDASolver * f_solver;
 bool f_on = false;
 ThinComputeInterface::EventID f_id = 0;
 int f_saved = 0;
 };

}  // anonymous namespace

/*--------------------------------------------------------------------------*/
/*---------------- METHODS OF LagrangianDualRelaxationSolver ---------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_Block( Block * block )
{
 // The base class destroys the old LagBFunction and their terms. Forget
 // the non-owning pointers first; keep their multipliers alive until the
 // old dual and its inner Solver have been detached and destroyed.
 map_varToLF.clear();
 f_restart_pending = false;
 f_requested_center.reset();
 f_fresh_center.reset();
 f_pending.reset();
 f_estimates.clear();
 f_bounds.clear();
 f_pen_L = std::numeric_limits< double >::quiet_NaN();
 f_pc = PseudoCosts();
 f_n_branch = 0;
 f_ncompute = 0;
 f_solved = false;
 f_NoHeur = false;
 f_depth = 0;
 f_root_bound = std::numeric_limits< double >::quiet_NaN();
 f_heur_by_gap = false;
 PrimalProximalHeur::set_Block( block );
 f_bound_multipliers.clear();
 }

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
  case( intLagFixing ):
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intLagFixing is neither 0 nor 1" ) );
   f_lag_fixing = value;
   break;
  case( intHeurDepth ): f_heur_depth = value; break;
  case( intKeepBundle ):
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intKeepBundle is neither 0 nor 1" ) );
   f_keep_bundle = value;
   break;
  case( intCutOffReread ):
   if( value < 0 )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intCutOffReread is negative" ) );
   f_cutoff_reread = value;
   break;
  case( intReliability ):
   if( value < 0 )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intReliability is negative" ) );
   f_reliability = value;
   break;
  case( intBranchBound ):
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "LagrangianDualRelaxationSolver::set_par: "
				  "intBranchBound is neither 0 nor 1" ) );
   f_branch_bound = value;
   break;
  default: PrimalProximalHeur::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_par( idx_type par , double value )
{
 if( par == dblHeurGap )
  f_heur_gap = value;
 else
  PrimalProximalHeur::set_par( par , value );
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
  case( intLagFixing ):      return( 0 );
  case( intHeurDepth ):      return( -1 );
  case( intKeepBundle ):     return( 0 );
  case( intCutOffReread ):   return( 0 );
  case( intReliability ):    return( 0 );
  case( intBranchBound ):    return( 0 );
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
  case( intLagFixing ):      return( f_lag_fixing );
  case( intHeurDepth ):      return( f_heur_depth );
  case( intKeepBundle ):     return( f_keep_bundle );
  case( intCutOffReread ):   return( f_cutoff_reread );
  case( intReliability ):    return( f_reliability );
  case( intBranchBound ):    return( f_branch_bound );
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
 if( name == "intLagFixing" )
  return( intLagFixing );
 if( name == "intHeurDepth" )
  return( intHeurDepth );
 if( name == "intKeepBundle" )
  return( intKeepBundle );
 if( name == "intCutOffReread" )
  return( intCutOffReread );
 if( name == "intReliability" )
  return( intReliability );
 if( name == "intBranchBound" )
  return( intBranchBound );
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
 static const std::string lagfix = "intLagFixing";
 static const std::string hdepth = "intHeurDepth";
 static const std::string keep = "intKeepBundle";
 static const std::string reread = "intCutOffReread";
 static const std::string reliab = "intReliability";
 static const std::string bbound = "intBranchBound";
 switch( par ) {
  case( intApplyStrategy ):  return( apply );
  case( intBranchStrategy ): return( branch );
  case( intStrongCands ):    return( cands );
  case( intBranchEstimate ): return( estimate );
  case( intHeurFreq ):       return( heur );
  case( intLagFixing ):      return( lagfix );
  case( intHeurDepth ):      return( hdepth );
  case( intKeepBundle ):     return( keep );
  case( intCutOffReread ):   return( reread );
  case( intReliability ):    return( reliab );
  case( intBranchBound ):    return( bbound );
  default: return( PrimalProximalHeur::int_par_idx2str( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::get_dflt_dbl_par( idx_type par ) const
{
 if( par == dblHeurGap )
  return( -1 );
 return( PrimalProximalHeur::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::get_dbl_par( idx_type par ) const
{
 if( par == dblHeurGap )
  return( f_heur_gap );
 return( PrimalProximalHeur::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolver::dbl_par_str2idx(
					      const std::string & name ) const
{
 if( name == "dblHeurGap" )
  return( dblHeurGap );
 return( PrimalProximalHeur::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::dbl_par_idx2str(
						     idx_type par ) const
{
 static const std::string gap = "dblHeurGap";
 if( par == dblHeurGap )
  return( gap );
 return( PrimalProximalHeur::dbl_par_idx2str( par ) );
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
 f_scalars.reset();
 f_inc_cell = f_bound_cell = nullptr;
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

int LagrangianDualRelaxationSolver::compute_inner_solver( bool changedvars )
{
 prepare_center();
 RestartParameter restart( InnerSolver , f_restart_pending );
 // the incumbent is read again during the Lagrangian Dual of the original
 // objective only, the penalized ones having no cutoff
 const bool reread = f_cutoff_reread && ( ! f_Penalized );
 CutOffReread cutoff( InnerSolver , reread ? f_cutoff_reread : 0 ,
		      reread ? atomic_cell( GlobalInformation::str_Incumbent ,
					    f_inc_cell ) : nullptr , f_max );
 return( PrimalProximalHeur::compute_inner_solver( changedvars ) );
 }

/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolver::compute( bool changedvars )
{
 // the heuristic runs at the first call and at one every f_heur_freq, at
 // the nodes down to the depth f_heur_depth, and at those whose bound is
 // close to the best one, which after_Lagrangian_dual() tells [see
 // dblHeurGap]
 const bool by_count = ( f_ncompute == 0 ) ||
                       ( f_heur_freq && ( f_ncompute % f_heur_freq == 0 ) );
 const bool by_depth = ( f_heur_depth >= 0 ) &&
                       ( f_depth <= Index( f_heur_depth ) );
 f_NoHeur = ! ( by_count || by_depth );
 f_heur_by_gap = f_NoHeur && ( f_heur_gap >= 0 );
 if( f_heur_by_gap )
  f_NoHeur = false;
 ++f_ncompute;
 f_solved = false;
 f_lag_fix.clear();
 f_lag_pen.clear();
 f_pen_L = std::numeric_limits< double >::quiet_NaN();
 int status;
 if( NumStatVar == 0 ) {
  // PPH delegates this single solve to LagrangianDualSolver, bypassing the
  // hook used by its proximal iterations.
  prepare_center();
  RestartParameter restart( InnerSolver , f_restart_pending );
  CutOffReread cutoff( InnerSolver , f_cutoff_reread ,
		       f_cutoff_reread ?
		       atomic_cell( GlobalInformation::str_Incumbent ,
				    f_inc_cell ) : nullptr , f_max );
  status = PrimalProximalHeur::compute( changedvars );
  }
 else
  status = PrimalProximalHeur::compute( changedvars );

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
 f_estimates.clear();
 f_bounds.clear();
 if( f_solved )
  return( std::vector< Change * >() );
 auto cands = candidates();
 ++f_n_branch;
 const auto all = f_branch_estimate ? cands : std::vector< Cand >();
 const Index c = choose( cands );   // it may shorten cands
 auto changes = branchings( cands[ c ].var , cands[ c ].value );

 // the fixings by reduced cost go to both children, the branching variable
 // apart [see intLagFixing]
 if( f_lag_fixing && ( f_apply_strategy == eSubproblem ) ) {
  auto fix = std::move( f_lag_fix );
  f_lag_fix.clear();
  std::erase_if( fix , [ & ]( const auto & f ) {
   return( f.first == cands[ c ].var ); } );
  if( ! fix.empty() ) {
   f_n_lag_fixed += fix.size();
   for( auto & ch : changes ) {
    auto lc = static_cast< LagrangianChange * >( ch );
    auto paths = lc->get_paths();
    auto data = lc->get_data();
    for( const auto & f : fix ) {
     paths.emplace_back( f.first , f_Block );
     data.push_back( f.second );
     }
    auto nc = new LagrangianChange( lc->get_type() , std::move( data ) ,
				    std::move( paths ) );
    nc->set_center( lc->shared_center() );
    nc->set_depth( lc->get_depth() );
    delete lc;
    ch = nc;
    }
   }
  if( f_log && ( logVerb >= 1 ) )
   *f_log << "  branch " << f_n_branch << ": " << fix.size()
          << " fixed by reduced cost, " << f_n_lag_fixed << " so far"
          << std::endl;
  }
 if( learn_pseudocosts() ) {
  // the children record what their bound will tell of the pseudocosts,
  // the "down" one first [see branchings()]
  const Cand & cb = cands[ c ];
  static_cast< LagrangianChange * >( changes[ 0 ] )->set_branching(
		 { cb.k , 0 , cb.value - std::floor( cb.value ) , valid_bound } );
  static_cast< LagrangianChange * >( changes[ 1 ] )->set_branching(
		 { cb.k , 1 , std::ceil( cb.value ) - cb.value , valid_bound } );
  if( f_branch_estimate )
   f_estimates = estimates( all , cb.k );
  }

 // the bounds of the children by the penalties of fixing the variable,
 // the "down" one first [see intBranchBound]
 if( f_branch_bound && ( ! std::isnan( f_pen_L ) ) )
  if( const auto it = f_lag_pen.find( cands[ c ].var ) ;
      it != f_lag_pen.end() )
   for( int d = 0 ; d < 2 ; ++d )
    f_bounds.push_back( f_max ? f_pen_L - it->second[ d ]
		              : f_pen_L + it->second[ d ] );
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
  if( const auto it = f_lag_pen.find( cd.var ) ; it != f_lag_pen.end() ) {
   bc.back().down = it->second[ 0 ];
   bc.back().up = it->second[ 1 ];
   }
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
 const auto parent_center = f_requested_center ? f_requested_center : center();
 ScopeExit restore_parent(
   [ & ]() noexcept { f_requested_center = parent_center; } );

 // the children of a candidate: its two LagrangianChange, each applied,
 // solved and undone, within the time left to the node [see dblMaxTime]
 std::vector< std::array< bool , 2 > > tried( cands.size() ,
					      { false , false } );
 auto children = [ & ]( const BranchCandidate & bc ) {
  const auto & cd = cands[ bc.id ];
  auto changes = branchings( cd.var , cd.value );
  std::array< std::unique_ptr< Change > , 2 > owned{
   std::unique_ptr< Change >( changes[ 0 ] ) ,
   std::unique_ptr< Change >( changes[ 1 ] )
   };
  for( const auto & ch : owned )
   static_cast< LagrangianChange * >( ch.get() )->set_center( parent_center );
  std::array< double , 2 > bounds{ valid_bound , valid_bound };
  for( int d = 0 ; d < 2 ; ++d ) {
   // the trials are given what is left of dblMaxTime, a child not tried
   // counting as no improvement, and as no observation of the pseudocosts
   const double inner_time = InnerSolver->get_dbl_par( Solver::dblMaxTime );
   ScopeExit restore_time( [ & ]() {
    InnerSolver->set_par( Solver::dblMaxTime , inner_time ); } );
   if( ! give_time_left( InnerSolver ) )
    break;
   bounds[ d ] = child_bound( owned[ d ].get() );
   tried[ bc.id ][ d ] = true;
   }
  return( bounds );
  };

 auto evals = StrongBranching::evaluate( describe( cands ) , children ,
					 valid_bound , f_max ,
					 ( f_log && ( logVerb >= 1 ) ) ?
					 f_log : nullptr );

 // the bounds of the children are observations of the pseudocosts too
 if( learn_pseudocosts() )
  for( std::size_t i = 0 ; i < cands.size() ; ++i ) {
   const double v = cands[ i ].value;
   for( int d = 0 ; d < 2 ; ++d )
    if( tried[ i ][ d ] )
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
LagrangianDualRelaxationSolver::reliability_branching(
						 std::vector< Cand > & cands )
{
 // the product of the improvements of the two children, each at least
 // delta, a variable that intLagFixing fixes being no candidate
 std::unordered_set< ColVariable * > fixed;
 for( const auto & f : f_lag_fix )
  fixed.insert( f.first );
 const double delta = 1e-6 * std::max( std::abs( valid_bound ) , 1.0 );
 auto at_least = [ delta ]( double p ) {
  return( std::isnan( p ) ? delta : std::max( p , delta ) );
  };
 const auto eta = static_cast< unsigned long >( f_reliability );
 const double NaN = std::numeric_limits< double >::quiet_NaN();

 // a reliable candidate is scored by its pseudocosts, one that is not by
 // its penalties if known, or else (with eStrongBranching) by strong
 // branching, for at most intStrongCands of them, the most fractional
 // first, and by its pseudocosts if it gets neither
 std::vector< double > score( cands.size() , -1 );
 std::vector< Cand > strong;
 std::vector< Index > strong_c;
 for( Index c = 0 ; c < cands.size() ; ++c ) {
  const Cand & cd = cands[ c ];
  if( fixed.contains( cd.var ) )
   continue;
  if( ( f_pc.observations( cd.k , 0 ) < eta ) ||
      ( f_pc.observations( cd.k , 1 ) < eta ) ) {
   if( const auto it = f_lag_pen.find( cd.var ) ;
       ( it != f_lag_pen.end() ) && ( ! std::isnan( it->second[ 0 ] ) ) &&
       ( ! std::isnan( it->second[ 1 ] ) ) ) {
    score[ c ] = at_least( it->second[ 0 ] ) * at_least( it->second[ 1 ] );
    continue;
    }
   if( ( f_branch_strategy == eStrongBranching ) &&
       ( strong.size() < Index( std::max( f_strong_cands , 0 ) ) ) ) {
    strong.push_back( cd );
    strong_c.push_back( c );
    continue;
    }
   }
  score[ c ] = at_least( f_pc.get( cd.k , 0 , NaN ) *
			 ( cd.value - std::floor( cd.value ) ) ) *
               at_least( f_pc.get( cd.k , 1 , NaN ) *
			 ( std::ceil( cd.value ) - cd.value ) );
  }
 if( ! strong.empty() ) {
  const auto evals = strong_branching( strong );
  for( Index i = 0 ; i < strong.size() ; ++i )
   score[ strong_c[ i ] ] = evals[ i ][ 2 ];
  }

 // the best score, ties to the most fractional
 Index best_c = 0;
 for( Index c = 1 ; c < cands.size() ; ++c )
  if( score[ c ] > score[ best_c ] )
   best_c = c;
 return( best_c );

 }  // end( LagrangianDualRelaxationSolver::reliability_branching )

/*--------------------------------------------------------------------------*/

LagrangianDualRelaxationSolver::Index
LagrangianDualRelaxationSolver::choose( std::vector< Cand > & cands )
{
 if( f_branch_strategy == eMostFractional )
  return( 0 );

 if( f_reliability && ( ( f_branch_strategy == eStrongBranching ) ||
			( f_branch_strategy == ePenalties ) ) )
  return( reliability_branching( cands ) );

 if( f_branch_strategy == ePenalties ) {
  // the product of the penalties, each at least delta, a variable that
  // intLagFixing fixes being no candidate; ties to the most fractional
  std::unordered_set< ColVariable * > fixed;
  for( const auto & f : f_lag_fix )
   fixed.insert( f.first );
  const double delta = 1e-6 * std::max( std::abs( valid_bound ) , 1.0 );
  auto at_least = [ delta ]( double p ) {
   return( std::isnan( p ) ? delta : std::max( p , delta ) );
   };
  Index best_c = 0;
  double best_s = - 1;
  for( Index c = 0 ; c < cands.size() ; ++c ) {
   if( fixed.contains( cands[ c ].var ) )
    continue;
   const auto it = f_lag_pen.find( cands[ c ].var );
   const double s = it == f_lag_pen.end() ? delta * delta :
                    at_least( it->second[ 0 ] ) * at_least( it->second[ 1 ] );
   if( s > best_s ) {
    best_s = s;
    best_c = c;
    }
   }
  return( best_c );
  }

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
 std::vector< bool > hit( entries.size() , false );
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
       hit[ h ] = true;
       }
     ++pos;
     } );
   }

 // an entry that has no variable of the sub-Blocks, its group having none
 // in its range or no group having its name, branches upon nothing, which
 // is a mistake of the configuration rather than a node with nothing
 // fractional
 for( Index h = 0 ; h < entries.size() ; ++h )
  if( ! hit[ h ] )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::"
				 "branchable: entry " + f_branch_groups[ h ] +
				 " of vstrBranchGroups has no variable, " +
				 ( named[ h ] ? "the range being past the "
				                "end of the group"
				              : "no group of the sub-Blocks "
				                "being named " +
				                entries[ h ].name ) ) );

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

 // Children share the multipliers of their parent, also across best-first
 // jumps to another branch.
 const auto cntr = f_requested_center ? f_requested_center : center();
 for( auto ch : changes ) {
  static_cast< LagrangianChange * >( ch )->set_center( cntr );
  static_cast< LagrangianChange * >( ch )->set_depth( f_depth + 1 );
  }
 // After strong branching the parent centre is queued, not loaded yet.
 if( ! f_requested_center )
  f_fresh_center = cntr;

 return( changes );

 }  // end( LagrangianDualRelaxationSolver::branchings )

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::after_Lagrangian_dual( void )
{
 // the bound of the root, and whether the heuristic is run at this node
 // [see dblHeurGap]
 if( f_depth == 0 )
  f_root_bound = valid_bound;
 if( f_heur_by_gap ) {
  f_heur_by_gap = false;
  f_NoHeur = ! close_to_best_bound();
  }

 // the penalties are read once, for the fixings, for ePenalties, for the
 // bounds of the children [see intBranchBound] and for the features of the
 // candidates [see strBranchModel]
 const bool fixing = f_lag_fixing && ( f_apply_strategy == eSubproblem );
 const bool scoring = f_branch_bound ||
                      ( f_branch_strategy == ePenalties ) ||
                      ( f_branch_strategy == eLearned ) ||
                      ( f_branch_strategy == eOnline ) ||
                      ( ( f_branch_strategy == eStrongBranching ) &&
			( ( ! f_strong_log.empty() ) || f_reliability ) );
 if( ! ( fixing || scoring ) )
  return;

 penalties pen;
 lagrangian_penalties( pen );
 for( const auto & [ v , p0 , p1 ] : pen )
  f_lag_pen[ v ] = { p0 , p1 };
 if( fixing )
  lagrangian_fixings( f_lag_fix , & pen );
 if( f_branch_bound && ( ! pen.empty() ) ) {
  double L;
  if( lagrangian_value( L ) )
   f_pen_L = L;
  }
 }

/*--------------------------------------------------------------------------*/



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
  const LagrangianChange::Center & cntr )
{
 // Without a way of telling the inner Solver to read its Variable, merely
 // writing them would leave its internal centre unchanged.
 if( InnerSolver->int_par_str2idx( "intRstAlg" ) == Inf< idx_type >() )
  return( false );

 for( Index h = 0 ; h < v_LBF.size() ; ++h )
  if( auto lbf = v_LBF[ h ] ) {
   const auto values = h < cntr.size() ? &cntr[ h ] : nullptr;
   for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i )
    if( auto y = static_cast< ColVariable * >( lbf->get_active_var( i ) ) )
     y->set_value( ( values && ( i < values->size() ) ) ? ( *values )[ i ]
		                                            : 0 );
   }

 f_restart_pending = true;
 f_fresh_center.reset();
 return( true );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::prepare_center( void )
{
 // with intKeepBundle the inner Solver starts from where it is
 if( f_requested_center ) {
  if( ( ! f_keep_bundle ) && ( f_requested_center != f_fresh_center ) )
   restore_center( *f_requested_center );
  f_requested_center.reset();
  }
 // The solve about to start may move the multipliers, including on failure.
 f_fresh_center.reset();
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::child_bound( Change * change )
{
 std::unique_ptr< Change > undo( apply( change , true ) );
 ScopeExit restore_branch( [ & ] { apply( undo.get() , false ); } );
 prepare_center();
 RestartParameter restart( InnerSolver , f_restart_pending );
 const int res = LagrangianDualSolver::compute( false );
 return( res == kInfeasible ? ( f_max ? -Inf< double >() : Inf< double >() )
			     : ( f_max ? LagrangianDualSolver::get_ub()
				       : LagrangianDualSolver::get_lb() ) );

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
 if( learn_pseudocosts() )
  f_pending = c->get_branching();

 // Navigation only replaces a shared pointer. The last forward Change
 // supplies the selected node's parent centre; undo Changes carry none.
 // Loading it (and zeroing newer coordinates) waits until the next solve.
 if( c->get_center() )
  f_requested_center = c->shared_center();

 // the depth of the node the Change leads to, that of the node it leaves
 // being the one of its undo [see intHeurDepth]
 const Index depth = f_depth;
 if( c->get_depth() < Inf< Index >() )
  f_depth = c->get_depth();

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
  auto u = new LagrangianChange( utype , std::move( data ) ,
				 c->get_paths() );
  u->set_depth( depth );
  return( u );
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
   auto y = std::make_unique< ColVariable >();
   y->is_positive( true );
   f_bound_multipliers.push_back( std::move( y ) );
   lbf->add_dual_pairs(
     LagBFunction::v_dual_pair{ { f_bound_multipliers.back().get() , g } } );
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

 // one variable per path, the data holding the value each is fixed to; with
 // more than one, those already fixed are left as they are (and out of the
 // undo), as the fixings by reduced cost [see intLagFixing] only add to the
 // branching [see branch()]
 const auto & paths = c->get_paths();
 const bool many = ( paths.size() > 1 );
 if( ( type == AbstractChange::eFixX ) &&
     ( c->get_data().size() < paths.size() ) )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the fixing Change has fewer data than "
				"paths" ) );

 std::vector< AbstractPath > upaths;
 std::vector< double > udata;
 int undo_type = AbstractChange::eUnfixX;

 for( Index i = 0 ; i < paths.size() ; ++i ) {
  auto v = i ? paths[ i ].get_element< ColVariable >( f_Block ) : pv;
  if( ! v )
   throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				 "the variable is not in the Block" ) );
  auto vlbf = i ? LagBF_of( v ) : lbf;

  if( type == AbstractChange::eFixX ) {
   if( v->is_fixed() ) {
    if( many )
     continue;
    undo_type = AbstractChange::eFixX;  // the undo gives back its value
    udata.push_back( v->get_value() );
    }
   upaths.push_back( paths[ i ] );

   // the columns that the fixing purges from the global pool are kept, by
   // the variable, in the GlobalInformation: the handler is in place while
   // the fixing is done, as the LagBFunction purges them right then
   auto mvts = map_varToSol;
   const auto id = vlbf->set_event_handler( LagBFunction::eColumnPurged ,
					    [ vlbf , v , mvts ]( void ) -> int {
    auto el = vlbf->release_current_purged_solution();
    if( el.sol )
     mvts->write_with( str_PurgedColumns , [ & ]( PurgedColumn & pc ) {
      pc[ v ].push_back( std::move( el ) );
      } );
    return( eContinue );
    } );
   v->set_value( c->get_data()[ i ] );
   v->is_fixed( true );
   vlbf->reset_event_handler( LagBFunction::eColumnPurged , id );
   continue;
   }

  // eUnfixX: the columns purged by the fixing of the variable go back
  if( many && ( ! v->is_fixed() ) )
   continue;
  undo_type = AbstractChange::eFixX;
  upaths.push_back( paths[ i ] );
  udata.push_back( v->get_value() );
  v->is_fixed( false );
  std::vector< LagBFunction::gpool_el > purged;
  map_varToSol->write_with( str_PurgedColumns , [ & ]( PurgedColumn & pc ) {
   if( auto it = pc.find( v ) ; it != pc.end() ) {
    purged = std::move( it->second );
    pc.erase( it );
    }
   } );
  vlbf->restore_purged_solutions( std::move( purged ) );
  }

 if( ( ! doUndo ) || upaths.empty() )
  return( nullptr );
 if( undo_type == AbstractChange::eUnfixX )
  udata.clear();
 auto u = new LagrangianChange( undo_type , std::move( udata ) ,
				std::move( upaths ) );
 u->set_depth( depth );
 return( u );

 }  // end( LagrangianDualRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

std::atomic< double > * LagrangianDualRelaxationSolver::atomic_cell(
			      const char * key , std::atomic< double > * & cell )
{
 if( ( ! cell ) && f_global_information ) {
  if( ! f_scalars )
   f_scalars = f_global_information->get_from_Universe<
     std::atomic< double > >( GlobalInformation::str_AtomicScalars );
  if( f_scalars && f_scalars->contains( key ) )
   cell = &( ( *f_scalars )[ key ] );
  }
 return( cell );

 }  // end( LagrangianDualRelaxationSolver::atomic_cell )

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolver::global_incumbent( void )
{
 double zbar = incumbent();
 if( auto cell = atomic_cell( GlobalInformation::str_Incumbent ,
			      f_inc_cell ) )
  if( const double v = cell->load() ; std::isfinite( v ) )
   zbar = f_max ? std::max( zbar , v ) : std::min( zbar , v );
 return( zbar );

 }  // end( LagrangianDualRelaxationSolver::global_incumbent )

/*--------------------------------------------------------------------------*/

bool LagrangianDualRelaxationSolver::close_to_best_bound( void )
{
 // with no incumbent, or no bound to measure the node against, the
 // heuristic is run
 const double zbar = global_incumbent();
 double best = f_root_bound;
 if( auto cell = atomic_cell( GlobalInformation::str_BestBound ,
			      f_bound_cell ) )
  if( const double b = cell->load() ; std::isfinite( b ) )
   best = b;
 if( ( ! std::isfinite( zbar ) ) || ( ! std::isfinite( best ) ) ||
     ( ! std::isfinite( valid_bound ) ) )
  return( true );

 const double gap = ( f_max ? best - valid_bound : valid_bound - best ) /
                    std::max( std::abs( zbar ) , 1.0 );
 return( gap <= f_heur_gap );

 }  // end( LagrangianDualRelaxationSolver::close_to_best_bound )

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
