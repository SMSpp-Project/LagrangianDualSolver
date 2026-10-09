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
  if( par == strBranchModel ) {
   f_branch_model = std::move( value );
   f_w1.clear();  // read again at the next use
   }
 else
  PrimalProximalHeur::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::get_dflt_str_par(
						     idx_type par ) const
{
 static const std::string empty;
 if( ( par == strStrongLog ) || ( par == strBranchModel ) )
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
 return( PrimalProximalHeur::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolver::str_par_idx2str(
						     idx_type par ) const
{
 static const std::string log = "strStrongLog";
 static const std::string model = "strBranchModel";
 if( par == strStrongLog )
  return( log );
 if( par == strBranchModel )
  return( model );
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
 const int status = PrimalProximalHeur::compute( changedvars );

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

 // the Lagrangian Dual not solved (stopped by a budget, or kLowPrecision)
 // gives a valid bound, and a convexified solution that need not satisfy
 // the relaxed constraints: the node can still be branched on a variable
 // fractional in it, but if there is none it can be neither branched nor
 // fenced, and what the inner Solver returned is said
 if( ld_status != kOK ) {
  const auto mask = branchable();
  bool frac = false;
  for( Index k = 0 ; ( ! frac ) && ( k < sol.size() ) ; ++k )
   frac = mask[ k ] && ( std::abs( sol[ k ] - std::round( sol[ k ] ) )
			 > FracEps );
  if( ! frac )
   return( ld_status );
  }

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
 auto cands = candidates();
 ++f_n_branch;
 const Index c = choose( cands );
 return( branchings( cands[ c ].var , cands[ c ].value ) );

 }  // end( LagrangianDualRelaxationSolver::branch )

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
  const auto mask = branchable();
  Index k = 0;
  for( Index sb = 0 ; sb < idx_to_var_sbi1.size() ; ++sb )
   for( Index pos = 0 ; pos < idx_to_var_sbi1[ sb ].size() ; ++pos , ++k ) {
    const auto & dv = idx_to_var_sbi1[ sb ][ pos ];
    const double v = sol[ k ];
    const double frac = std::abs( v - std::round( v ) );
    if( mask[ k ] && ( frac > FracEps ) )
     cands.push_back( { dv.second , v , frac , k , sb , pos , dv.first } );
    }
  }

 // no fractional variable: an integral Lagrangian solution of a
 // well-terminated Lagrangian Dual satisfies the relaxed constraints, i.e.,
 // it solves the node, which is then fenced by bound and never branched;
 // were it to be branched, the two children would be the same node, hence
 // this is an error, of the termination of the inner Solver (say, dblNZEps
 // too large), or of vstrBranchGroups leaving out a group that the others
 // do not determine
 if( cands.empty() )
  throw( std::logic_error( "LagrangianDualRelaxationSolver::candidates: no "
			   "variable branched upon is fractional, either the "
			   "relaxed constraints are violated beyond the "
			   "termination tolerance of the Lagrangian Dual "
			   "(dblNZEps) or vstrBranchGroups leaves out a "
			   "fractional group" ) );

 std::stable_sort( cands.begin() , cands.end() ,
		   []( const Cand & a , const Cand & b ) {
		    return( a.frac > b.frac ); } );
 return( cands );

 }  // end( LagrangianDualRelaxationSolver::candidates )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::Features >
LagrangianDualRelaxationSolver::features( const std::vector< Cand > & cands )
 const
{
 Index n_fixed = 0;
 for( const auto & sbd : idx_to_var_sbi1 )
  for( const auto & dv : sbd )
   if( dv.second->is_fixed() )
    ++n_fixed;
 const double fixed = NumStatVar ? double( n_fixed ) / NumStatVar : 0.0;
 double maxc = 0;
 for( const auto & cd : cands )
  maxc = std::max( maxc , std::abs( cd.cost ) );
 if( maxc <= 0 )
  maxc = 1;
 const auto & lagr = get_Lagrangian_initial_solution();
 const Index n = cands.size();

 std::vector< Features > x( n );
 for( Index c = 0 ; c < n ; ++c ) {
  const auto & cd = cands[ c ];
  const auto n_sb = idx_to_var_sbi1[ cd.sb ].size();
  const double lg = cd.k < lagr.size() ? lagr[ cd.k ] : cd.value;
  x[ c ] = { cd.frac , cd.value , lg , std::abs( cd.value - lg ) ,
	     cd.cost / maxc ,
	     n_sb > 1 ? double( cd.pos ) / ( n_sb - 1 ) : 0.0 , fixed ,
	     n > 1 ? double( c ) / ( n - 1 ) : 0.0 };
  }
 return( x );

 }  // end( LagrangianDualRelaxationSolver::features )

/*--------------------------------------------------------------------------*/

std::vector< LagrangianDualRelaxationSolver::StrongEval >
LagrangianDualRelaxationSolver::strong_branching(
					       const std::vector< Cand > & cands )
{
 // the improvement of a child bound over that of the node, on the side of
 // the relaxation; at least eps, so that a candidate with one child not
 // improving still ranks by the other one
 const double node = valid_bound;
 const double eps = 1e-6 * std::max( std::abs( node ) , 1.0 );
 auto gain = [ & ]( double child ) {
  const double g = f_max ? node - child : child - node;
  return( std::max( g , eps ) );
  };

 std::vector< StrongEval > evals( cands.size() );
 for( Index c = 0 ; c < cands.size() ; ++c ) {
  auto changes = branchings( cands[ c ].var , cands[ c ].value );
  const double down = child_bound( changes[ 0 ] );
  const double up = child_bound( changes[ 1 ] );
  for( auto ch : changes )
   delete ch;

  evals[ c ] = { down , up , gain( down ) * gain( up ) };
  if( f_log && ( logVerb >= 1 ) )
   *f_log << "LagrangianDualRelaxationSolver::strong_branching: candidate "
	  << c << " value " << cands[ c ].value << " down " << down
	  << " up " << up << " score " << evals[ c ][ 2 ] << std::endl;
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

 const bool is_new = ! std::filesystem::exists( f_strong_log );
 std::ofstream out( f_strong_log , std::ios::app );
 if( ! out )
  throw( std::runtime_error( "LagrangianDualRelaxationSolver::"
			     "write_strong_log: cannot open " +
			     f_strong_log ) );
 if( is_new )
  out << "call,node,fixed,rank,sb,pos,value,lagr,frac,cost,down,up,"
      << "score,chosen" << std::endl;
 const auto x = features( cands );
 out << std::setprecision( 10 );
 for( Index c = 0 ; c < cands.size() ; ++c ) {
  const auto & cd = cands[ c ];
  out << f_n_branch << "," << valid_bound << "," << x[ c ][ 6 ] << "," << c
      << "," << cd.sb << "," << x[ c ][ 5 ] << "," << cd.value << ","
      << x[ c ][ 2 ] << "," << cd.frac << "," << cd.cost << ","
      << evals[ c ][ 0 ] << "," << evals[ c ][ 1 ] << "," << evals[ c ][ 2 ]
      << "," << ( c == chosen ? 1 : 0 ) << std::endl;
  }

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

 Index best_c = 0;
 if( f_branch_strategy == eStrongBranching ) {
  const auto evals = strong_branching( cands );
  for( Index c = 1 ; c < cands.size() ; ++c )
   if( evals[ c ][ 2 ] > evals[ best_c ][ 2 ] )
    best_c = c;
  write_strong_log( cands , evals , best_c );
  return( best_c );
  }

 // eLearned: the candidate the model of strBranchModel scores best
 if( f_w1.empty() )
  load_branch_model();
 const auto x = features( cands );
 double best_score = - Inf< double >();
 for( Index c = 0 ; c < cands.size() ; ++c ) {
  double score = f_b2;
  for( Index h = 0 ; h < f_w1.size() ; ++h ) {
   double a = f_b1[ h ];
   for( Index i = 0 ; i < NFeatures ; ++i )
    a += f_w1[ h ][ i ] * x[ c ][ i ];
   score += f_w2[ h ] * std::tanh( a );
   }
  if( score > best_score ) {
   best_score = score;
   best_c = c;
   }
  }
 return( best_c );

 }  // end( LagrangianDualRelaxationSolver::choose )

/*--------------------------------------------------------------------------*/

std::vector< bool > LagrangianDualRelaxationSolver::branchable( void ) const
{
 std::vector< bool > mask( NumStatVar , f_branch_groups.empty() );
 if( f_branch_groups.empty() )
  return( mask );

 // the variables of the named groups of each sub-Block, then the flags in
 // the order in which PrimalProximalHeur has registered them
 std::unordered_set< const ColVariable * > vars;
 for( const auto & sbi : f_Block->get_nested_Blocks() )
  for( const auto & group : sbi->get_static_variable_groups() )
   if( group && ( std::find( f_branch_groups.begin() , f_branch_groups.end() ,
			     group->get_name() ) != f_branch_groups.end() ) )
    group->for_each_as< ColVariable >( [ & ]( ColVariable & var ) {
      vars.insert( & var );
      } );

 Index k = 0;
 for( const auto & sbv : idx_to_var_sbi1 )
  for( const auto & dv : sbv )
   mask[ k++ ] = vars.contains( dv.second );

 return( mask );

 }  // end( LagrangianDualRelaxationSolver::branchable )

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
 auto cntr = std::make_shared< LagrangianChange::Center >();
 for( auto lbf : v_LBF )
  for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i )
   if( auto y = static_cast< ColVariable * >( lbf->get_active_var( i ) ) )
    cntr->emplace_back( y , y->get_value() );
 return( cntr );

 }  // end( LagrangianDualRelaxationSolver::center )

/*--------------------------------------------------------------------------*/

bool LagrangianDualRelaxationSolver::restore_center(
				     const LagrangianChange::Center & center )
{
 const auto par = InnerSolver->int_par_str2idx( "intRstAlg" );
 if( par == Inf< idx_type >() )
  return( false );

 std::unordered_map< ColVariable * , double > value;
 for( auto lbf : v_LBF )
  for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i )
   if( auto y = static_cast< ColVariable * >( lbf->get_active_var( i ) ) )
    value[ y ] = 0;
 for( const auto & [ y , v ] : center )
  if( auto it = value.find( y ) ; it != value.end() )
   it->second = v;
 for( const auto & [ y , v ] : value )
  y->set_value( v );

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

void LagrangianDualRelaxationSolver::load_branch_model( void )
{
 static const std::string _prfx =
  "LagrangianDualRelaxationSolver::load_branch_model: ";

 if( f_branch_model.empty() )
  throw( std::invalid_argument( _prfx + "eLearned needs strBranchModel" ) );
 std::ifstream in( f_branch_model );
 if( ! in )
  throw( std::runtime_error( _prfx + "cannot open " + f_branch_model ) );

 std::string tag;
 Index nin , nhid;
 if( ! ( in >> tag >> nin >> nhid ) || ( tag != "LBRModel" ) || ( nin != 8 )
     || ( nhid == 0 ) )
  throw( std::invalid_argument( _prfx + f_branch_model + " is not a model "
				"of 8 features" ) );

 f_w1.assign( nhid , std::vector< double >( nin ) );
 f_b1.resize( nhid );
 f_w2.resize( nhid );
 for( auto & row : f_w1 )
  for( auto & w : row )
   in >> w;
 for( auto & b : f_b1 )
  in >> b;
 for( auto & w : f_w2 )
  in >> w;
 in >> f_b2;
 if( ! in ) {
  f_w1.clear();
  throw( std::invalid_argument( _prfx + f_branch_model + " is truncated" ) );
  }

 }  // end( LagrangianDualRelaxationSolver::load_branch_model )

/*--------------------------------------------------------------------------*/

Change * LagrangianDualRelaxationSolver::apply( Change * change ,
						bool doUndo )
{
 auto c = dynamic_cast< LagrangianChange * >( change );
 if( ! c )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolver::apply: "
				"the Change is not a LagrangianChange" ) );

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
