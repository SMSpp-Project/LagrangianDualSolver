/*--------------------------------------------------------------------------*/
/*---------------------------- File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of LagrangianDualSolver on Block built in memory.
 *
 * What is checked is the Lagrangian Dual that LagrangianDualSolver forms,
 * as it is handed to the inner Solver. This is ProbeCDASolver, defined
 * here: it solves nothing, and at each compute() it records, for each
 * Lagrangian multiplier, the coefficient in the Objective of the Lagrangian
 * Dual, the sign constraint, and the coefficient of the Variable of the
 * first component in the Lagrangian term. These have to be the ones that
 * the class comments of LagrangianDualSolver.h prescribe after (9), for an
 * equality, a <=, a >= and a free row, in a minimization and in a
 * maximization problem, with int_LDSlv_NNMult 0 and 1 (the two
 * LagrangianDualSolver of BSPar.txt).
 *
 * The cases are the edge ones of the relaxation of a row:
 *
 * - a ranged row, which has to be refused when the Solver is attached and
 *   when it is added as a dynamic row, with either value of
 *   int_LDSlv_NNMult;
 *
 * - changes of the sides that keep the type of a row, sequences whose
 *   intermediate state is ranged included, which change the coefficient of
 *   the multiplier alone;
 *
 * - changes that change the type of a row (a finite LHS given to a <= row,
 *   as setting a minimum budget does to a pollutant row of a UCBlock, the
 *   two sides of an equality moved apart, a side set to infinity, a row
 *   relaxed), which have to be refused by compute(), and again by each
 *   following call until the row gets its type back, after which the
 *   Lagrangian Dual has to be the right one;
 *
 * - dynamic rows added and removed, after which a change of a remaining
 *   one has to reach its own multiplier;
 *
 * - with intRecursive, the rows of the Block below the root, whose changes
 *   have to reach the Lagrangian Dual as those of the root do;
 *
 * - the constant term of the Objective of the root, which has to be that of
 *   the Objective of the Lagrangian Dual, also after it changes;
 *
 * - the conditional bound given to the Lagrangian Dual, i.e., the value of
 *   the components on their box with the opposite sense, which has to
 *   count the constant term of the Objective of the root (once, as that of
 *   a component), of either sign, also after it changes, and be infinite
 *   if a component is unbounded on its box;
 *
 * - intPushCostToOwner together with vintWhichPushCost, unsorted and with
 *   an index given twice, which has to give intPushCostToOwner 0 to the
 *   LagBFunction of the components listed and 1 to all the others.
 * - intRowScaling 1 (BSPar-scale.txt), with which the coefficients of each
 *   multiplier have to be those without scaling times the factor of its
 *   row, the multiplier the dual of the row divided by it, and the duals
 *   written back by get_dual_solution() the ones it was made from.
 * - intRowScaling 2 and 3 (BSPar-scale.txt), with Variable whose bounds
 *   come from a BoxConstraint, from a relaxed row and from nothing, with
 *   which the factors have to be sqrt( t_i / R_i ) and 1 / R_i normalized
 *   to geometric mean 1, and 1 for a row with R_i = 0.
 *
 * The test needs nothing but the core SMS++ and this module; the
 * Configuration files are read from the working directory.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractBlock.h"
#include "BlockSolverConfig.h"
#include "OneVarConstraint.h"
#include "CDASolver.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LagBFunction.h"
#include "LagrangianDualSolver.h"
#include "LinearFunction.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static const double INF = Inf< double >();

static int failures = 0;  // number of failed checks

/// what the inner Solver sees of one Lagrangian multiplier
struct Multiplier {
 double coef;  ///< the coefficient in the Objective of the Lagrangian Dual
 int sign;     ///< +1 if it is >= 0, -1 if it is <= 0, 0 if it is free
 double term;  ///< the coefficient of the Variable of the first component
               ///< in its Lagrangian term, 0 if it has none
 double value = 0;  ///< the value of the multiplier
 std::vector< double > terms;  ///< the coefficient of the first Variable
                               ///< of each component in its Lagrangian
                               ///< term, 0 if it has none
 };

/// what each compute() of ProbeCDASolver has seen, the last one at the end
static std::vector< std::vector< Multiplier > > Records;

/// the constant term of the Objective of the Lagrangian Dual at each
/// compute() of ProbeCDASolver, the last one at the end
static std::vector< double > Constants;

/// the conditional bound of the Lagrangian Dual at each compute() of
/// ProbeCDASolver, the upper one if it is a maximization problem and the
/// lower one otherwise, the last one at the end
static std::vector< double > Bounds;

/// intPushCostToOwner of the LagBFunction of each component at each
/// compute() of ProbeCDASolver, the last one at the end
static std::vector< std::vector< int > > PushCosts;

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ProbeCDASolver --------------------------*/
/*--------------------------------------------------------------------------*/
/// the inner Solver of the LagrangianDualSolver: records what it is given

class ProbeCDASolver : public CDASolver {
 public:
 ProbeCDASolver( void ) : CDASolver() {}
 ~ProbeCDASolver() override = default;

 int compute( bool changedvars = true ) override {
  std::vector< Multiplier > rec;
  auto obj = static_cast< FRealObjective * >( f_Block->get_objective() );
  auto lf = static_cast< LinearFunction * >( obj->get_function() );
  for( const auto & el : lf->get_v_var() ) {
   Multiplier m{ el.second , 0 , 0 };
   m.value = static_cast< ColVariable * >( el.first )->get_value();
   if( el.first->is_positive() )
    m.sign = 1;
   else
    if( el.first->is_negative() )
     m.sign = -1;
   for( auto sb : f_Block->get_nested_Blocks() ) {
    auto lbf = static_cast< LagBFunction * >( static_cast< FRealObjective * >(
					   sb->get_objective() )->get_function() );
    double term = 0;
    const auto i = lbf->is_active( el.first );
    if( i < Inf< Index >() )
     if( auto t = static_cast< LinearFunction * >(
					   lbf->get_Lagrangian_term( i ) ) )
      if( t->get_num_active_var() )
       term = t->get_v_var().front().second;
    m.terms.push_back( term );
    }
   m.term = m.terms.front();
   rec.push_back( m );
   }

  std::vector< int > push;
  for( auto sb : f_Block->get_nested_Blocks() )
   push.push_back(
     static_cast< LagBFunction * >(
       static_cast< FRealObjective * >( sb->get_objective() )->get_function() )
       ->get_int_par( LagBFunction::intPushCostToOwner ) );
  PushCosts.push_back( std::move( push ) );

  Records.push_back( std::move( rec ) );
  Constants.push_back( lf->get_constant_term() );
  Bounds.push_back( obj->get_sense() == Objective::eMax ?
		    f_Block->get_valid_upper_bound( true ) :
		    f_Block->get_valid_lower_bound( true ) );
  return( kOK );
  }

 bool has_var_solution( void ) override { return( false ); }
 void get_var_solution( Configuration * solc = nullptr ) override {}
 bool has_dual_solution( void ) override { return( false ); }
 void get_dual_solution( Configuration * solc = nullptr ) override {}
 void add_Modification( sp_Mod & mod ) override {}

 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_0( ProbeCDASolver );

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/
/// a relaxed row: lhs <= a0 x_0 + a1 x_1 <= rhs, with x_0 the Variable of
/// the first component, which is not in the row if a0 == 0; relaxed if it
/// is so (Constraint::relax())

struct Row {
 double lhs , rhs , a0 , a1;
 bool relaxed = false;
 };

/*--------------------------------------------------------------------------*/
/// the multiplier that the class comments of LagrangianDualSolver.h
/// prescribe for a row in a minimization (max == false) or maximization
/// problem, with int_LDSlv_NNMult nnmult

static Multiplier expected( const Row & d , bool max , bool nnmult )
{
 if( d.relaxed || ( ( d.lhs == -INF ) && ( d.rhs == INF ) ) )
  return( Multiplier{ 0 , 0 , 0 } );  // a free row, with no Lagrangian term

 if( d.lhs == d.rhs )       // an equality, with a free multiplier
  return( Multiplier{ - d.rhs , 0 , d.a0 } );

 const bool le = ( d.lhs == -INF );
 const double r = le ? d.rhs : d.lhs;
 // the natural sign: >= 0 for <= in a minimization and for >= in a
 // maximization, <= 0 for the other two
 const int sign = ( le != max ) ? 1 : -1;
 if( nnmult && ( sign < 0 ) )  // the row is reversed
  return( Multiplier{ r , 1 , - d.a0 } );

 return( Multiplier{ - r , sign , d.a0 } );
 }

/*--------------------------------------------------------------------------*/
/// the address of the Variable of a component

static ColVariable * var_of( Block * b )
{
 return( & ( *b->get_static_variable_v< ColVariable >( "x" ) )[ 0 ] );
 }

/*--------------------------------------------------------------------------*/
/// a component: one ColVariable with cost c

static AbstractBlock * component( double c , bool max )
{
 auto b = new AbstractBlock();
 auto x = new std::vector< ColVariable >( 1 );
 b->add_static_variable( *x , "x" );
 auto obj = new FRealObjective();
 obj->set_function( new LinearFunction(
		    LinearFunction::v_coeff_pair( { { & ( *x )[ 0 ] , c } } ) ) );
 obj->set_sense( max ? Objective::eMax : Objective::eMin , eNoMod );
 b->set_objective( obj );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/// writes the row d over the Variable of the Block x0 and x1 in r

static void set_row( FRowConstraint & r , const Row & d , Block * x0 ,
		     Block * x1 )
{
 LinearFunction::v_coeff_pair vp;
 if( d.a0 != 0 )
  vp.push_back( { var_of( x0 ) , d.a0 } );
 vp.push_back( { var_of( x1 ) , d.a1 } );
 r.set_function( new LinearFunction( std::move( vp ) ) , eNoMod );
 r.set_lhs( d.lhs , eNoMod );
 r.set_rhs( d.rhs , eNoMod );
 if( d.relaxed )
  r.relax( true , eNoMod );
 }

/*--------------------------------------------------------------------------*/
/// an Objective of the given sense with no Variable and the constant term
/// ct, as the root has to have

static void empty_objective( Block * b , bool max , double ct = 0 )
{
 auto obj = new FRealObjective();
 obj->set_function( new LinearFunction( {} , ct ) );
 obj->set_sense( max ? Objective::eMax : Objective::eMin , eNoMod );
 b->set_objective( obj );
 }

/*--------------------------------------------------------------------------*/
/// the root Block: two components and the rows over them, static in the
/// group "link" and dynamic (none at the beginning) in the group "dlink"

static AbstractBlock * root( const std::vector< Row > & rows , bool max ,
			     double ct = 0 )
{
 auto b = new AbstractBlock();
 b->add_nested_Block( component( 1 , max ) );
 b->add_nested_Block( component( 2 , max ) );

 auto link = new std::vector< FRowConstraint >( rows.size() );
 for( Index i = 0 ; i < rows.size() ; ++i )
  set_row( ( *link )[ i ] , rows[ i ] , b->get_nested_Block( 0 ) ,
	   b->get_nested_Block( 1 ) );
 b->add_static_constraint( *link , "link" );
 b->add_dynamic_constraint( *( new std::list< FRowConstraint >() ) ,
			    "dlink" );
 empty_objective( b , max , ct );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/// the static rows of the root

static std::vector< FRowConstraint > & link( Block * b )
{
 return( *b->get_static_constraint_v< FRowConstraint >( "link" ) );
 }

/*--------------------------------------------------------------------------*/
/// registers the Solver of the given BlockSolverConfig file to b, and
/// returns the (clear()-ed) BlockSolverConfig that unregisters them

static BlockSolverConfig * attach( Block * b , const std::string & file )
{
 auto c = Configuration::deserialize( file );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc ) {
  delete c;
  throw( std::invalid_argument( "cannot read " + file ) );
  }
 bsc->apply( b );
 bsc->clear();
 return( bsc );
 }

/*--------------------------------------------------------------------------*/
/// unregisters and deletes the Solver, then deletes the Block

static void detach( Block * b , BlockSolverConfig * bsc )
{
 bsc->apply( b );
 delete bsc;
 delete b;
 }

/*--------------------------------------------------------------------------*/
/// calls compute() of each Solver of b and compares what its inner Solver
/// sees with the rows rows; nnmult[ k ] is int_LDSlv_NNMult of the k-th
/// Solver (by default, those of BSPar.txt), and with term == false the
/// Lagrangian terms are not checked, the first component not being the
/// one of the rows

static void check_dual( Block * b , const std::vector< Row > & rows ,
			bool max , const std::string & what ,
			bool term = true ,
			const std::vector< bool > & nnmult_of = { false ,
								  true } )
{
 Index k = 0;
 for( auto s : b->get_registered_solvers() ) {
  const bool nnmult = nnmult_of[ k++ ];
  const std::string w = what + " (int_LDSlv_NNMult " +
                        std::to_string( nnmult ) + ")";
  try {
   s->compute();
   }
  catch( const std::exception & e ) {
   check( false , w + ": compute() throws " + e.what() );
   continue;
   }

  const auto & rec = Records.back();
  check( rec.size() == rows.size() , w + ": " +
	 std::to_string( rec.size() ) + " multipliers instead of " +
	 std::to_string( rows.size() ) );
  for( Index i = 0 ; ( i < rec.size() ) && ( i < rows.size() ) ; ++i ) {
   const auto e = expected( rows[ i ] , max , nnmult );
   const auto & m = rec[ i ];
   check( ( m.coef == e.coef ) && ( m.sign == e.sign ) &&
	  ( ( ! term ) || ( m.term == e.term ) ) ,
	  w + ": multiplier " + std::to_string( i ) + " has coefficient " +
	  std::to_string( m.coef ) + ", sign " + std::to_string( m.sign ) +
	  " and term " + std::to_string( m.term ) + " instead of " +
	  std::to_string( e.coef ) + ", " + std::to_string( e.sign ) +
	  " and " + std::to_string( e.term ) );
   }
  }
 }

/*--------------------------------------------------------------------------*/
/// calls compute() of each Solver of b, each of which has to throw an
/// exception of type E

template< class E >
static void check_refused( Block * b , const std::string & what )
{
 Index k = 0;
 for( auto s : b->get_registered_solvers() ) {
  const std::string w = what + " (Solver " + std::to_string( k++ ) + ")";
  bool refused = false;
  try {
   s->compute();
   }
  catch( const E & ) {
   refused = true;
   }
  catch( const std::exception & e ) {
   check( false , w + ": compute() throws the wrong exception " +
	  e.what() );
   continue;
   }
  check( refused , w + ": the change is not refused" );
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- TESTS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// a ranged row is refused when the Solver is attached, the first Solver
/// of BSPar.txt, which has int_LDSlv_NNMult 0, included

static void test_ranged_at_construction( void )
{
 auto b = root( { { 5 , 5 , 1 , 1 } , { 1 , 3 , 1 , 1 } } , false );
 bool refused = false;
 try {
  attach( b , "BSPar.txt" );
  }
 catch( const std::invalid_argument & ) {
  refused = true;
  }
 check( refused && b->get_registered_solvers().empty() ,
	"a ranged row is accepted when the Solver is attached" );

 // the Solver that threw is not registered, and the Block is not given
 // back to it: the test leaves both alone
 }

/*--------------------------------------------------------------------------*/
/// rows of each type, changed in their sides without changing type, in a
/// minimization and in a maximization problem

static void test_changes_keeping_type( bool max )
{
 const std::string p = max ? "max: " : "min: ";
 std::vector< Row > rows = { { 5 , 5 , 1 , 1 } ,       // =
			     { -INF , 7 , 1 , 2 } ,    // <=
			     { -3 , INF , 1 , -1 } ,   // >=
			     { -INF , INF , 3 , 1 } ,  // free
			     { -INF , 4 , 2 , 1 , true } };  // relaxed
 auto b = root( rows , max );
 auto bsc = attach( b , "BSPar.txt" );
 check_dual( b , rows , max , p + "the rows as they are born" );

 auto & l = link( b );

 // one side at a time, and both
 l[ 0 ].set_both( 6 );
 rows[ 0 ].lhs = rows[ 0 ].rhs = 6;
 l[ 1 ].set_rhs( 8 );
 rows[ 1 ].rhs = 8;
 l[ 2 ].set_lhs( -2 );
 rows[ 2 ].lhs = -2;
 l[ 4 ].set_rhs( 9 );  // a relaxed row stays without a Lagrangian term
 rows[ 4 ].rhs = 9;
 check_dual( b , rows , max , p + "one side changed" );

 // through a ranged state that is gone by the time compute() is called:
 // the two sides of the equality one after the other, a finite LHS given
 // to the <= row and taken away, a finite RHS given to the >= row and
 // taken away
 l[ 0 ].set_rhs( 9 );
 l[ 0 ].set_lhs( 9 );
 rows[ 0 ].lhs = rows[ 0 ].rhs = 9;
 l[ 1 ].set_lhs( 1 );
 l[ 1 ].set_rhs( 10 );
 l[ 1 ].set_lhs( -INF );
 rows[ 1 ].rhs = 10;
 l[ 2 ].set_rhs( 5 );
 l[ 2 ].set_lhs( -4 );
 l[ 2 ].set_rhs( INF );
 rows[ 2 ].lhs = -4;
 check_dual( b , rows , max , p + "sides changed through a ranged row" );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// each change of type is refused, again at the second call, and the
/// Lagrangian Dual is the right one once the type is back

static void test_changes_of_type( bool max )
{
 const std::string p = max ? "max: " : "min: ";
 std::vector< Row > rows = { { 5 , 5 , 1 , 1 } ,       // =
			     { -INF , 7 , 1 , 2 } ,    // <=
			     { -3 , INF , 1 , -1 } ,   // >=
			     { -INF , INF , 3 , 1 } };  // free
 auto b = root( rows , max );
 auto bsc = attach( b , "BSPar.txt" );
 check_dual( b , rows , max , p + "the rows as they are born" );

 auto & l = link( b );

 // a function that refuses the change twice, undoes it with undo, and
 // checks that the Lagrangian Dual is the one of rows
 const auto refused = [ & ]( const std::string & what , auto undo ) {
  check_refused< std::logic_error >( b , p + what );
  check_refused< std::logic_error >( b , p + what + ", second call" );
  undo();
  check_dual( b , rows , max , p + what + ", undone" );
  };

 // a finite LHS given to the <= row, which becomes ranged
 l[ 1 ].set_lhs( 2 );
 refused( "<= row made ranged" , [ & ]() { l[ 1 ].set_lhs( -INF ); } );

 // the same row made an equality
 l[ 1 ].set_lhs( 7 );
 refused( "<= row made an equality" , [ & ]() { l[ 1 ].set_lhs( -INF ); } );

 // the same row turned into a >= one
 l[ 1 ].set_lhs( 3 );
 l[ 1 ].set_rhs( INF );
 refused( "<= row turned into a >= one" ,
	  [ & ]() { l[ 1 ].set_rhs( 7 ); l[ 1 ].set_lhs( -INF ); } );

 // the same row made free
 l[ 1 ].set_rhs( INF );
 refused( "<= row made free" , [ & ]() { l[ 1 ].set_rhs( 7 ); } );

 // the same row relaxed
 l[ 1 ].relax( true );
 refused( "<= row relaxed" , [ & ]() { l[ 1 ].relax( false ); } );

 // a finite RHS given to the >= row, which becomes ranged
 l[ 2 ].set_rhs( 1 );
 refused( ">= row made ranged" , [ & ]() { l[ 2 ].set_rhs( INF ); } );

 // the two sides of the equality moved apart
 l[ 0 ].set_lhs( 4 );
 refused( "equality made ranged" , [ & ]() { l[ 0 ].set_lhs( 5 ); } );

 // a side of the equality set to infinity
 l[ 0 ].set_rhs( INF );
 refused( "equality made a >= row" , [ & ]() { l[ 0 ].set_rhs( 5 ); } );

 // a side given to the free row
 l[ 3 ].set_rhs( 1 );
 refused( "free row made a <= one" , [ & ]() { l[ 3 ].set_rhs( INF ); } );

 // and finally a change that is not undone, but done again in the other
 // direction: an equality that is not the original one
 l[ 0 ].set_lhs( 2 );
 check_refused< std::logic_error >( b , p + "equality made ranged again" );
 l[ 0 ].set_rhs( 2 );
 rows[ 0 ].lhs = rows[ 0 ].rhs = 2;
 check_dual( b , rows , max , p + "equality moved to another value" );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// the constant term of the Objective of the root (which has no Variable)
/// is that of the Objective of the Lagrangian Dual, when the Solver is
/// attached and after it changes

static void test_objective_constant( bool max )
{
 const std::string p = max ? "max: " : "min: ";
 std::vector< Row > rows = { { -INF , 7 , 1 , 2 } , { 5 , 5 , 1 , 1 } };
 auto b = root( rows , max , 7 );
 auto bsc = attach( b , "BSPar.txt" );

 const auto check_constant = [ & ]( double ct , const std::string & what ) {
  check_dual( b , rows , max , p + what );
  // the last two compute() are those of the two Solver
  for( Index k = 0 ; k < 2 ; ++k ) {
   const auto c = Constants[ Constants.size() - 2 + k ];
   check( c == ct , p + what + " (Solver " + std::to_string( k ) +
	  "): the constant term of the Lagrangian Dual is " +
	  std::to_string( c ) + " instead of " + std::to_string( ct ) );
   }
  };

 check_constant( 7 , "the constant term of the root" );

 static_cast< LinearFunction * >( static_cast< FRealObjective * >(
		   b->get_objective() )->get_function() )->set_constant_term( 9 );
 link( b )[ 0 ].set_rhs( 8 );
 rows[ 0 ].rhs = 8;
 check_constant( 9 , "the constant term of the root changed" );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// dynamic rows: removed one at a time, a range of them and a subset of
/// them, one of which has no term in the first component, and changed after
/// that; a ranged one is refused, and once it is removed it is as if it had
/// never been

static void test_dynamic_rows( void )
{
 std::vector< Row > rows = { { -INF , 7 , 1 , 2 } };
 auto b = root( rows , false );
 auto bsc = attach( b , "BSPar.txt" );
 auto dlink = b->get_dynamic_constraint< FRowConstraint >( "dlink" );

 const auto add = [ & ]( const std::vector< Row > & nr ) {
  std::list< FRowConstraint > nc( nr.size() );
  auto ncit = nc.begin();
  for( const auto & d : nr )
   set_row( *( ncit++ ) , d , b->get_nested_Block( 0 ) ,
	    b->get_nested_Block( 1 ) );
  b->add_dynamic_constraints( *dlink , nc );
  rows.insert( rows.end() , nr.begin() , nr.end() );
  };

 // five dynamic rows, the third one without the first component, which
 // therefore has no Lagrangian term there
 add( { { 4 , 4 , 1 , 3 } ,        // =
	{ -1 , INF , 2 , 1 } ,     // >=
	{ -INF , 6 , 0 , 1 } ,     // <=, second component only
	{ -INF , 5 , -1 , 1 } ,    // <=
	{ 2 , INF , 3 , 1 } } );   // >=
 check_dual( b , rows , false , "five dynamic rows added" );

 // the first two removed, which is a range of multipliers, and the fourth
 // one changed in its RHS
 b->remove_dynamic_constraint( *dlink , dlink->begin() );
 b->remove_dynamic_constraint( *dlink , dlink->begin() );
 std::next( dlink->begin() )->set_rhs( 8 );
 rows = { rows[ 0 ] , rows[ 3 ] , rows[ 4 ] , rows[ 5 ] };
 rows[ 2 ].rhs = 8;
 check_dual( b , rows , false , "a range of dynamic rows removed" );

 // the one without the first component removed alone, and the last one
 // changed in its LHS
 b->remove_dynamic_constraint( *dlink , dlink->begin() );
 dlink->back().set_lhs( 1 );
 rows = { rows[ 0 ] , rows[ 2 ] , rows[ 3 ] };
 rows[ 2 ].lhs = 1;
 check_dual( b , rows , false , "a dynamic row with no term removed" );

 // two more added, then the first and the third of the four removed, which
 // is not a range
 add( { { 0 , 0 , 1 , -1 } ,       // =
	{ -INF , 9 , 2 , 0.5 } } );  // <=
 check_dual( b , rows , false , "two more dynamic rows added" );
 b->remove_dynamic_constraint( *dlink , dlink->begin() );
 b->remove_dynamic_constraint( *dlink , std::next( dlink->begin() ) );
 rows = { rows[ 0 ] , rows[ 2 ] , rows[ 4 ] };
 check_dual( b , rows , false , "a subset of dynamic rows removed" );

 // a change of type of a dynamic row is refused as that of a static one
 dlink->front().set_rhs( 3 );
 check_refused< std::logic_error >( b , "a >= dynamic row made ranged" );
 dlink->front().set_rhs( INF );
 check_dual( b , rows , false , "the >= dynamic row back as it was" );

 // a ranged dynamic row, refused with either value of int_LDSlv_NNMult
 std::list< FRowConstraint > rc( 1 );
 set_row( rc.front() , { 1 , 3 , 1 , 1 } , b->get_nested_Block( 0 ) ,
	  b->get_nested_Block( 1 ) );
 b->add_dynamic_constraints( *dlink , rc );
 check_refused< std::invalid_argument >( b , "a ranged dynamic row" );

 // removed, it is as if it had never been
 b->remove_dynamic_constraint( *dlink , std::prev( dlink->end() ) );
 check_dual( b , rows , false , "the ranged dynamic row removed" );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// with intRecursive, the rows of the two sons of the root, each over the
/// two components below it, are relaxed after those of the root, and their
/// changes reach the Lagrangian Dual

static void test_recursive( void )
{
 auto b = new AbstractBlock();
 for( Index k = 0 ; k < 2 ; ++k ) {
  auto s = new AbstractBlock();
  s->add_nested_Block( component( 1 , false ) );
  s->add_nested_Block( component( 2 , false ) );
  b->add_nested_Block( s );
  }

 auto leaf = [ & ]( Index k , Index h ) {
  return( b->get_nested_Block( k )->get_nested_Block( h ) );
  };

 std::vector< Row > rows = { { -INF , 8 , 1 , 1 } ,   // the root
			     { 2 , INF , 1 , 1 } ,    // the first son
			     { 3 , 3 , 1 , 1 } };     // the second son

 auto rl = new std::vector< FRowConstraint >( 1 );
 set_row( ( *rl )[ 0 ] , rows[ 0 ] , leaf( 0 , 0 ) , leaf( 1 , 0 ) );
 b->add_static_constraint( *rl , "link" );
 empty_objective( b , false );

 for( Index k = 0 ; k < 2 ; ++k ) {
  auto s = b->get_nested_Block( k );
  auto sl = new std::vector< FRowConstraint >( 1 );
  set_row( ( *sl )[ 0 ] , rows[ 1 + k ] , leaf( k , 0 ) , leaf( k , 1 ) );
  static_cast< AbstractBlock * >( s )->add_static_constraint( *sl ,
							      "link" );
  }

 auto bsc = attach( b , "BSPar-rec.txt" );
 check_dual( b , rows , false , "recursive: the rows as they are born" ,
	     false , { true } );

 // a change of a row of each son
 link( b->get_nested_Block( 0 ) )[ 0 ].set_lhs( 1 );
 rows[ 1 ].lhs = 1;
 link( b->get_nested_Block( 1 ) )[ 0 ].set_both( 4 );
 rows[ 2 ].lhs = rows[ 2 ].rhs = 4;
 check_dual( b , rows , false , "recursive: the rows of the sons changed" ,
	     false , { true } );

 // and a change of type there is refused
 link( b->get_nested_Block( 0 ) )[ 0 ].set_rhs( 6 );
 check_refused< std::logic_error >( b ,
				    "recursive: a row of a son made ranged" );
 link( b->get_nested_Block( 0 ) )[ 0 ].set_rhs( INF );
 check_dual( b , rows , false , "recursive: the row of the son back" ,
	     false , { true } );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// with intPushCostToOwner 0 and vintWhichPushCost { 3 , 1 , 3 } (those of
/// BSPar-push.txt), the LagBFunction of the components 1 and 3 of 5 have
/// intPushCostToOwner 0 and the others 1

static void test_push_cost( void )
{
 auto b = new AbstractBlock();
 for( Index h = 0 ; h < 5 ; ++h )
  b->add_nested_Block( component( 1 + h , false ) );

 auto rl = new std::vector< FRowConstraint >( 1 );
 set_row( ( *rl )[ 0 ] , { -INF , 8 , 1 , 1 } , b->get_nested_Block( 0 ) ,
	  b->get_nested_Block( 1 ) );
 b->add_static_constraint( *rl , "link" );
 empty_objective( b , false );

 auto bsc = attach( b , "BSPar-push.txt" );
 PushCosts.clear();
 b->get_registered_solvers().front()->compute();
 const std::vector< int > expect = { 1 , 0 , 1 , 0 , 1 };
 check( ( PushCosts.size() == 1 ) && ( PushCosts.back() == expect ) ,
	"intPushCostToOwner with vintWhichPushCost" );

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// a component whose Variable is in [ 0 , ub ] by a BoxConstraint (none if
/// ub is infinite), with cost c and constant term ct

static AbstractBlock * box_component( double c , bool max , double ub ,
				      double ct )
{
 auto b = new AbstractBlock();
 auto x = new std::vector< ColVariable >( 1 );
 b->add_static_variable( *x , "x" );
 if( ub < INF ) {
  auto box = new std::vector< BoxConstraint >( 1 );
  ( *box )[ 0 ].set_variable( & ( *x )[ 0 ] );
  ( *box )[ 0 ].set_lhs( 0 , eNoMod );
  ( *box )[ 0 ].set_rhs( ub , eNoMod );
  b->add_static_constraint( *box , "box" );
  }
 auto obj = new FRealObjective();
 obj->set_function( new LinearFunction(
	       LinearFunction::v_coeff_pair( { { & ( *x )[ 0 ] , c } } ) , ct ) );
 obj->set_sense( max ? Objective::eMax : Objective::eMin , eNoMod );
 b->set_objective( obj );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/// the conditional bound of the Lagrangian Dual: the components x_0 in
/// [ 0 , 3 ] with cost 1 and constant term 2 and x_1 in [ 0 , 3 ] (or
/// unbounded) with cost 2, the root with the constant term ct and the rows
/// x_0 + x_1 <= 4 and x_0 - x_1 = 1. The bound is the value of the
/// components on their box with the sense opposite to that of the problem,
/// i.e., 2 + 3 + 6 for a minimization and 2 + 0 + 0 for a maximization,
/// plus ct, which is not in the components; the constant term 2 of x_0 is
/// counted once. It follows a change of ct, and is infinite (+ for a
/// minimization, - for a maximization) if x_1 is unbounded.

static void test_box_bound_constant( bool max )
{
 const std::string p = max ? "max: " : "min: ";
 std::vector< Row > rows = { { -INF , 4 , 1 , 1 } , { 1 , 1 , 1 , -1 } };

 for( double ct : { 7.0 , -4.0 , 0.0 } )
  for( bool bounded : { true , false } ) {
   const auto w = p + "box bound, root constant " + std::to_string( ct ) +
                  ( bounded ? "" : ", x_1 unbounded" );
   auto b = new AbstractBlock();
   b->add_nested_Block( box_component( 1 , max , 3 , 2 ) );
   b->add_nested_Block( box_component( 2 , max , bounded ? 3 : INF , 0 ) );
   auto l = new std::vector< FRowConstraint >( rows.size() );
   for( Index i = 0 ; i < rows.size() ; ++i )
    set_row( ( *l )[ i ] , rows[ i ] , b->get_nested_Block( 0 ) ,
	     b->get_nested_Block( 1 ) );
   b->add_static_constraint( *l , "link" );
   empty_objective( b , max , ct );
   auto bsc = attach( b , "BSPar.txt" );

   const double box = max ? 2 : 2 + 3 + 6;
   const double none = max ? - INF : INF;
   const auto check_bound = [ & ]( double e , const std::string & what ) {
    check_dual( b , rows , max , what );
    // the last two compute() are those of the two Solver
    for( Index k = 0 ; k < 2 ; ++k ) {
     const auto v = Bounds[ Bounds.size() - 2 + k ];
     check( v == e , what + " (Solver " + std::to_string( k ) + "): the "
	    "bound of the Lagrangian Dual is " + std::to_string( v ) +
	    " instead of " + std::to_string( e ) );
     }
    };

   check_bound( bounded ? box + ct : none , w );

   static_cast< LinearFunction * >( static_cast< FRealObjective * >(
		b->get_objective() )->get_function() )->set_constant_term(
								     ct + 5 );
   check_bound( bounded ? box + ct + 5 : none ,
		w + ", the root constant changed by 5" );

   detach( b , bsc );
   }
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// with intRowScaling 1 (the second of the four Solver of BSPar-scale.txt,
/// which have intRowScaling 0 ... 3) the
/// coefficient of each multiplier in the Objective and in the Lagrangian
/// terms is the one without scaling (the first Solver) times the factor of
/// the row, and the multiplier is the dual of the row divided by it. The
/// components are one with the Variable a and b and the row a + b <= 3, one
/// with c and the row c <= 0.5, and one each with d, e and g; the relaxed
/// rows are
///   a + c = 4, b + d <= 2, a - c = 0, b - d >= 0, e - g = 0
/// whose groups are { a , c }, { b , d } and { e , g }. The first two rows
/// have the factor 1 / |b| of their own; a - c = 0 takes the largest of
/// those of the rows that touch its group, 2 of c <= 0.5; b - d >= 0 the
/// 2 / 3 of a + b <= 3, which touches two groups and beats the 1 / 2 of
/// b + d <= 2; e - g = 0 has no row with a nonzero side in its group, and 1

static void test_row_scaling( void )
{
 const std::vector< double > alpha = { 0.25 , 0.5 , 2 , 2.0 / 3 , 1 };
 const std::vector< double > duals = { 1 , 2 , -3 , -4 , 5 };

 auto b = new AbstractBlock();
 auto comp0 = new AbstractBlock();
 auto ab = new std::vector< ColVariable >( 2 );
 comp0->add_static_variable( *ab , "x" );
 auto r0 = new std::vector< FRowConstraint >( 1 );
 ( *r0 )[ 0 ].set_function( new LinearFunction( LinearFunction::v_coeff_pair(
			  { { & ( *ab )[ 0 ] , 1 } , { & ( *ab )[ 1 ] , 1 } } ) ) ,
			    eNoMod );
 ( *r0 )[ 0 ].set_lhs( -INF , eNoMod );
 ( *r0 )[ 0 ].set_rhs( 3 , eNoMod );
 comp0->add_static_constraint( *r0 , "row" );
 auto obj = new FRealObjective();
 obj->set_function( new LinearFunction( LinearFunction::v_coeff_pair(
			  { { & ( *ab )[ 0 ] , 1 } , { & ( *ab )[ 1 ] , 1 } } ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 comp0->set_objective( obj );
 b->add_nested_Block( comp0 );

 for( Index h = 1 ; h < 5 ; ++h )
  b->add_nested_Block( component( double( h ) , false ) );

 auto comp1 = static_cast< AbstractBlock * >( b->get_nested_Block( 1 ) );
 auto r1 = new std::vector< FRowConstraint >( 1 );
 ( *r1 )[ 0 ].set_function( new LinearFunction( LinearFunction::v_coeff_pair(
					      { { var_of( comp1 ) , 1 } } ) ) ,
			    eNoMod );
 ( *r1 )[ 0 ].set_lhs( -INF , eNoMod );
 ( *r1 )[ 0 ].set_rhs( 0.5 , eNoMod );
 comp1->add_static_constraint( *r1 , "row" );

 const auto a = & ( *ab )[ 0 ];
 const auto bv = & ( *ab )[ 1 ];
 const auto c = var_of( comp1 );
 const auto d = var_of( b->get_nested_Block( 2 ) );
 const auto e = var_of( b->get_nested_Block( 3 ) );
 const auto g = var_of( b->get_nested_Block( 4 ) );

 struct LRow { ColVariable * x , * y; double sy , lhs , rhs; };
 const std::vector< LRow > rows = {
  { a , c , 1 , 4 , 4 } , { bv , d , 1 , -INF , 2 } , { a , c , -1 , 0 , 0 } ,
  { bv , d , -1 , 0 , INF } , { e , g , -1 , 0 , 0 } };

 auto lnk = new std::vector< FRowConstraint >( rows.size() );
 for( Index i = 0 ; i < rows.size() ; ++i ) {
  auto & r = ( *lnk )[ i ];
  r.set_function( new LinearFunction( LinearFunction::v_coeff_pair(
			{ { rows[ i ].x , 1 } , { rows[ i ].y , rows[ i ].sy } } ) ) ,
		  eNoMod );
  r.set_lhs( rows[ i ].lhs , eNoMod );
  r.set_rhs( rows[ i ].rhs , eNoMod );
  r.set_dual( duals[ i ] );
  }
 b->add_static_constraint( *lnk , "link" );
 empty_objective( b , false );

 auto bsc = attach( b , "BSPar-scale.txt" );
 std::vector< std::vector< Multiplier > > rec;
 for( auto s : b->get_registered_solvers() ) {
  s->compute();
  rec.push_back( Records.back() );
  }

 check( ( rec.size() == 4 ) && ( rec[ 0 ].size() == alpha.size() ) &&
	( rec[ 1 ].size() == alpha.size() ) ,
	"row scaling: not four Solver with 5 multipliers each" );

 const auto near = []( double x , double y ) {
  return( std::abs( x - y ) <= 1e-12 * std::max( 1.0 , std::abs( y ) ) );
  };

 for( Index i = 0 ; ( rec.size() == 4 ) && ( i < rec[ 0 ].size() ) &&
		    ( i < rec[ 1 ].size() ) ; ++i ) {
  const auto & u = rec[ 0 ][ i ];
  const auto & sc = rec[ 1 ][ i ];
  const std::string w = "row scaling: multiplier " + std::to_string( i );
  check( near( sc.coef , alpha[ i ] * u.coef ) ,
	 w + " has coefficient " + std::to_string( sc.coef ) +
	 " instead of " + std::to_string( alpha[ i ] * u.coef ) );
  check( sc.sign == u.sign , w + " has another sign constraint" );
  for( Index h = 0 ; h < u.terms.size() ; ++h )
   check( near( sc.terms[ h ] , alpha[ i ] * u.terms[ h ] ) ,
	  w + " has term " + std::to_string( sc.terms[ h ] ) +
	  " in component " + std::to_string( h ) + " instead of " +
	  std::to_string( alpha[ i ] * u.terms[ h ] ) );
  check( near( alpha[ i ] * sc.value , u.value ) ,
	 w + " has value " + std::to_string( sc.value ) + " instead of " +
	 std::to_string( u.value / alpha[ i ] ) );
  }

 // the duals written back are those the multipliers were made from
 for( auto s : b->get_registered_solvers() ) {
  for( auto & r : ( *lnk ) )
   r.set_dual( 0 );
  static_cast< CDASolver * >( s )->get_dual_solution();
  for( Index i = 0 ; i < ( *lnk ).size() ; ++i )
   check( near( ( *lnk )[ i ].get_dual() , duals[ i ] ) ,
	  "row scaling: the dual of row " + std::to_string( i ) +
	  " written back is " + std::to_string( ( *lnk )[ i ].get_dual() ) +
	  " instead of " + std::to_string( duals[ i ] ) );
  }

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/
/// with intRowScaling 2 and 3 (the third and fourth Solver of
/// BSPar-scale.txt) the factor of row i is sqrt( t_i / R_i ) and 1 / R_i,
/// normalized to geometric mean 1 over the rows that have one. The
/// components have one Variable each: x1 in [ 0 , 2 ] of cost 1, x2 in
/// [ 0 , 1 ] of cost 3, x3 in [ 0 , 4 ] of cost 2 and x4 free of cost 5;
/// the relaxed rows are
///   x1 + x2 = 1, 2 x1 - x3 <= 0, x2 + x4 >= 0, x3 - x4 = 0, x4 <= 10
/// and the first one implies x1 <= 1. So R is 2, 6, 1, 4 (x4 adds nothing)
/// and t is 2, 5 / 2, 2, 7, the spreads of the breakpoints { 1 , 3 },
/// { 1 / 2 , -2 }, { 3 , 5 } and { 2 , -5 }; the last row has R = 0, and 1

static void test_curvature_scaling( void )
{
 const std::vector< double > R = { 2 , 6 , 1 , 4 , 0 };
 const std::vector< double > t = { 2 , 2.5 , 2 , 7 , 0 };

 const auto normalized = []( std::vector< double > v ) {
  double lg = 0;
  Index n = 0;
  for( const auto x : v )
   if( x > 0 ) {
    lg += std::log( x );
    ++n;
    }
  const double gm = std::exp( lg / n );
  for( auto & x : v )
   x = ( x > 0 ) ? x / gm : 1;
  return( v );
  };
 std::vector< double > a2( R.size() , 0 ) , a3( R.size() , 0 );
 for( Index i = 0 ; i < R.size() ; ++i )
  if( R[ i ] > 0 ) {
   a2[ i ] = std::sqrt( t[ i ] / R[ i ] );
   a3[ i ] = 1 / R[ i ];
   }
 const std::vector< std::vector< double > > alpha = { normalized( a2 ) ,
						      normalized( a3 ) };

 auto b = new AbstractBlock();
 b->add_nested_Block( box_component( 1 , false , 2 , 0 ) );
 b->add_nested_Block( box_component( 3 , false , 1 , 0 ) );
 b->add_nested_Block( box_component( 2 , false , 4 , 0 ) );
 b->add_nested_Block( box_component( 5 , false , INF , 0 ) );
 std::vector< ColVariable * > x;
 for( Index h = 0 ; h < 4 ; ++h )
  x.push_back( var_of( b->get_nested_Block( h ) ) );

 struct LRow { Index j , k; double aj , ak , lhs , rhs; };
 const std::vector< LRow > rows = {
  { 0 , 1 , 1 , 1 , 1 , 1 } , { 0 , 2 , 2 , -1 , -INF , 0 } ,
  { 1 , 3 , 1 , 1 , 0 , INF } , { 2 , 3 , 1 , -1 , 0 , 0 } ,
  { 3 , 3 , 1 , 0 , -INF , 10 } };

 auto lnk = new std::vector< FRowConstraint >( rows.size() );
 for( Index i = 0 ; i < rows.size() ; ++i ) {
  const auto & d = rows[ i ];
  LinearFunction::v_coeff_pair vp = { { x[ d.j ] , d.aj } };
  if( d.ak != 0 )
   vp.push_back( { x[ d.k ] , d.ak } );
  auto & r = ( *lnk )[ i ];
  r.set_function( new LinearFunction( std::move( vp ) ) , eNoMod );
  r.set_lhs( d.lhs , eNoMod );
  r.set_rhs( d.rhs , eNoMod );
  r.set_dual( double( i + 1 ) );
  }
 static_cast< AbstractBlock * >( b )->add_static_constraint( *lnk ,
							      "link" );
 empty_objective( b , false );

 auto bsc = attach( b , "BSPar-scale.txt" );
 std::vector< std::vector< Multiplier > > rec;
 for( auto s : b->get_registered_solvers() ) {
  s->compute();
  rec.push_back( Records.back() );
  }

 check( ( rec.size() == 4 ) && ( rec[ 0 ].size() == rows.size() ) &&
	( rec[ 2 ].size() == rows.size() ) &&
	( rec[ 3 ].size() == rows.size() ) ,
	"curvature scaling: not four Solver with 5 multipliers each" );

 const auto near = []( double y , double z ) {
  return( std::abs( y - z ) <= 1e-12 * std::max( 1.0 , std::abs( z ) ) );
  };

 for( Index m = 0 ; ( rec.size() == 4 ) && ( m < 2 ) ; ++m )
  for( Index i = 0 ; ( i < rows.size() ) && ( i < rec[ 0 ].size() ) &&
		     ( i < rec[ m + 2 ].size() ) ; ++i ) {
   const auto & u = rec[ 0 ][ i ];
   const auto & sc = rec[ m + 2 ][ i ];
   const double al = alpha[ m ][ i ];
   const std::string w = "intRowScaling " + std::to_string( m + 2 ) +
			 ": multiplier " + std::to_string( i );
   check( near( sc.coef , al * u.coef ) ,
	  w + " has coefficient " + std::to_string( sc.coef ) +
	  " instead of " + std::to_string( al * u.coef ) );
   for( Index h = 0 ; h < u.terms.size() ; ++h )
    check( near( sc.terms[ h ] , al * u.terms[ h ] ) ,
	   w + " has term " + std::to_string( sc.terms[ h ] ) +
	   " in component " + std::to_string( h ) + " instead of " +
	   std::to_string( al * u.terms[ h ] ) );
   check( near( al * sc.value , u.value ) ,
	  w + " has value " + std::to_string( sc.value ) + " instead of " +
	  std::to_string( u.value / al ) );
   }

 detach( b , bsc );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 try {
  test_ranged_at_construction();
  test_changes_keeping_type( false );
  test_changes_keeping_type( true );
  test_changes_of_type( false );
  test_changes_of_type( true );
  test_objective_constant( false );
  test_objective_constant( true );
  test_dynamic_rows();
  test_recursive();
  test_box_bound_constant( false );
  test_box_bound_constant( true );
  test_push_cost();
  test_row_scaling();
  test_curvature_scaling();
  }
 catch( const std::exception & e ) {
  std::cout << "FAILED: exception " << e.what() << std::endl;
  return( 1 );
  }

 if( failures ) {
  std::cout << failures << " checks FAILED" << std::endl;
  return( 1 );
  }

 std::cout << "all checks passed" << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
