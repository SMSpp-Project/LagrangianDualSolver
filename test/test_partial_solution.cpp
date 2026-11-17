/*--------------------------------------------------------------------------*/
/*--------------------- File test_partial_solution.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the part of the solution that get_var_solution() and
 * get_dual_solution() of LagrangianDualSolver write when they are given the
 * "meta" Configuration, which chooses Block by classname() [see
 * Solver::get_var_solution()]. The Block is
 *
 *     R    the row x + y + w <= 2, relaxed
 *     |- T    x , x <= 1
 *     |  `- L1    w , w <= 1
 *     `- L2    y , y <= 1
 *
 * with x , y , w >= 0 and the Objective of each Block maximizing its
 * Variable, where R is an AbstractBlock, T a TwigBlock and L1 and L2 are
 * LeafBlock, two AbstractBlock with a classname() of their own defined
 * here. The sub-Block of R, T and L2, are the components of the Lagrangian
 * Dual, each solved whole. After a first solve whose whole solution is the
 * reference, every value and every dual is overwritten with a sentinel and
 * then asked for again, and what is asked for has to be the reference and
 * the rest still the sentinel:
 *
 * - get_var_solution() with { TwigBlock : nullptr } writes T whole, x and
 *   w, with { LeafBlock : nullptr } T too, since L1 is in it, and L2, with
 *   { AbstractBlock : nullptr } or { * : nullptr } everything, with
 *   { NoSuchBlock : nullptr } nothing;
 *
 * - get_dual_solution() with { TwigBlock : nullptr } writes the rows of T
 *   and L1, with { LeafBlock : nullptr } those of L1 and L2, the Solver of
 *   T being passed the map and choosing L1 by itself, with
 *   { relaxed : nullptr } the relaxed row alone, with
 *   { AbstractBlock : nullptr } everything, with { NoSuchBlock : nullptr }
 *   nothing.
 *
 * The LagrangianDualSolver is the one of BSPar-sol.txt.
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
#include <map>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "CDASolver.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Meta = SimpleConfiguration< std::map< std::string , Configuration * > >;

static constexpr double SENTINEL = -77;

/*--------------------------------------------------------------------------*/
/// an AbstractBlock with a classname() of its own, for the "meta" format

class TwigBlock : public AbstractBlock {
 public:
 explicit TwigBlock( Block * father = nullptr ) : AbstractBlock( father ) {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_1( TwigBlock );

/// another one, for the leaves

class LeafBlock : public AbstractBlock {
 public:
 explicit LeafBlock( Block * father = nullptr ) : AbstractBlock( father ) {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_1( LeafBlock );

/*--------------------------------------------------------------------------*/
/// the Block and what the test reads of it, named by a letter

struct Problem {
 AbstractBlock * root;
 std::vector< ColVariable * > vars;         ///< x , w , y
 std::vector< FRowConstraint * > rows;      ///< x <= 1 , w <= 1 , y <= 1 ,
                                            ///< the relaxed one
 std::string var_name = "xwy";
 std::string row_name = "xwyr";
 };

/*--------------------------------------------------------------------------*/
/// adds to b a Variable x >= 0 , the row x <= 1 and max x

static void add_item( Problem & pb , AbstractBlock * b ,
                      const std::string & name )
{
 auto x = new ColVariable;
 x->set_type( ColVariable::kContinuous );
 x->is_positive( true , eNoMod );
 b->add_static_variable( * x , name.c_str() );

 auto row = new FRowConstraint;
 row->set_function( new LinearFunction(
                       LinearFunction::v_coeff_pair{ { x , 1.0 } } ) );
 row->set_lhs( - Inf< double >() );
 row->set_rhs( 1 );
 b->add_static_constraint( * row , ( name + "_cap" ).c_str() );

 auto of = new FRealObjective( b , new LinearFunction(
                       LinearFunction::v_coeff_pair{ { x , 1.0 } } ) );
 of->set_sense( Objective::eMax );
 b->set_objective( of );

 pb.vars.push_back( x );
 pb.rows.push_back( row );
 }

/*--------------------------------------------------------------------------*/

static Problem make_problem( void )
{
 Problem pb;
 pb.root = new AbstractBlock;
 auto T = new TwigBlock( pb.root );
 auto L1 = new LeafBlock( T );
 auto L2 = new LeafBlock( pb.root );
 pb.root->add_nested_Block( T );
 pb.root->add_nested_Block( L2 );
 T->add_nested_Block( L1 );

 add_item( pb , T , "x" );
 add_item( pb , L1 , "w" );
 add_item( pb , L2 , "y" );

 auto link = new FRowConstraint;
 link->set_function( new LinearFunction( LinearFunction::v_coeff_pair{
	{ pb.vars[ 0 ] , 1.0 } , { pb.vars[ 1 ] , 1.0 } ,
	{ pb.vars[ 2 ] , 1.0 } } ) );
 link->set_lhs( - Inf< double >() );
 link->set_rhs( 2 );
 pb.root->add_static_constraint( * link , "link" );
 pb.rows.push_back( link );

 auto of = new FRealObjective( pb.root , new LinearFunction() );
 of->set_sense( Objective::eMax );
 pb.root->set_objective( of );

 return( pb );
 }

/*--------------------------------------------------------------------------*/

static void overwrite( Problem & pb )
{
 for( auto v : pb.vars )
  v->set_value( SENTINEL );
 for( auto r : pb.rows )
  r->set_dual( SENTINEL );
 }

/*--------------------------------------------------------------------------*/
/// the names of the items that are no longer the sentinel

static std::string written( const std::string & names ,
                            const std::vector< double > & now )
{
 std::string w;
 for( std::size_t i = 0 ; i < now.size() ; ++i )
  if( now[ i ] != SENTINEL )
   w += names[ i ];
 return( w );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 auto bsc = dynamic_cast< BlockSolverConfig * >(
		    Configuration::deserialize( "BSPar-sol.txt" ) );
 if( ! bsc ) {
  std::cout << "BSPar-sol.txt does not hold a BlockSolverConfig"
            << std::endl;
  return( 1 );
  }

 auto pb = make_problem();
 bsc->apply( pb.root );
 bsc->clear();  // kept, clear()-ed, for the final cleanup

 auto solver = static_cast< CDASolver * >(
                                pb.root->get_registered_solvers().front() );
 int failures = 0;
 const int status = solver->compute();
 if( ( status != Solver::kOK ) || ( ! solver->has_var_solution() ) ||
     ( ! solver->has_dual_solution() ) ) {
  std::cout << "no solution to start from, status " << status << std::endl;
  ++failures;
  }
 else {
  auto vars = [ & ]() {
   std::vector< double > v;
   for( auto x : pb.vars )
    v.push_back( x->get_value() );
   return( v );
   };
  auto duals = [ & ]() {
   std::vector< double > v;
   for( auto r : pb.rows )
    v.push_back( r->get_dual() );
   return( v );
   };

  const std::vector< std::pair< std::string , std::string > > var_parts = {
   { "TwigBlock" , "xw" } , { "LeafBlock" , "xwy" } ,
   { "AbstractBlock" , "xwy" } , { "*" , "xwy" } , { "NoSuchBlock" , "" } };
  for( const auto & [ key , w ] : var_parts ) {
   Meta solc( { { key , nullptr } } );
   overwrite( pb );
   solver->get_var_solution( & solc );
   const auto got = written( pb.var_name , vars() );
   if( got != w ) {
    std::cout << "get_var_solution { " << key << " }: " << got
              << " written, not " << w << std::endl;
    ++failures;
    }
   }

  const std::vector< std::pair< std::string , std::string > > dual_parts = {
   { "TwigBlock" , "xw" } , { "LeafBlock" , "wy" } , { "relaxed" , "r" } ,
   { "AbstractBlock" , "xwyr" } , { "NoSuchBlock" , "" } };
  for( const auto & [ key , w ] : dual_parts ) {
   Meta solc( { { key , nullptr } } );
   overwrite( pb );
   solver->get_dual_solution( & solc );
   const auto got = written( pb.row_name , duals() );
   if( got != w ) {
    std::cout << "get_dual_solution { " << key << " }: " << got
              << " written, not " << w << std::endl;
    ++failures;
    }
   }
  }

 bsc->apply( pb.root );
 delete bsc;
 delete pb.root;

 if( failures )
  std::cout << failures << " checks failed" << std::endl;
 else
  std::cout << "All tests passed!!" << std::endl;
 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_partial_solution.cpp ------------------*/
/*--------------------------------------------------------------------------*/
