/*--------------------------------------------------------------------------*/
/*-------------- File LagrangianDualRelaxationSolverML.cpp -----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the class LagrangianDualRelaxationSolverML.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include <fstream>

#include <iomanip>

#include "LagrangianDualRelaxationSolverML.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// register LagrangianDualRelaxationSolverML to the Solver factory

SMSpp_insert_in_factory_cpp_0( LagrangianDualRelaxationSolverML );

/*--------------------------------------------------------------------------*/
/*------------- METHODS OF LagrangianDualRelaxationSolverML ----------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolverML::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intLearnNodes ):  f_learn_nodes = value; break;
  case( intTrainEpochs ): f_train_epochs = value; break;
  case( intHidden ):      f_hidden = value; break;
  case( intLearnSeed ):   f_learn_seed = value; break;
  default: LagrangianDualRelaxationSolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolverML::set_par( idx_type par , double value )
{
 if( par == dblLearnRate )
  f_learn_rate = value;
 else
  LagrangianDualRelaxationSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolverML::set_par( idx_type par ,
						std::string && value )
{
 if( par == strOnlineModel )
  f_online_model = std::move( value );
 else
  LagrangianDualRelaxationSolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolverML::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intLearnNodes ):  return( 20 );
  case( intTrainEpochs ): return( 300 );
  case( intHidden ):      return( 16 );
  case( intLearnSeed ):   return( 0 );
  default: return( LagrangianDualRelaxationSolver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

int LagrangianDualRelaxationSolverML::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intLearnNodes ):  return( f_learn_nodes );
  case( intTrainEpochs ): return( f_train_epochs );
  case( intHidden ):      return( f_hidden );
  case( intLearnSeed ):   return( f_learn_seed );
  default: return( LagrangianDualRelaxationSolver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolverML::int_par_str2idx(
					       const std::string & name ) const
{
 if( name == "intLearnNodes" )
  return( intLearnNodes );
 if( name == "intTrainEpochs" )
  return( intTrainEpochs );
 if( name == "intHidden" )
  return( intHidden );
 if( name == "intLearnSeed" )
  return( intLearnSeed );
 return( LagrangianDualRelaxationSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolverML::int_par_idx2str(
						       idx_type par ) const
{
 static const std::string nodes = "intLearnNodes";
 static const std::string epochs = "intTrainEpochs";
 static const std::string hidden = "intHidden";
 static const std::string seed = "intLearnSeed";
 switch( par ) {
  case( intLearnNodes ):  return( nodes );
  case( intTrainEpochs ): return( epochs );
  case( intHidden ):      return( hidden );
  case( intLearnSeed ):   return( seed );
  default: return( LagrangianDualRelaxationSolver::int_par_idx2str( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolverML::get_dflt_dbl_par( idx_type par )
 const
{
 if( par == dblLearnRate )
  return( 1e-2 );
 return( LagrangianDualRelaxationSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualRelaxationSolverML::get_dbl_par( idx_type par ) const
{
 if( par == dblLearnRate )
  return( f_learn_rate );
 return( LagrangianDualRelaxationSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolverML::dbl_par_str2idx(
					       const std::string & name ) const
{
 if( name == "dblLearnRate" )
  return( dblLearnRate );
 return( LagrangianDualRelaxationSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolverML::dbl_par_idx2str(
						       idx_type par ) const
{
 static const std::string rate = "dblLearnRate";
 if( par == dblLearnRate )
  return( rate );
 return( LagrangianDualRelaxationSolver::dbl_par_idx2str( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolverML::get_dflt_str_par(
						       idx_type par ) const
{
 static const std::string empty;
 if( par == strOnlineModel )
  return( empty );
 return( LagrangianDualRelaxationSolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolverML::get_str_par(
						       idx_type par ) const
{
 if( par == strOnlineModel )
  return( f_online_model );
 return( LagrangianDualRelaxationSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type LagrangianDualRelaxationSolverML::str_par_str2idx(
					       const std::string & name ) const
{
 if( name == "strOnlineModel" )
  return( strOnlineModel );
 return( LagrangianDualRelaxationSolver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualRelaxationSolverML::str_par_idx2str(
						       idx_type par ) const
{
 static const std::string model = "strOnlineModel";
 if( par == strOnlineModel )
  return( model );
 return( LagrangianDualRelaxationSolver::str_par_idx2str( par ) );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolverML::set_global_information(
						    GlobalInformation * gi )
{
 LagrangianDualRelaxationSolver::set_global_information( gi );
 f_online = nullptr;
 if( ! gi )
  return;

 if( ! gi->exists( str_OnlineBranching ) )
  gi->add_to_Universe< OnlineData >( str_OnlineBranching );
 f_online = gi->get_from_Universe< OnlineData >( str_OnlineBranching );
 if( f_online && ( ! f_online->contains( str_OnlineData ) ) )
  f_online->write( str_OnlineData , OnlineData{} );

 }  // end( LagrangianDualRelaxationSolverML::set_global_information )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

LagrangianDualRelaxationSolverML::Index
LagrangianDualRelaxationSolverML::choose( std::vector< Cand > & cands )
{
 if( f_branch_strategy != eOnline )
  return( LagrangianDualRelaxationSolver::choose( cands ) );

 if( f_strong_cands <= 0 )
  throw( std::invalid_argument( "LagrangianDualRelaxationSolverML::choose: "
				"intStrongCands must be positive" ) );
 if( cands.size() > Index( f_strong_cands ) )
  cands.resize( f_strong_cands );
 const auto x = features( cands );

 // the model, if it is already fitted: the candidate it scores best
 torch::nn::Sequential net{ nullptr };
 with_data( [ & ]( OnlineData & d ) { net = d.net; } );
 if( ! net.is_empty() ) {
  torch::NoGradGuard no_grad;
  auto in = torch::empty( { long( x.size() ) , long( NFeatures ) } );
  auto acc = in.accessor< float , 2 >();
  for( Index c = 0 ; c < x.size() ; ++c )
   for( Index i = 0 ; i < NFeatures ; ++i )
    acc[ c ][ i ] = float( x[ c ][ i ] );
  return( Index( net->forward( in ).squeeze( 1 ).argmax().item< long >() ) );
  }

 // strong branching, whose choice and data are kept: the targets are the
 // scores over the largest one, an infinite score (a child found
 // infeasible) counting as ten times the largest finite one
 const auto evals = strong_branching( cands );
 Index best_c = 0;
 for( Index c = 1 ; c < cands.size() ; ++c )
  if( evals[ c ][ 2 ] > evals[ best_c ][ 2 ] )
   best_c = c;
 write_strong_log( cands , evals , best_c );

 if( cands.size() < 2 )  // nothing to rank
  return( best_c );

 double maxf = 0;
 for( const auto & e : evals )
  if( std::isfinite( e[ 2 ] ) )
   maxf = std::max( maxf , e[ 2 ] );
 if( maxf <= 0 )
  maxf = 1;
 std::vector< double > y( cands.size() );
 double maxy = 0;
 for( Index c = 0 ; c < cands.size() ; ++c ) {
  y[ c ] = std::isfinite( evals[ c ][ 2 ] ) ? evals[ c ][ 2 ] : 10 * maxf;
  maxy = std::max( maxy , y[ c ] );
  }
 for( auto & v : y )
  v /= maxy;

 // the node is added to the data, and the first Solver that brings them to
 // intLearnNodes nodes fits the model (out of the lock, on a copy of them)
 bool fit_now = false;
 std::vector< std::vector< Features > > xs;
 std::vector< std::vector< double > > ys;
 with_data( [ & ]( OnlineData & d ) {
   d.x.push_back( x );
   d.y.push_back( std::move( y ) );
   if( d.net.is_empty() && ( ! d.fitting ) &&
       ( d.x.size() >= Index( std::max( f_learn_nodes , 1 ) ) ) ) {
    d.fitting = fit_now = true;
    xs = d.x;
    ys = d.y;
    }
   } );

 if( fit_now ) {
  auto fitted = fit( xs , ys );
  if( ! f_online_model.empty() )
   write_model( fitted );
  with_data( [ & ]( OnlineData & d ) {
    d.net = fitted;
    d.fitting = false;
    } );
  }

 return( best_c );

 }  // end( LagrangianDualRelaxationSolverML::choose )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

template< class F >
void LagrangianDualRelaxationSolverML::with_data( F && f )
{
 if( f_online )
  f_online->write_with( str_OnlineData , std::forward< F >( f ) );
 else
  f( f_local );
 }

/*--------------------------------------------------------------------------*/

torch::nn::Sequential LagrangianDualRelaxationSolverML::fit(
			     const std::vector< std::vector< Features > > & x ,
			     const std::vector< std::vector< double > > & y )
{
 torch::manual_seed( f_learn_seed );
 const auto hid = long( std::max( f_hidden , 1 ) );
 torch::nn::Sequential net( torch::nn::Linear( long( NFeatures ) , hid ) ,
			    torch::nn::Tanh() ,
			    torch::nn::Linear( hid , 1 ) );
 torch::optim::Adam opt( net->parameters() ,
			 torch::optim::AdamOptions( f_learn_rate ) );

 // one pair of tensors per node with at least two candidates
 std::vector< std::pair< torch::Tensor , torch::Tensor > > data;
 for( Index g = 0 ; g < x.size() ; ++g ) {
  const auto n = long( x[ g ].size() );
  if( n < 2 )
   continue;
  auto in = torch::empty( { n , long( NFeatures ) } );
  auto t = torch::empty( { n } );
  auto acc = in.accessor< float , 2 >();
  auto tac = t.accessor< float , 1 >();
  for( long c = 0 ; c < n ; ++c ) {
   for( Index i = 0 ; i < NFeatures ; ++i )
    acc[ c ][ i ] = float( x[ g ][ c ][ i ] );
   tac[ c ] = float( y[ g ][ c ] );
   }
  data.emplace_back( in , t );
  }

 double last = 0;
 for( int epoch = 0 ; ( epoch < f_train_epochs ) && ( ! data.empty() ) ;
      ++epoch ) {
  opt.zero_grad();
  auto loss = torch::zeros( { 1 } );
  // listwise: the cross entropy between the distribution of the scores of
  // the model on the candidates of a node and that of their targets
  for( auto & [ in , t ] : data ) {
   auto p = torch::log_softmax( net->forward( in ).squeeze( 1 ) , 0 );
   auto q = torch::softmax( t * 10 , 0 );
   loss = loss - ( q * p ).sum();
   }
  loss = loss / double( data.size() );
  loss.backward();
  opt.step();
  last = loss.item< double >();
  }

 if( f_log && ( logVerb >= 1 ) )
  *f_log << "LagrangianDualRelaxationSolverML::fit: " << data.size()
	 << " nodes, " << f_train_epochs << " epochs, loss " << last
	 << std::endl;

 return( net );

 }  // end( LagrangianDualRelaxationSolverML::fit )

/*--------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolverML::write_model(
					    torch::nn::Sequential & net ) const
{
 std::ofstream out( f_online_model );
 if( ! out )
  throw( std::runtime_error( "LagrangianDualRelaxationSolverML::"
			     "write_model: cannot open " + f_online_model ) );

 torch::NoGradGuard no_grad;
 auto l1 = net[ 0 ]->as< torch::nn::Linear >();
 auto l2 = net[ 2 ]->as< torch::nn::Linear >();
 auto w1 = l1->weight.to( torch::kDouble );
 auto b1 = l1->bias.to( torch::kDouble );
 auto w2 = l2->weight.to( torch::kDouble );
 auto b2 = l2->bias.to( torch::kDouble );
 const auto hid = w1.size( 0 );

 out << "LBRModel " << NFeatures << " " << hid << std::endl
     << std::setprecision( 9 );
 for( long h = 0 ; h < hid ; ++h ) {
  for( Index i = 0 ; i < NFeatures ; ++i )
   out << ( i ? " " : "" ) << w1[ h ][ long( i ) ].item< double >();
  out << std::endl;
  }
 for( long h = 0 ; h < hid ; ++h )
  out << ( h ? " " : "" ) << b1[ h ].item< double >();
 out << std::endl;
 for( long h = 0 ; h < hid ; ++h )
  out << ( h ? " " : "" ) << w2[ 0 ][ h ].item< double >();
 out << std::endl << b2[ 0 ].item< double >() << std::endl;

 }  // end( LagrangianDualRelaxationSolverML::write_model )

/*--------------------------------------------------------------------------*/
/*------------ End File LagrangianDualRelaxationSolverML.cpp ---------------*/
/*--------------------------------------------------------------------------*/
