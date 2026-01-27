/*--------------------------------------------------------------------------*/
/*---------------------- File PrimalProximalHeur.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework to implement the
 * (so far, binary version only) of the Lagrangian-based Primal Proximal
 * heuristic proposed in
 *
 *  A. Daniilidis, C. Lemarechal "On a primal-proximal heuristic in
 *  discrete optimization" Mathematical Programming 104, 105-128, 2005
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PrimalProximalHeur
 #define __PrimalProximalHeur
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LagrangianDualSolver.h"

#include "Solution.h"

#include <queue>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class FRowConstraint;  // forward definition of FRowConstraint
  
/*--------------------------------------------------------------------------*/
/*-------------------- CLASS PrimalProximalHeur ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver solving the Lagrangian Dual of a "generic" Block
/** The PrimalProximalHeur class implements the CDASolver interface within
 * the SMS++ framework to implement the (so far, binary version only) of the
 * Lagrangian-based Primal Proximal heuristic proposed in
 *
 *  A. Daniilidis, C. Lemarechal "On a primal-proximal heuristic in
 *  discrete optimization" Mathematical Programming 104, 105-128, 2005
 *
 * The PrimalProximalHeur class is derived from the LagrangianDualSolver
 * class, since the Primal Proximal heuristic is based on adding a
 * quadratic regularization term that tries to push the convexified
 * solution produced by the Lagrangian dual to become integer (and hence
 * feasible), and then linearizes it (exploiting the binary nature of the
 * variables) so that it does a two-level change of the objective function
 * coefficients. */

class PrimalProximalHeur :  public LagrangianDualSolver 
{
  
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 using p_DQF = DQuadFunction *;
 using p_LF = LinearFunction *;

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that PrimalProximalHeur  has in addition to these of
  * LagrangianDualSolver. The value intLastPPHPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum int_par_type_PPH {
  intMaxIterPPH = intLastLDSlvPar ,  ///< maximum number of PPH iterations

  intLastPPHPar    ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of int algorithmic parameters. */

  };  // end( int_par_type_PPH )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PrimalProximalHeur has in addition to these of 
  * LagrangianDualSolver. The value  is provided so that the
  * list can be easily further extended by derived classes. */

 enum dbl_par_type_PPH {
  dbl_penaltyFactor = dblLastLDSlvPar , ///< PPH penalty factor

  dblLastLPPHPar ///< first allowed new double parameter for derived classes
                 /**< Convenience value for easily allow derived classes to
		  * extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_PPH )

/** @} ---------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING PrimalProximalHeur ------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PrimalProximalHeur
 *  @{ */

 /// constructor: ensure every field is initialized

 PrimalProximalHeur( void ) : LagrangianDualSolver() {
  logVerb = get_dflt_int_par( intLogVerb );
  maxIter = get_dflt_int_par( intMaxIterPPH );
  R = get_dflt_dbl_par( dbl_penaltyFactor );
  }

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~PrimalProximalHeur() {
  // not necessary, that of LagrangianDualSolver does it
  // set_Block( nullptr );
  guts_of_destructor();
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the PrimalProximalHeur has to solve
 /** Set (or changes) the Block that the PrimalProximalHeur has to solve. As
  * customary, set_Block() is a very important method where \p block is
  * thoroughly scanned (after being lock()-ed) and all relevant information is
  * extracted that the PrimalProximalHeur uses; this is even more crucial
  * here, since
  */

 void initialize( void );

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the int parameters of PrimalProximalHeur
 /** Set the int parameters specific of PrimalProximalHeur, and call the
  * version of LagrangianDualSolver for setting its own (which allows to
  * also set those of the inner Solver used to solve the Lagrangian Dual,
  * see the comments to set_ComputeConfig() for details). The following
  * parameters from Solver are managed:
  *
  * - intMaxIterPPH [Inf< int >()]: maximum number of Proximal Point Heuristic
  *                              iterations
  *
  * - intMaxSol [1]: maximum number of feasible solutions generated by the
  *                  Proximal Point Heuristic that are kept; the best ones
  *                  (in terms of objective value, of course) are kept
  *
  * - intLogVerb [0]: masks the first two bits of \p value ( & 3 ) and sets
  *                   the  verbosity of PrimalProximalHeur as
  *   = 0 : no log
  *   = 1 : detailed iteration-by-iteration log
  *   = 2 : even more detailed debug log
  *   then passes \p value >> 2 (shifted right 2 places, i.e., killing the
  *   first two bits) to LagrangianDualSolver
  *
  * Also, the following PrimalProximalHeur-specific parameters are managed:
  *
  * - intMaxIter [Inf< int >()]: passed as intMaxIter to
  *                                LagrangianDualSolver to set the maximum
  *                                number of iterations for each call of
  *   the inner Solver used to solve the Lagrangian Dual. */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of PrimalProximalHeur
 /** Set the double parameters specific of PrimalProximalHeur, and call the
  * version of LagrangianDualSolver for setting its own (which allows to
  * also set those of the inner Solver used to solve the Lagrangian Dual,
  * see the comments to set_ComputeConfig() for details). The parameters
  * are:
  *
  * - dbl_penaltyFactor [0]: the factor for the penalty term added to the
  *                          objective to try to force integrality. */

 void set_par( idx_type par , double value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the PrimalProximalHeur log

 void set_log( std::ostream * log_stream = nullptr ) override {
  f_log = log_stream;
  InnerSolver->set_log( f_log );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name running the Primal Proximal Heuristic
 *  @{ */

 /// run the Primal Proximal Heuristic on the given Block

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 OFValue get_lb( void ) override {
  return( f_max ? best_bound : - Inf< double >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override {
  return( f_max ? Inf< double >() : best_bound );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_var_solution( void ) override { 
  return( v_best_sol.empty() ? false : true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_var_solution( void ) override {
  if( v_best_sol.empty() )        // no current solution
   return( false );                   // no next solution either

  // v_best_sol is a heap with the best Solution at the top, i.e.,
  // the one with largest value if f_max and smallest otherwise; thus
  // use a custom comparator (it returns true if a is worse than b)
  std::pop_heap( v_best_sol.begin() , v_best_sol.end() ,
		 [ this ]( const sol_value & a , const sol_value & b ) {
		  return( f_max ? ( a.second < b.second ) :
			          ( a.second > b.second ) );
		 } );
  delete v_best_sol.back().first;   // delete current best solution
  v_best_sol.pop_back();            // remove current best solution
  if( v_best_sol.empty() ) {
   best_bound = f_max ? - Inf< double >() : Inf< double >();
   worst_bound = - best_bound;
   return( false );
   }

  best_bound = v_best_sol.front().second;
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_var_value( void ) override {
  if( v_best_sol.empty() )
   return( f_max ? - Inf< double >() : Inf< double >() );
  else
   return( v_best_sol.front().second );
  }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( v_best_sol.empty() )
   throw( std::logic_error( "PrimalProximalHeur::get_var_solution() called "
			    "with no available solution" ) );

  return( v_best_sol.front().first->write( f_Block ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_funct_value( void ) { 
  double value_funct = 0.0;
  remove_penalty_terms();

  if( !InnerSolver->has_var_solution() ){
    return( f_max ? - Inf< double >() : Inf< double >());
  } else {
    for( const auto & sbi : f_Block->get_nested_Blocks() ) {
      auto Funct_sbi_idx = static_cast< Function * >( 
        static_cast< FRealObjective * >( sbi->get_objective() )->get_function() );
      Funct_sbi_idx->compute( true );
      auto v = Funct_sbi_idx->get_value();
      ///! if( !isnan(v) )
        value_funct += v;
    }
  }
  return( value_funct ); 
 }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the PrimalProximalHeur
 *
 * Handle the few parameters of PrimalProximalHeur and redirect the others
 * to LagrangianDualSolver. Recall that the latter allows to change via its
 * standard parameter interface all the parameters of the "inner Solver"
 * used to actually solve the Lagrangian Dual, but the somewhat complicated
 * mechanism is completely transparent to PrimalProximalHeur. */

  [[nodiscard]] idx_type get_num_int_par( void ) const override {
   return( LagrangianDualSolver::get_num_int_par() + 1 );
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( LagrangianDualSolver::get_num_dbl_par() + 1 );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par == intMaxIterPPH )
   return( Inf< int >() );
  else
   return( LagrangianDualSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  if( par ==  dbl_penaltyFactor )
   return( 0 );
  else
   return( LagrangianDualSolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( intMaxSol ):  return( f_MaxSol );
   case( intLogVerb ):
    return( logVerb + ( LagrangianDualSolver::get_int_par( par ) << 2 ) );
   case( intMaxIterPPH ): return( maxIter );
   }
  return( LagrangianDualSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  if( par ==  dbl_penaltyFactor )
   return( R );
  else
   return( LagrangianDualSolver::get_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intMaxIterPPH" )
   return( intMaxIterPPH );
  else
   return( LagrangianDualSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  if( name == "dbl_penaltyFactor" )
   return( dbl_penaltyFactor );
  else
  return( LagrangianDualSolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string _ret = "intMaxIterPPH";
  if( idx == intMaxIterPPH )
   return( _ret );
  else
   return( LagrangianDualSolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {
  static const std::string _ret = "dbl_penaltyFactor";
  if( idx == dbl_penaltyFactor )
   return( _ret );
  else
   return( LagrangianDualSolver::dbl_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING THE State OF THE PrimalProximalHeur ---------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the State of the PrimalProximalHeur
 *  @{ */

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 using sol_value = std::pair< Solution * , double >;

 using double_var = std::pair< double , ColVariable * >;

 using var_col_int = std::tuple< ColVariable * , Index , Index >;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 int logVerb;
 ///< verbosity of PrimalProximalHeur is verbosity of InnerSolver + 2 

 Index maxIter;   ///< maximum number of iterations for PrimalProximalHeur
 Index f_MaxSol;  ///< maximum number of Solution to keep

 double R;
 double penalty = 0;      // penalty factor
 double addterm;

 double value_FUNCTION;   // objective function value for the current solution

 double best_bound;       // best bound of the feasibile solutions  
 double worst_bound;      // worst bound of the feasibile solutions  

 bool changed_penalties = false;  // true if penalty terms changed

 Index NumStatVar;      ///< (current) number of static variables

 Index pos_id;

 std::vector< int > pos_id_sbi;
 std::vector< double > previous_sol; 
 
 std::vector< Function * > Funct_sbi;
                        ///< vector of objective functions for sub-Block sbi

 std::vector<std::vector< double_var >> idx_to_var_sbi1;
 ///< from index to static variable (only linear terms)
 std::vector<std::vector< double_var >> idx_to_var_sbi2;
 ///< from index to static variable (only quadratic terms)

 std::vector< sol_value > v_best_sol;  ///< best feasible solutions
 /**< v_best_sol is managed as a binary heap */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/

 void add_penalty_terms( void );

/*--------------------------------------------------------------------------*/

 void remove_penalty_terms( void );

/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class PrimalProximalHeur )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* PrimalProximalHeur.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File PrimalProximalHeur.h -----------------------*/
/*--------------------------------------------------------------------------*/
