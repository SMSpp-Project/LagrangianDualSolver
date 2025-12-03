/*--------------------------------------------------------------------------*/
/*---------------------- File PrimalProximalHeur.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework.
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

#include <queue>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class FRowConstraint;  // forward definition of FRowConstraint
  
/*--------------------------------------------------------------------------*/
/*-------------------- CLASS PrimalProximalHeur --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver solving the Lagrangian Dual of a "generic" Block
/** The PrimalProximalHeur class implements the CDASolver interface within
 * the SMS++ framework.
 */

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

 // "import" basic types from Block
 using Index = Block::Index;
 using c_Index = Block::c_Index;

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;

 using p_DQF = DQuadFunction *;

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PrimalProximalHeur has in addition to these of CDASolver. The
  * value dblLastLDSSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum dbl_par_type_LDSlv {
  dblLastLDSlvPar = dblLastParCDAS ,
  dbl_penaltyFactor ,
  maxIterPP
   ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set
   * of double algorithmic parameters. */

  };  // end( dbl_par_type_LDSlv )

/** @} ---------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING PrimalProximalHeur ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PrimalProximalHeur
 *  @{ */

 /// constructor: ensure every field is initialized

 PrimalProximalHeur( void ) : LagrangianDualSolver() , best_bound( Inf< double >() ) ,
    best_solutions( ) ,
    R( 0.0 ) , logVerb( 2 ) , maxIter( 10.0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the PrimalProximalHeur from the Block

 virtual ~PrimalProximalHeur() { set_Block( nullptr ); }

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of PrimalProximalHeur / the inner Solver
 /** Set the double parameters specific of PrimalProximalHeur, and allow to
  * directly set those of the inner Solver used to solve the Lagrangian Dual;
  * see the comments to set_ComputeConfig() for details. */

 void set_par( idx_type par , double value ) override;
 void set_par( idx_type par , int value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the PrimalProximalHeur log

 void set_log( std::ostream * log_stream = nullptr ) override {
  f_log = log_stream;
  InnerSolver->set_log( f_log );
  /*!!
  for( auto lbf : v_LBF )
   for( auto s : lbf->get_inner_block()->get_registered_solvers() )
    s->set_log( f_log );
    !!*/
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Lagrangian Dual of the given Block
 *  @{ */

 /// (try to) solve the Lagrangian Dual of the given Block

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 OFValue get_lb( void ) override { return( f_max ? best_bound : - Inf<double>() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { return( f_max ? Inf<double>() : best_bound ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_var_solution( void ) override { 
  return( best_solutions.empty() ? false : true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_var_solution( void ) override { return( has_new_solution ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_var_value( void ) override {
  if( best_solutions.empty() ) return( best_solutions.top().second ); }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration * solc = nullptr ) override { return( 
    best_solutions.top().first->write( f_Block )); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_funct_value( void ) { return( value_FUNCTION ); }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the PrimalProximalHeur
 *
 * While PrimalProximalHeur itself has comparatively few parameters, it
 * allows to change via its standard parameter interface all the parameters
 * of the "inner Solver" used to actually solve the Lagrangian Dual. The
 * indices of these parameters are automatically translated (see *_par_is()
 * and the comments to set_ComputeConfig() for details) so that they can be
 * automatically set and queried as if they were "natural" parameters of
 * PrimalProximalHeur itself.
 * @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  switch( par ) {
    case( dbl_penaltyFactor ):     return( R );
    case( maxIterPP ):     return( maxIter );
  }
  return( InnerSolver->get_dbl_par( dbl_par_lds( par ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > dbl_pars_map = {
    { "dbl_penaltyFactor" , PrimalProximalHeur::dbl_penaltyFactor } ,
    { "maxIterPP" , PrimalProximalHeur::maxIterPP } ,
  };

  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );

  return( dbl_par_is( InnerSolver->dbl_par_str2idx( name ) ) );
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

 typedef std::pair< Solution * , double > sol_value;
 typedef std::pair< double , ColVariable * > double_var;
 typedef std::tuple< ColVariable * , Index , Index > var_col_int;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 int logVerb;
 double maxIter;

 double penalty = 0.0;
 double R;
 double value_FUNCTION;
 double value_FUNCTION1;
 double best_bound;

 bool changed_penalties = false;
 bool has_new_solution = false;

 Index NumStatVar;      ///< (current) number of static variables
 Index pos_id;

 std::vector<int> pos_id_sbi;
 std::vector< double > previous_sol; 

 std::vector< var_col_int > var_to_idx;   ///< from static variable to index
 std::vector< double_var > idx_to_var1;   ///< from index to static variable
 std::vector< double_var > idx_to_var2;   ///< from index to static variable
 std::vector<p_DQF> Funct_sbi;            ///< vector of objective functions for sub-block sbi

 std::vector<std::vector< double_var >> idx_to_var_sbi1;   ///< from index to static variable (linear term)
 std::vector<std::vector< double_var >> idx_to_var_sbi2;   ///< from index to static variable (quadratic)

 std::priority_queue< sol_value > best_solutions; ///< best feasible solutions 

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

 void add_penalty_terms();

/*--------------------------------------------------------------------------*/

 void remove_penalty_terms();

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
/*-------------------- End File PrimalProximalHeur.h ---------------------*/
/*--------------------------------------------------------------------------*/
