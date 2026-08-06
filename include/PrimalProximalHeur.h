/*--------------------------------------------------------------------------*/
/*---------------------- File PrimalProximalHeur.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the PrimalProximalHeur class, which implements the
 * CDASolver interface within the SMS++ framework as a primal-proximal
 * heuristic for binary optimization problems built on top of a
 * LagrangianDualSolver.
 *
 * The heuristic is described in
 *
 *   A. Daniilidis, C. Lemarechal "On a primal-proximal heuristic in
 *   discrete optimization" Mathematical Programming 104, 105-128, 2005
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

#include "ColVariable.h"

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
/// CDASolver implementing a Lagrangian-based primal-proximal heuristic
/** The PrimalProximalHeur class implements the CDASolver interface within
 * the SMS++ framework as a Lagrangian-based primal-proximal heuristic for
 * problems with binary variables, as described in
 *
 *   A. Daniilidis, C. Lemarechal "On a primal-proximal heuristic in
 *   discrete optimization" Mathematical Programming 104, 105-128, 2005
 *
 * The class is derived from LagrangianDualSolver since the heuristic is
 * built around it: at every iteration the (linear) Lagrangian objective
 * is augmented with a quadratic proximal term that pushes the convexified
 * solution produced by the Lagrangian dual towards an integer (hence
 * feasible) point; the quadratic term is then linearized exploiting the
 * 0-1 nature of the variables, so the end result is a two-level update of
 * the objective coefficients that can be passed verbatim to the inner
 * LagrangianDualSolver. Feasible integer solutions produced along the
 * way are accumulated and made available through new_var_solution().
 *
 * PrimalProximalHeur exposes very few specific parameters of its own;
 * everything else (in particular all the parameters of the inner Solver
 * solving the Lagrangian Dual) is inherited from LagrangianDualSolver. */

class PrimalProximalHeur : public LagrangianDualSolver
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
 using p_LF  = LinearFunction *;

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that PrimalProximalHeur has in addition to these of
  * LagrangianDualSolver. The value intLastPPHPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum int_par_type_PPH {

  intMaxIterPPH = intLastLDSlvPar ,  ///< maximum number of PPH iterations

  intLastPPHPar   ///< first allowed new int parameter for derived classes
                  /**< Convenience value for easily allow derived classes
                   * to extend the set of int algorithmic parameters. */

  };  // end( int_par_type_PPH )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PrimalProximalHeur has in addition to these of
  * LagrangianDualSolver. The value dblLastPPHPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum dbl_par_type_PPH {

  dbl_penaltyFactor = dblLastLDSlvPar ,  ///< PPH quadratic penalty factor

  dbl_LDSRelAcc ,  ///< relative accuracy required to the inner Solver
                   /**< The accuracy the Lagrangian Dual is solved with,
                    * i.e., what is passed to the inner Solver as its
                    * dblRelAcc. It is a different thing from the dblRelAcc
                    * of PrimalProximalHeur, which is the accuracy required
                    * of the solution the heuristic produces, cf. compute()
                    * and gap_closed(). */

  dblLastPPHPar   ///< first allowed new double parameter for derived classes
                  /**< Convenience value for easily allow derived classes
                   * to extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_PPH )

/** @} ---------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING PrimalProximalHeur ------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PrimalProximalHeur
 *  @{ */

 /// constructor: ensure every field is initialized

 PrimalProximalHeur( void ) : LagrangianDualSolver() {
  logVerb  = get_dflt_int_par( intLogVerb );
  maxIter  = get_dflt_int_par( intMaxIterPPH );
  f_MaxSol = get_dflt_int_par( intMaxSol );
  R        = get_dflt_dbl_par( dbl_penaltyFactor );
  RelAcc   = Solver::get_dflt_dbl_par( dblRelAcc );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: clean up via guts_of_destructor()

 virtual ~PrimalProximalHeur() {
  // ~LagrangianDualSolver() will take care of the base-class state;
  // here we only need to release the PrimalProximalHeur-specific resources
  guts_of_destructor();
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the PrimalProximalHeur has to solve
 /** Set (or change) the Block that the PrimalProximalHeur has to solve. The
  * method delegates to LagrangianDualSolver::set_Block() and then scans
  * the inner sub-Blocks to build the per-sub-Block dictionaries needed by
  * the proximal step (see initialize() for the details). */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// build the per-sub-Block dictionaries used by the proximal step
 /** Scan the sub-Blocks of f_Block and build, for each, the dictionaries
  * mapping the index of each binary "active" Variable in the inner
  * objective to (a) its linear coefficient and (b) its quadratic
  * coefficient, plus a count of those variables. These dictionaries are
  * used by add_penalty_terms() / remove_penalty_terms() to apply and
  * later strip the quadratic proximal term. Throws if no binary Variable
  * is found. */

 void initialize( void );

/*--------------------------------------------------------------------------*/
 /// set the int parameters of PrimalProximalHeur
 /** Set the int parameters specific of PrimalProximalHeur, and call the
  * version of LagrangianDualSolver for setting its own (which allows to
  * also set those of the inner Solver used to solve the Lagrangian Dual,
  * see the comments to set_ComputeConfig() for details). The PPH-specific
  * int parameters are:
  *
  * - intMaxIterPPH [Inf< int >()]: maximum number of PPH iterations.
  *
  * - intMaxSol [1]: maximum number of feasible solutions generated by the
  *                  heuristic that are kept; the best ones (in terms of
  *                  objective value) are kept.
  *
  * - intLogVerb: masks the first two bits of \p value ( & 3 ) and sets the
  *               verbosity of PrimalProximalHeur as
  *               = 0 : no log;
  *               = 1 : detailed iteration-by-iteration log;
  *               = 2 : even more detailed debug log;
  *               then passes \p value >> 2 (shifted right 2 places, i.e.,
  *               killing the first two bits) to LagrangianDualSolver. */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of PrimalProximalHeur
 /** Set the double parameters specific of PrimalProximalHeur, and call the
  * version of LagrangianDualSolver for setting its own (which allows to
  * also set those of the inner Solver used to solve the Lagrangian Dual,
  * see the comments to set_ComputeConfig() for details). The PPH-specific
  * double parameters are:
  *
  * - dbl_penaltyFactor [0]: factor for the quadratic proximal term added
  *                          to the inner objective to force integrality. */

 void set_par( idx_type par , double value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the PrimalProximalHeur log

 void set_log( std::ostream * log_stream = nullptr ) override {
  f_log = log_stream;
  InnerSolver->set_log( f_log );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Running the primal-proximal heuristic
 *  @{ */

 /// run the primal-proximal heuristic on the registered Block

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// lower bound: the best feasible value (max) or the dual bound (min)
 /** For a minimization problem this is the bound of the Lagrangian Dual of
  * the original objective: that of the inner Solver when no proximal
  * penalty is ever applied, i.e., when R == 0 or the Block has no static
  * binary Variable to penalise (in both cases PrimalProximalHeur
  * degenerates into a warm-started LagrangianDualSolver), and the bound of
  * the unpenalized first iteration otherwise, since the inner Solver then
  * bounds the penalized function. */

 OFValue get_lb( void ) override {
  if( f_max )
   return( best_bound );
  return( ( R == 0 ) || ( NumStatVar == 0 ) ?
	  LagrangianDualSolver::get_lb() : valid_bound );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override {
  if( ! f_max )
   return( best_bound );
  return( ( R == 0 ) || ( NumStatVar == 0 ) ?
	  LagrangianDualSolver::get_ub() : valid_bound );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_var_solution( void ) override {
  return( ! v_best_sol.empty() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_var_solution( void ) override {
  if( v_best_sol.empty() )  // no current solution
   return( false );         // no next solution either

  // v_best_sol is a heap with the best Solution at the top, i.e., the one
  // with largest value if f_max and smallest otherwise; use a custom
  // comparator (returns true if a is worse than b in the chosen sense)
  std::pop_heap( v_best_sol.begin() , v_best_sol.end() ,
                 [ this ]( const sol_value & a , const sol_value & b ) {
                  return( f_max ? ( a.second < b.second )
                                : ( a.second > b.second ) );
                  } );
  delete v_best_sol.back().first;  // delete current best solution
  v_best_sol.pop_back();           // remove current best solution
  if( v_best_sol.empty() ) {
   best_bound  =   f_max ? - Inf< double >() : Inf< double >();
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
  return( v_best_sol.front().second );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( v_best_sol.empty() )
   throw( std::logic_error(
                "PrimalProximalHeur::get_var_solution: no Solution stored"
                          ) );

  v_best_sol.front().first->write( f_Block );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the best solution found is within dblRelAcc of the valid bound
 /** The optimum lies between the bound of the unpenalized iteration and the
  * value of the best feasible solution found, so this is the criterion by
  * which the heuristic has delivered the accuracy it was required to. */

 bool gap_closed( void ) const;

/*--------------------------------------------------------------------------*/
 /// evaluate the original (unpenalized) objective at the current point
 /** Compute and return the value of the original objective function of
  * f_Block at the current values of its sub-Block variables, ignoring
  * any active proximal term. Returns +/- Inf< double >() (depending on
  * f_max) if the inner Solver does not currently have a variable
  * Solution. */

 OFValue get_funct_value( void ) {
  InnerSolver->get_var_solution();
  if( ! InnerSolver->has_var_solution() )
   return( f_max ? - Inf< double >() : Inf< double >() );

  double value = 0;
  Index idx = 0;
  for( const auto & sbi : f_Block->get_nested_Blocks() ) {
   if( is_linear[ idx ] ) {
    Funct_sbi[ idx ].compute( true );
    value += Funct_sbi[ idx ].get_value();
    }
   else {
    Funct_sbi_quad[ idx ].compute( true );
    value += Funct_sbi_quad[ idx ].get_value();
    }
   ++idx;
   }
  return( value );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of PrimalProximalHeur
 *
 * Handle the few parameters of PrimalProximalHeur and redirect everything
 * else to LagrangianDualSolver. Recall that the latter allows to change
 * via its standard parameter interface all the parameters of the "inner
 * Solver" used to actually solve the Lagrangian Dual, but the somewhat
 * convoluted mechanism is completely transparent to PrimalProximalHeur.
 *
 *  @{ */

 [[nodiscard]] idx_type int_par_first_is( void ) const override {
  return( intLastPPHPar );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_first_is( void ) const override {
  return( dblLastPPHPar );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  switch( par ) {
   case( intMaxIterPPH ): return( Inf< int >() );
   case( intLogVerb ):    return( 2 );
   }
  return( LagrangianDualSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  switch( par ) {
   case( dbl_penaltyFactor ): return( 0 );
   case( dbl_LDSRelAcc ):
    return( InnerSolver->get_dflt_dbl_par( dblRelAcc ) );
   }
  return( LagrangianDualSolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( intMaxSol ):     return( f_MaxSol );
   case( intLogVerb ):    return( logVerb +
                                  ( LagrangianDualSolver::get_int_par( par )
                                    << 4 ) );
   case( intMaxIterPPH ): return( maxIter );
   }
  return( LagrangianDualSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  switch( par ) {
   case( dbl_penaltyFactor ): return( R );
   case( dblRelAcc ):         return( RelAcc );
   case( dbl_LDSRelAcc ):     return( InnerSolver->get_dbl_par( dblRelAcc ) );
   }
  return( LagrangianDualSolver::get_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intMaxIterPPH" )
   return( intMaxIterPPH );
  return( LagrangianDualSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  if( name == "dbl_penaltyFactor" )
   return( dbl_penaltyFactor );
  if( name == "dbl_LDSRelAcc" )
   return( dbl_LDSRelAcc );
  return( LagrangianDualSolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string _ret = "intMaxIterPPH";
  if( idx == intMaxIterPPH )
   return( _ret );
  return( LagrangianDualSolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {
  static const std::string _pf = "dbl_penaltyFactor";
  static const std::string _ra = "dbl_LDSRelAcc";
  if( idx == dbl_penaltyFactor )
   return( _pf );
  if( idx == dbl_LDSRelAcc )
   return( _ra );
  return( LagrangianDualSolver::dbl_par_idx2str( idx ) );
  }

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
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters- - - - - - - - - - - - - - - - - - - - - - - - - -

 int logVerb;
 ///< verbosity of PrimalProximalHeur (verbosity of InnerSolver is + 2)

 Index maxIter;   ///< maximum number of PPH iterations

 Index f_MaxSol;  ///< maximum number of Solution to keep

 double R;        ///< proximal penalty factor (dbl_penaltyFactor)

 double RelAcc;   ///< accuracy required of the solution found (dblRelAcc)

 // working state - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double penalty  = 0;  ///< current value of the quadratic penalty term

 double addterm;       ///< current linearized penalty additive contribution

 double value_FUNCTION;
 ///< original objective value at the current candidate solution

 double best_bound;   ///< best objective value over the feasible solutions

 double valid_bound;
 ///< bound on the original problem, from the unpenalized first iteration
 /**< The value of the Lagrangian Dual of the *original* objective, which
  * the first iteration computes before any proximal penalty is applied:
  * a valid lower bound for a minimization problem (upper for a
  * maximization one), whereas the bound of every later iteration is on
  * the penalized objective and says nothing about the original one. */

 double worst_bound;  ///< worst objective value over the feasible solutions

 bool changed_penalties = false;
 ///< true if the penalty terms have been changed since the last clear

 Index NumStatVar;    ///< (current) number of binary static Variables

 Index pos_id;        ///< running counter of binary static Variables per
                      ///< sub-Block (working variable during initialize())

 // per-sub-Block dictionaries- - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< int > pos_id_sbi;
 ///< number of binary static Variables in each sub-Block

 std::vector< double > previous_sol;
 ///< value of each binary static Variable at the previous PPH iteration

 std::vector< LinearFunction > Funct_sbi;
 ///< copy of the (linear) inner objective Function of each sub-Block,
 ///< used to evaluate the unpenalized objective via get_funct_value()

 std::vector< DQuadFunction > Funct_sbi_quad;
 ///< copy of the (quadratic) inner objective Function of each sub-Block

 std::vector< bool > is_linear;
 ///< flags telling whether each sub-Block's inner objective is linear
 ///< (true) or quadratic (false)

 std::vector< std::vector< double_var > > idx_to_var_sbi1;
 ///< linear coefficient and pointer for each binary static Variable in
 ///< each sub-Block

 std::vector< std::vector< double_var > > idx_to_var_sbi2;
 ///< quadratic coefficient and pointer for each binary static Variable
 ///< in each sub-Block

 std::vector< sol_value > v_best_sol;
 ///< best feasible solutions found, managed as a binary heap

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
 /// instantiate an auxiliary Solver from a BlockSolverConfig .txt file
 /** Load the BlockSolverConfig from the given file, a sibling of this class'
  * source file, instantiate the (first) named Solver via the factory and
  * apply its ComputeConfig (if any). The Solver is NOT registered on
  * f_Block, so it cannot interfere with PrimalProximalHeur itself (which is
  * already attached to f_Block); the caller owns it and must delete it. */

 CDASolver * new_aux_solver( const std::string & cfgname );

/*--------------------------------------------------------------------------*/
 /// recover a feasible completion of the current (integer) point
 /** A Lagrangian point in general violates the coupling constraints of
  * f_Block (the ones the Lagrangian Dual dualizes), so its objective value
  * is not a valid bound. This method fixes the binary Variables driven by
  * the proximal term at their current (rounded) values and solves the
  * restricted problem on f_Block with the auxiliary Solver picked from
  * RecoveryCfg.txt, which enforces the coupling constraints: its optimum is
  * a genuinely feasible completion of the integer point. On success returns
  * true, writes the completion cost in \p cost and leaves the completion in
  * the Variables of f_Block (so that a Solution snapshot picks it up); if
  * the restricted problem is infeasible (the integer point admits no
  * feasible completion) returns false. The fixing is done with eNoMod and
  * undone before returning, so it is invisible to any other Solver. */

 bool recover_primal( double & cost );

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
