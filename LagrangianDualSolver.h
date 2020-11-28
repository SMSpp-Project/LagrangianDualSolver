/*--------------------------------------------------------------------------*/
/*---------------------- File LagrangianDualSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the LagrangianDualSolver class, which implements the
 * CDASolver interface within the SMS++ framework for a "generic"
 * Lagrangian-based Solver.
 *
 * \version 0.01
 *
 * \date 11 - 11 - 2020
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * Copyright &copy by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LagrangianDualSolver
 #define __LagrangianDualSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class AbstractBlock;  // forward declaration of AbstractBlock

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup LagrangianDualSolver_CLASSES Classes in LagrangianDualSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BundleSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver soving the Lagrangian Dual of a "generic" Block
/** The LagrangianDualSolver class implements the CDASolver interface within
 * the SMS++ framework for a "generic" Lagrangian-based Solver.
 *
 * This can "solve" (see below for the reason of the scare quotes) any Block
 * (B) with the following structure:
 *
 * - no Variable in (B)
 *
 * - (B) does not depend on any "external" Variable, i.e., a Variable that
 *   does not belong to (B) (or any of its sub-Block, recurively)
 *
 * - (B) has at least one sub-Block (necessarily, for otherwise it would be
 *   "completely empty")
 *
 * - if there is more than one sub-Block, the Constraint in (B) are all and
 *   only the ones that link its sub-Block; that is, no sub-Block must depend
 *   on any "external" Variable, i.e., a Variable that does not belong to the
 *   sub-Block (or any of its sub-Block, recurively)
 *
 * - All the Constraint in (B) are "linear constraint", i.e., FRowConstraint
 *   with a LinearFunction inside. Note that OneVarConstraint are "linear
 *   constraint" as well, but since they only concern one variable they
 *   cannot be "linking constraints". Although it may in principle be that
 *   one may want to deal with them in a Lagrangian fashion, in most of the
 *   cases including them in the subproblem is better, and therefore
 *   LagrangianDualSolver currently do not support them (although this may
 *   change later if a serious use case arises).
 *
 * - Each sub-Block may never make any assumption on which type (B) is or
 *   make any direct reference to any of its data.
 *
 * The reason for the last requirement is that LagrangianDualSolver "cheats"
 * on (B): it stealthily constructs a new Block corresponding to its
 * Lagrangian Dual, "physically moving" the sub-Block of (B) inside it while
 * not changing the pointers in (B). That is, the sub-Block of (B)
 * (temporarily) change father Block to a new Block that remains hidden 
 * inside the LagrangianDualSolver (this is undone when the
 * LagrangianDualSolver is unregistered from (B)), while (B) still
 * "believes" that they remain its sub-Block. For this, consistency is kept:
 * any Modification coming from the sub-Block is forwarded to (B).
 *
 * An appropriate Solver is then registered to the Lagrangian Dual Block, and
 * it is used to solve it. The solution it used as the dual solution for (B),
 * while a primal solution is constructed by convexification. Here comes the
 * reason for the scare quotes: if (B) does not represent a convex program
 * (say, the sub-Block have integer variables), then the Lagrangian Dual
 * Block is not equivalent to (B) but to its "convexified relaxation", and
 * this is what is solved.
 *
 * A different issue is that (B) may represent a convex program which is
 * "nonlinear enough" so that strong duality does not hold; say, the primal
 * problem may not have finite optimum (and not be unbounded), or the dual
 * problem may be infeasible even if the primal does have an optimal solution.
 * We assume that these cases either do not occur or are dealt with by the
 * user of LagrangianDualSolver. */

 class LagrangianDualSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from Function and C05Function.
 *
 *  @{ */

 using Index = Block::Index;
 using c_Index = Block::c_Index;

 using Range = Block::Range;
 using c_Range = Block::c_Range;

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;

/*----------------------------- CONSTANTS ----------------------------------*/

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "int" type that LagrangianDualSolver has in addition to these of
  * CDASolver. The value intLastLDSSlvPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum int_par_type_LDSlv {

 intLPar1 = CDASolver::intLastParCDAS ,
 ///< if the R3Block has be used for the father block

 intLastLDSSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_LDSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "double" type that LagrangianDualSolver has in addition to these of
  * CDASolver. The value dblLastLDSSlvPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum dbl_par_type_LDSlv {
  dblLastLDSlvPar = dblLastParCDAS ,
   ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set
   * of double algorithmic parameters. */

  };  // end( dbl_par_type_LDSlv )

/*@} -----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING LagrangianDualSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing LagrangianDualSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 LagrangianDualSolver( void ) : CDASolver() , NumVar( 0 ) ,
  f_LDBConfig( nullptr ) , InnrSlv( nullptr )
 {
  // ensure all parameters are properly given their default value
  LogVerb = CDASolver::get_dflt_int_par( intLogVerb );
  LPar1 = dflt_int_par[ intLPar1 - intLastParCDAS ];
  }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the LagrangianDualSolver from the Block

 virtual ~LagrangianDualSolver() { set_Block( nullptr ); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the "int" paramaters of LagrangianDualSolver
 /** Set the "int" paramaters specific of LagrangianDualSolver, together with
  * the paramaters of CDASolver that LagrangianDualSolver actually "listens to":
  *
  * - intLPar1 [0]: if the R3Block has be used for the father block
  */

 void set_par( const idx_type par , const int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the "double" paramaters of LagrangianDualSolver
 /** Set the "double" paramaters specific of LagrangianDualSolver, together
  * with the paramaters of CDASolver that LagrangianDualSolver actually
  * "listens to":
  *
  */

 void set_par( const idx_type par , const double value ) override;

/*--------------------------------------------------------------------------*/
 /// set the whole set of parameters in one blow
 /** This method sets the whole set of parameters in one blow using a
  * ComputeConfig object.
  *
  * The method of LagrangianDualSolver calls
  * ThinComputeInterface::set_ComputeConfig() to do the bulk of the work, and
  * then manages f_extra_Configuration. The field must be non-nullptr at
  * least once before compute() and any follow up-method is called, and it
  * can contain:
  *
  * - a BlockSolverConfig *
  *
  * - a SimpleConfiguration< std::pair< Configuration * , Configuration * > >
  *   there .first is a BlockSolverConfig * and .second is a
  *   BlockConfig *
  *
  * The BlockSolverConfig is apply()-ed to the Lagrangian Dual Block; it
  * has to register *at least* one (appropriate) CDASolver to it, and it
  * also has the chance to register Solver to the inner Block (inside the
  * LagBFunction inside the FRealObjective of each sub-Block). The *first*
  * CDASolver registered to the Lagrangian Dual Block will be used as the
  * "main" Solver by LagrangianDualSolver, and there typically must be at
  * least one Solver registered to each sub-Block (unless the "main" Solver
  * allows not to because it deals with some LagBFunction in specialised
  * ways), the first of which will be used by the LagBFunction to compute()
  * themselves.
  *
  * If a BlockConfig is provided, it is also apply()-ed to the Lagrangian
  * Dual Block (and therefore it has the chance to BlockConfig-ure also
  * all the inner Block inside the LagBFunction inside the FRealObjective of
  * each sub-Block).
  *
  * The BlockSolverConfig becomes "property" of the LagrangianDualSolver;
  * it is clear()-ed and used when the LagrangianDualSolver is unregistered
  * from the Block to clear away all the Solver that it has attached.
  *
  * Note that, if the LagrangianDualSolver is not registered to a Block
  * when set_ComputeConfig() is called, the pointers are stored and the
  * configuration of the Block is done in set_Block(). For the
  * BlockSolverConfig, this could be repeated any number of times since a
  * BlockSolverConfig is not "consumed" when it is apply()-ed. However, a
  * BlockConfig is "consumed" instead. Thus, for consistency also the
  * BlockSolverConfig is immediately clear()-ed: registering the
  * LagrangianDualSolver to a new Block requires a new call to
  * set_ComputeConfig().
  *
  * Note that if multiple calls to set_ComputeConfig() occur before the
  * LagrangianDualSolver is registered to a Block, all the BlockSolverConfig
  * / BlockConfig of all calls save the last one are lost and have no effect
  * (but at least are properly deleted). */

 void set_ComputeConfig( ComputeConfig * scfg = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the LagrangianDualSolver log

 void set_log( std::ostream *log_stream = nullptr ) override {
  if( f_log != log_stream ) {
   f_log = log_stream;
   if( InnrSlv )
    InnrSlv->set_log( f_log );
   }
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Lagrangian Dual of the given Block
 *  @{ */

 /// (try to) solve the Lagrangian Dual of the given Block

 int compute( bool changedvars = true ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 OFValue get_lb( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->get_ub() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->get_lb() );
  }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->has_dual_solution() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_var_solution( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->has_dual_solution() );
  }

/*--------------------------------------------------------------------------*/

 bool is_var_feasible( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->is_dual_feasible() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool is_dual_feasible( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->is_var_feasible() );
  }

/*--------------------------------------------------------------------------*/
 /// write the "current" solution

 void get_var_solution( Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution

 void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->new_dual_solution() )
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_dual_solution( void ) override {
  #ifndef NDEBUG
   if( ! InnrSlv )
    throw( std::logic_error( "inner CDASolver not initialised yet" ) );
  #endif
  return( InnrSlv->new_var_solution() );
  }

/*--------------------------------------------------------------------------*/
/*
  void set_unbounded_threshold( const VarValue thr ) override { }

  bool has_var_direction( void ) override { return( true ); }

  bool has_dual_direction( void ) override { return( true ); }

  void get_var_direction( Configuration *dirc = nullptr ) override {}

  void get_dual_direction( Configuration *dirc = nullptr ) override {}

  virtual bool new_var_direction( void ) override { return( false ); }
  
  virtual bool new_dual_direction( void ) override{ return( false ); }
*/

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the LagrangianDualSolver
 *
 *  @{ */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastLdsSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastLdsSlvPar ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_dflt_int_par( const idx_type par ) const override {
  if( ( par >= intLastParCDAS ) && ( par < intLastLdsSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );
  else
   return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( const idx_type par ) const override {
  if( ( par >= dblLastParCDAS ) && ( par < dblLastLdsSlvPar ) )
   return( dflt_dbl_par[ par - dblLastParCDAS ] );
  else
   return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( const idx_type par ) const override;

/*--------------------------------------------------------------------------*/

 idx_type int_par_str2idx( const std::string & name ) const override {
  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( const idx_type idx ) const override {
  if( ( idx >= intLastParCDAS ) && ( idx < intLastLdsSlvPar ) )
   return( int_pars_str[ idx - intLPar1 ] );
  else
   return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( const idx_type idx ) const override {
  if( ( idx >= dblLastParCDAS ) && ( idx < dblLastLdsSlvPar ) )
   return( dbl_pars_str[ idx - dblLastParCDAS ] );
  else
   return( CDASolver::dbl_par_idx2str( idx ) );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void Log1( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void Log2( void );

/*--------------------------------------------------------------------------*/

 void set_default_inner_BlockSolverConfig( void );

/*--------------------------------------------------------------------------*/

 void configure_LagrangianDualBlock( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorthmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 int LogVerb;       ///< "verbosity" of the log

 int LPar1;         ///< if the R3Block conversion has to be done

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumVar;      ///< (current) number of variables

 //!! std::vector< ColVariable * > LamVcblr;  ///< map Lambda -> ColVariable

 AbstractBlock * LagrDual;

 Configuration * f_LDBConfig;   ///< the Configuration for LagrDual
 
 CDASolver * InnrSlv;

/*--------------------------------------------------------------------------*/

 const static std::vector<int> dflt_int_par;
 ///< the (static const) vector of int parameters default values

 const static std::vector<double> dflt_dbl_par;
 ///< the (static const) vector of double parameters default values

 const static std::vector< std::string > int_pars_str;
 ///< the (static const) vector of int parameters names

 const static std::vector< std::string > dbl_pars_str;
 ///< the (static const) vector of double parameters names

 const static std::map< std::string , idx_type > int_pars_map;
  ///< the (static const) map for int parameters names

 const static std::map< std::string , idx_type > dbl_pars_map;
 ///< the (static const) map for double parameters names

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

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class LagrangianDualSolver )

/*@}  end( group( LagrangianDualSolver_CLASSES ) ) -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* LagrangianDualSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File LagrangianDualSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
