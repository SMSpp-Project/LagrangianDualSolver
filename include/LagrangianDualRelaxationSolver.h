/*--------------------------------------------------------------------------*/
/*---------------- File LagrangianDualRelaxationSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class LagrangianDualRelaxationSolver, the
 * RelaxationSolver [see ChangeSolver.h] that uses the Lagrangian Dual of a
 * Block, and the primal solutions that PrimalProximalHeur recovers from it,
 * as the relaxation of the nodes of a Branch-and-Bound (as in
 * BranchAndXSolver), and for the class LagrangianChange, the Change that it
 * applies at the nodes.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LagrangianDualRelaxationSolver
 #define __LagrangianDualRelaxationSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include <unordered_map>

#include "AbstractChange.h"

#include "ChangeSolver.h"

#include "GlobalInformation.h"

#include "PrimalProximalHeur.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS LagrangianChange --------------------------*/
/*--------------------------------------------------------------------------*/
/// the AbstractChange applied at the nodes by LagrangianDualRelaxationSolver
/** A LagrangianChange is an AbstractChange whose paths are relative to the
 * LagrangianDualBlock of a LagrangianDualRelaxationSolver, which applies it
 * itself [see LagrangianDualRelaxationSolver::apply()]; applied to a Block,
 * this has to be an AbstractBlock. */

class LagrangianChange : public AbstractChange {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: an empty change

 LagrangianChange( void ) : AbstractChange() {}

 /// constructor: the type, the data and the paths of the change

 LagrangianChange( int type , std::vector< double > data ,
		   std::vector< AbstractPath > paths )
  : AbstractChange( type , std::move( data ) , std::move( paths ) ) {}

 /// destructor: does nothing

 ~LagrangianChange() override = default;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// apply the LagrangianChange to the given AbstractBlock

 Change * apply( Block * block , bool doUndo = false ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck ) override {
  if( ! dynamic_cast< AbstractBlock * >( block ) )
   throw( std::invalid_argument( "LagrangianChange::apply: the Block is "
				 "not an AbstractBlock" ) );
  return( AbstractChange::apply( block , doUndo , issueMod , issueAMod ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  // insert LagrangianChange in the factory

 };  // end( class( LagrangianChange ) )

/*--------------------------------------------------------------------------*/
/*---------------- CLASS LagrangianDualRelaxationSolver --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Lagrangian Dual as the relaxation of a Branch-and-Bound
/** LagrangianDualRelaxationSolver is the RelaxationSolver that solves at
 * each node of a Branch-and-Bound the Lagrangian Dual of the Block, with
 * the primal recovery of PrimalProximalHeur, which also gives the true
 * bounds and solution of the node. branch() picks a variable whose value in
 * the Lagrangian solution (the convexified one of the Lagrangian Dual, see
 * PrimalProximalHeur::get_Lagrangian_initial_solution()) is fractional [see
 * intBranchStrategy], and produces the two LagrangianChange that branch on
 * it, which apply() applies [see intApplyStrategy]:
 *
 * The LagrangianChange identify the variable by its AbstractPath relative
 * to the LagrangianDualBlock if the sub-Blocks are copied into it
 * (int_LDSlv_iBCopy), and to the Block being solved otherwise.
 *
 * - eMaster: in the Lagrangian Dual, i.e., the bound of each child (x <=
 *   floor( v ) or x >= ceil( v )) becomes a new dual pair of the
 *   LagBFunction of the sub-Block of the variable, a later bound on the
 *   same variable changing the constant term of its pair;
 *
 * - eSubproblem: in the sub-Block, fixing the variable to floor( v ) or to
 *   ceil( v ); the columns of the global pool of the LagBFunction that the
 *   fixing makes infeasible are kept in the GlobalInformation [see
 *   set_global_information()], and given back to the LagBFunction when the
 *   variable is unfixed. */

class LagrangianDualRelaxationSolver : public RelaxationSolver ,
				       public PrimalProximalHeur {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Block::Index;

 /// name of the Collection of the GlobalInformation with the purged columns

 static constexpr const char * str_VarToSol = "map_varToSol";

 /// key, in that Collection, of the purged columns

 static constexpr const char * str_PurgedColumns = "purgedColumns";

 /// the purged columns, by the variable whose fixing purged them

 using PurgedColumn = std::map< ColVariable * ,
				std::vector< LagBFunction::gpool_el > >;

 /// the possible values of intBranchStrategy

 enum branch_strategy {
  eMostFractional = 0  ///< the variable with the most fractional value
  };

 /// the possible values of intApplyStrategy

 enum apply_strategy {
  eMaster = 0 ,     ///< the bounds become dual pairs of the LagBFunction
  eSubproblem = 1   ///< the variable is fixed in the sub-Block
  };

 /// public enum for the int algorithmic parameters
 /** Public enum describing the algorithmic parameters of int type that
  * LagrangianDualRelaxationSolver has in addition to these of
  * PrimalProximalHeur:
  *
  * - intApplyStrategy [eMaster]: where the branching is applied [see
  *   apply_strategy and the general notes of the class];
  *
  * - intBranchStrategy [eMostFractional]: how the branching variable is
  *   chosen [see branch_strategy]. */

 enum int_par_type_LDRS {
  intApplyStrategy = intLastPPHPar ,  ///< where the branching is applied
  intBranchStrategy ,                 ///< how the branching variable is chosen
  intLastLDRSPar      ///< first allowed new int parameter for derived classes
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor

 LagrangianDualRelaxationSolver( void ) : PrimalProximalHeur() {}

 /// destructor

 ~LagrangianDualRelaxationSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 using PrimalProximalHeur::set_par;  // keep the other set_par() visible

 void set_par( idx_type par , int value ) override;

 [[nodiscard]] idx_type int_par_first_is( void ) const override {
  return( intLastLDRSPar );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// give the Solver the GlobalInformation of the search
 /** Besides recording it, makes sure that the GlobalInformation has the
  * Collection where the columns purged by the fixings are kept [see
  * str_VarToSol and str_PurgedColumns]. */

 void set_global_information( GlobalInformation * gi ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the Lagrangian Dual of the node, with the primal recovery
 /** As PrimalProximalHeur::compute(), but a status between kOK and kError
  * (or kLowPrecision) is kOK when the Lagrangian solution is there, as
  * branch() only needs that. */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override { return( PrimalProximalHeur::get_lb() ); }

 OFValue get_ub( void ) override { return( PrimalProximalHeur::get_ub() ); }

 OFValue get_true_lb( void ) override {
  return( PrimalProximalHeur::get_lb() );
  }

 OFValue get_true_ub( void ) override {
  return( PrimalProximalHeur::get_ub() );
  }

 bool has_true_var_solution( void ) override {
  return( PrimalProximalHeur::has_var_solution() );
  }

 bool new_true_var_solution( void ) override {
  return( PrimalProximalHeur::new_var_solution() );
  }

 /// the true solution: the primal solution of PrimalProximalHeur

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  PrimalProximalHeur::get_var_solution( solc );
  }

 /// the Solution of the relaxation: that of the Lagrangian Dual

 Solution * get_Solution( Configuration * solc = nullptr ) override;

 /// the true Solution: that of PrimalProximalHeur

 Solution * get_true_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// branch on a fractional variable [see intBranchStrategy]

 std::vector< Change * > branch( void ) override;

 /// apply a LagrangianChange [see intApplyStrategy]

 Change * apply( Change * change , bool doUndo = false ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the Block the paths of the LagrangianChange are relative to
 /** The variables of the sub-Blocks are in the LagrangianDualBlock when
  * these are copied into it (int_LDSlv_iBCopy), and in the Block being
  * solved otherwise. */

 Block * path_base( void ) const {
  return( iBCopy ? static_cast< Block * >( LagrDual ) : f_Block );
  }

/*--------------------------------------------------------------------------*/
 /// the LagBFunction of the sub-Block the variable belongs to

 LagBFunction * LagBF_of( Variable * var ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int f_branch_strategy = eMostFractional;  ///< intBranchStrategy
 int f_apply_strategy = eMaster;           ///< intApplyStrategy

 /// for each variable, the dual pairs of its upper (first) and lower
 /// (second) bound added by apply() with eMaster, nullptr if none

 std::unordered_map< Variable * ,
		     std::pair< LinearFunction * , LinearFunction * > >
  map_varToLF;

 /// the Collection of the GlobalInformation with the purged columns

 std::shared_ptr< Collection< PurgedColumn > > map_varToSol;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( LagrangianDualRelaxationSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* LagrangianDualRelaxationSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File LagrangianDualRelaxationSolver.h ---------------*/
/*--------------------------------------------------------------------------*/
