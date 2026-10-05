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

#include <algorithm>

#include <array>

#include <map>
#include <memory>
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
 /// the multipliers of the node that branched, which a child may start from

 using Center = std::vector< std::pair< ColVariable * , double > >;

 /// set the multipliers of the node that branched [see get_center()]

 void set_center( std::shared_ptr< const Center > center ) {
  f_center = std::move( center );
  }

 /// the multipliers of the node that branched, nullptr if none
 /** Returns the Lagrangian multipliers, each with its value, that the node
  * producing this LagrangianChange had when it branched, i.e., the stability
  * centre its Lagrangian Dual ended at, which the child starts from [see
  * LagrangianDualRelaxationSolver::apply()]. The LagrangianChange that undo
  * another one carry none. */

 [[nodiscard]] const Center * get_center( void ) const {
  return( f_center.get() );
  }

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

 std::shared_ptr< const Center > f_center;  ///< see get_center()

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
 * bounds and solution of the node. branch() picks a variable whose value
 * in the Lagrangian solution (the convexified one of the Lagrangian Dual,
 * see PrimalProximalHeur::get_Lagrangian_convexified_solution()) is
 * fractional [see intBranchStrategy], and produces the two LagrangianChange
 * that branch on it, which apply() applies [see intApplyStrategy]. The
 * Lagrangian Dual has to be solved accurately enough (see dblNZEps and
 * intWZNorm of the inner Solver) that an integral Lagrangian solution
 * satisfies the relaxed constraints, and hence solves the node: branch()
 * finding none fractional throws.
 *
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
  eMostFractional = 0 ,  ///< the variable with the most fractional value
  eStrongBranching = 1 , ///< the best of intStrongCands by strong branching
  eLearned = 2 ,         ///< the best of intStrongCands by strBranchModel
  eOnline = 3            ///< strong branching, then a model learned online
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
  *   chosen [see branch_strategy]: with eStrongBranching the intStrongCands
  *   most fractional variables are candidates, and for each of them the
  *   Lagrangian Dual of both children is solved (without the heuristic),
  *   the chosen one being that with the largest product of the bound
  *   improvements of its children over the node, a child found infeasible
  *   counting as an infinite improvement;
  *
  * - intStrongCands [10]: the number of candidates of eStrongBranching and
  *   of eLearned, which picks the one a model learned from the data of
  *   strong branching ranks first [see strBranchModel]; eOnline learns that
  *   model while the search goes on, and is only handled by
  *   LagrangianDualRelaxationSolverML, which needs Torch. */

 enum int_par_type_LDRS {
  intApplyStrategy = intLastPPHPar ,  ///< where the branching is applied
  intBranchStrategy ,                 ///< how the branching variable is chosen
  intStrongCands ,                    ///< candidates of strong branching
  intLastLDRSPar      ///< first allowed new int parameter for derived classes
  };

 /// public enum for the string algorithmic parameters
 /** Public enum describing the algorithmic parameters of string type that
  * LagrangianDualRelaxationSolver has in addition to these of
  * PrimalProximalHeur:
  *
  * - strStrongLog [""]: the file to which eStrongBranching appends one line
  *   per candidate, in CSV, so that a branching rule can be learned from
  *   strong branching: the progressive number of the call of branch(), the
  *   bound of the node, the fraction of the binary variables that are fixed,
  *   the rank of the candidate, the index of its sub-Block and its position
  *   in it (in [0, 1]), its value in the convexified and in the Lagrangian
  *   solution, its fractionality, its cost in the objective of the sub-Block,
  *   the bounds of the two children, the score, and 1 for the chosen one; a
  *   header line is written when the file is created, and with the empty
  *   string nothing is written;
  *
  * - strBranchModel [""]: the file of the model of eLearned, a perceptron
  *   with one tanh hidden layer that gives a score to each candidate out of
  *   its features (in this order: its fractionality, its value in the
  *   convexified and in the Lagrangian solution, the absolute difference of
  *   the two, its cost over the largest absolute cost of the candidates, its
  *   position in its sub-Block, the fraction of the binary variables that
  *   are fixed, its rank by fractionality over the number of candidates
  *   less one), written as the line "LBRModel <nin> <nhid>", the nhid rows
  *   of nin weights of the hidden layer, its nhid biases, the nhid weights
  *   of the output and its bias. */

 enum str_par_type_LDRS {
  strStrongLog = strLastPPHPar ,  ///< the file of the strong branching data
  strBranchModel ,                ///< the file of the model of eLearned
  strLastLDRSPar      ///< first allowed new str parameter for derived classes
  };

 /// public enum for the vector-of-string algorithmic parameters
 /** Public enum describing the algorithmic parameters of vector-of-string
  * type that LagrangianDualRelaxationSolver has in addition to these of
  * PrimalProximalHeur:
  *
  * - vstrBranchGroups [empty]: the names of the groups of static Variable of
  *   the sub-Blocks that are branched upon, the binary variables of any
  *   other group being never candidates; empty means all of them. A group
  *   whose variables are determined by these of the others (say, start-up
  *   variables by the commitment ones) need not be branched upon, and
  *   leaving it out also spares the Solver of the sub-Block fixings it may
  *   not support. */

 enum vstr_par_type_LDRS {
  vstrBranchGroups = vstrLastLDSlvPar ,  ///< the groups branched upon
  vstrLastLDRSPar  ///< first allowed new vector-of-string parameter
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: branch() needs the convexified Lagrangian solution

 LagrangianDualRelaxationSolver( void ) : PrimalProximalHeur() {
  f_save_conv_sol = true;
  }

 /// destructor

 ~LagrangianDualRelaxationSolver() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 using PrimalProximalHeur::set_par;  // keep the other set_par() visible

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , std::string && value ) override;

 void set_par( idx_type par , std::vector< std::string > && value ) override;

 [[nodiscard]] idx_type int_par_first_is( void ) const override {
  return( intLastLDRSPar );
  }

 [[nodiscard]] idx_type str_par_first_is( void ) const override {
  return( strLastLDRSPar );
  }

 [[nodiscard]] idx_type vstr_par_first_is( void ) const override {
  return( vstrLastLDRSPar );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type par )
  const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type par )
  const override;

 [[nodiscard]] const std::vector< std::string > & get_dflt_vstr_par(
						 idx_type par ) const override;

 [[nodiscard]] const std::vector< std::string > & get_vstr_par(
						 idx_type par ) const override;

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type par )
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
  * branch() only needs that; when it is not there, the status is that of
  * the Lagrangian Dual, or kError if that is kOK. */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// the bounds of the relaxation
 /** On the side of the relaxation (the lower one for a minimization
  * problem), that of PrimalProximalHeur, i.e., the bound of the Lagrangian
  * Dual of the original objective; with a proximal penalty
  * (dbl_penaltyFactor > 0) the last call of the inner Solver bounds the
  * penalized function, which is not a bound on the node. On the other side,
  * the best between the bound of the Lagrangian Dual on that side [see
  * PrimalProximalHeur::relaxation_other_bound()], which may be infinite,
  * and the value of the best feasible solution of the heuristic, which
  * bounds the relaxation too. */

 OFValue get_lb( void ) override {
  if( ! f_max )
   return( PrimalProximalHeur::get_lb() );
  return( std::max( PrimalProximalHeur::get_lb() ,
		    relaxation_other_bound() ) );
  }

 OFValue get_ub( void ) override {
  if( f_max )
   return( PrimalProximalHeur::get_ub() );
  return( std::min( PrimalProximalHeur::get_ub() ,
		    relaxation_other_bound() ) );
  }

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

 /// the true Solution: the best one of PrimalProximalHeur
 /** A copy of the best Solution that PrimalProximalHeur has saved, which the
  * caller owns, or nullptr if there is none; only a \p solc asking for a
  * part of it requires writing it into the Block and reading it back. */

 Solution * get_true_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// branch on a fractional variable [see intBranchStrategy]

 std::vector< Change * > branch( void ) override;

 /// apply a LagrangianChange [see intApplyStrategy]

 Change * apply( Change * change , bool doUndo = false ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// a candidate of the branching, a variable fractional in the convexified
 /// solution

 struct Cand {
  ColVariable * var;  ///< the variable
  double value;       ///< its value in the convexified solution
  double frac;        ///< its fractionality
  Index k;            ///< its index among the binary static variables
  Index sb;           ///< its sub-Block
  Index pos;          ///< its position among those of the sub-Block
  double cost;        ///< its cost in the objective of the sub-Block
  };

 /// the number of the features of a candidate [see strBranchModel]

 static constexpr Index NFeatures = 8;

 /// the features of a candidate [see strBranchModel]

 using Features = std::array< double , NFeatures >;

 /// what strong branching finds of a candidate: the bound of the "down"
 /// child, that of the "up" one and the score [see intBranchStrategy]

 using StrongEval = std::array< double , 3 >;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the candidates of the node, the most fractional first
 /** The binary variables of the sub-Blocks that may be branched upon [see
  * vstrBranchGroups] and are fractional in the convexified solution; throws
  * if there is none [see branch()]. */

 std::vector< Cand > candidates( void );

/*--------------------------------------------------------------------------*/
 /// the features of the candidates [see strBranchModel]
 /** They are relative to the candidates given (say, the cost is over the
  * largest one among them, the rank over their number), which are in the
  * order of candidates(). */

 std::vector< Features > features( const std::vector< Cand > & cands ) const;

/*--------------------------------------------------------------------------*/
 /// strong branching on the candidates [see intBranchStrategy]

 std::vector< StrongEval > strong_branching(
					      const std::vector< Cand > & cands );

/*--------------------------------------------------------------------------*/
 /// append the data of strong branching to strStrongLog (if any)

 void write_strong_log( const std::vector< Cand > & cands ,
			const std::vector< StrongEval > & evals ,
			Index chosen );

/*--------------------------------------------------------------------------*/
 /// the candidate to branch upon, by intBranchStrategy
 /** Returns the index in \p cands of the chosen candidate; \p cands may be
  * shortened (say, to the intStrongCands most fractional ones). A derived
  * class that handles more values of intBranchStrategy (as
  * LagrangianDualRelaxationSolverML does with eOnline) overrides this. */

 virtual Index choose( std::vector< Cand > & cands );

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 int f_branch_strategy = eMostFractional;  ///< intBranchStrategy
 int f_strong_cands = 10;                  ///< intStrongCands
 std::string f_strong_log;                 ///< strStrongLog
 Index f_n_branch = 0;                     ///< calls of branch() so far

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
 /// the two LagrangianChange branching on var at the fractional value
 /** The first one is the "down" child (var <= floor( value ), or var fixed
  * to it), the second the "up" one [see intApplyStrategy]. */

 std::vector< Change * > branchings( ColVariable * var , double value ) const;

/*--------------------------------------------------------------------------*/
 /// the bound of the node with the LagrangianChange applied
 /** Applies the Change, solves the Lagrangian Dual alone (not the
  * heuristic), and undoes the Change; returns the bound on the side of the
  * relaxation (the lower one for a minimization problem), +/- infinity if
  * the Lagrangian Dual finds the child infeasible. */

 double child_bound( Change * change );

/*--------------------------------------------------------------------------*/
 /// the binary static variables that may be branched upon
 /** One flag per variable, in the order of the Lagrangian solution [see
  * vstrBranchGroups]. */

 std::vector< bool > branchable( void ) const;

/*--------------------------------------------------------------------------*/
 /// read the model of eLearned from strBranchModel

 void load_branch_model( void );

/*--------------------------------------------------------------------------*/
 /// the Lagrangian multipliers, i.e., the active Variable of the LagBFunction

 std::shared_ptr< const LagrangianChange::Center > center( void ) const;

/*--------------------------------------------------------------------------*/
 /// put back the multipliers of center, those added since then at 0
 /** Writes in the Lagrangian multipliers the values in center, and 0 in
  * those that are not there, i.e., the multipliers of the dual pairs added
  * since center was taken; a Variable of center that is no longer a
  * multiplier is left alone. The inner Solver is asked, once, to restart
  * from them at its next compute(), emptying its bundle [see intRstAlg of
  * BundleSolver]; returns false, doing nothing, if it cannot be asked. */

 bool restore_center( const LagrangianChange::Center & center );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int f_apply_strategy = eMaster;           ///< intApplyStrategy
 std::string f_branch_model;               ///< strBranchModel
 std::vector< std::string > f_branch_groups; ///< vstrBranchGroups

 /// the model of eLearned, read from strBranchModel at its first use: the
 /// weights and the biases of the hidden layer, those of the output
 std::vector< std::vector< double > > f_w1;
 std::vector< double > f_b1;
 std::vector< double > f_w2;
 double f_b2 = 0;

 /// the multipliers taken by the last branch(), if no node has been solved
 /// since then, so that they are still those of the inner Solver [see
 /// apply()]; nullptr otherwise
 mutable const LagrangianChange::Center * f_fresh_center = nullptr;

 /// the intRstAlg of the inner Solver to put back after the compute() that
 /// restore_center() asked a restart for, -1 if none
 int f_rst_alg = -1;

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
