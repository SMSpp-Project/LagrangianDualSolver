/*--------------------------------------------------------------------------*/
/*--------------- File LagrangianDualRelaxationSolverML.h ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class LagrangianDualRelaxationSolverML, the
 * LagrangianDualRelaxationSolver that learns its branching rule while the
 * Branch-and-Bound goes on: on the first nodes it branches by strong
 * branching, keeping the features of the candidates and their scores, and
 * then it fits on those data, with the Torch C++ API, the scoring function
 * by which it chooses on all the other nodes (intBranchStrategy eOnline).
 * This is the idea of
 *
 *  E.B. Khalil, P. Le Bodic, L. Song, G. Nemhauser, B. Dilkina "Learning to
 *  Branch in Mixed Integer Programming" Proceedings of AAAI 2016
 *
 * with the scoring function and the loss of the eLearned rule of
 * LagrangianDualRelaxationSolver [see strBranchModel], which learns the same
 * function out of a file written by a previous run.
 *
 * The class is in a library of its own, which is only built where Torch is:
 * whoever wants it links SMS++::LagrangianDualSolverML.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LagrangianDualRelaxationSolverML
 #define __LagrangianDualRelaxationSolverML
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <torch/torch.h>

#include "LagrangianDualRelaxationSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------- CLASS LagrangianDualRelaxationSolverML ------------------*/
/*--------------------------------------------------------------------------*/
/// a LagrangianDualRelaxationSolver that learns its branching rule online
/** With intBranchStrategy == eOnline, LagrangianDualRelaxationSolverML
 * branches by strong branching on the intStrongCands most fractional
 * candidates of the first intLearnNodes nodes, and keeps for each node the
 * features of its candidates [see strBranchModel] and their scores, over the
 * largest one; then it fits on those data the perceptron of eLearned (one
 * tanh hidden layer of intHidden units) by the listwise loss of
 *
 *  Z. Cao, T. Qin, T.-Y. Liu, M.-F. Tsai, H. Li "Learning to Rank: From
 *  Pairwise Approach to Listwise Approach" Proceedings of ICML 2007
 *
 * i.e., the cross entropy between the distribution that the softmax of the
 * scores of the model gives on the candidates of a node and that of ten
 * times their relative strong branching scores, by intTrainEpochs epochs of
 * Adam with step dblLearnRate, and from then on it branches on the
 * candidate the model scores best. Any other value of intBranchStrategy is
 * that of LagrangianDualRelaxationSolver.
 *
 * The data and the model are kept in the Collection str_OnlineBranching of
 * the GlobalInformation of the search [see set_global_information()], so
 * that all the LagrangianDualRelaxationSolverML of a search (say, those of
 * the workers of a parallel BranchAndXSolver) gather their data together
 * and use the same model, which the first of them to reach intLearnNodes
 * nodes fits; without a GlobalInformation each one keeps its own. Once
 * fitted, the model can be written to strOnlineModel in the format of
 * strBranchModel, so that another run can use it by eLearned. */

class LagrangianDualRelaxationSolverML :
 public LagrangianDualRelaxationSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// name, in the GlobalInformation, of the Collection of the data and of
 /// the model of eOnline

 static constexpr const char * str_OnlineBranching = "onlineBranching";

 /// key, in that Collection, of the data and of the model

 static constexpr const char * str_OnlineData = "onlineData";

 /// the data of eOnline: the features and the targets of the candidates of
 /// each node of strong branching, and the model once it is fitted

 struct OnlineData {
  std::vector< std::vector< Features > > x;  ///< features, per node
  std::vector< std::vector< double > > y;    ///< targets, per node
  torch::nn::Sequential net{ nullptr };      ///< the model, once fitted
  bool fitting = false;                      ///< whether one is fitting it
  };

 /// public enum for the int algorithmic parameters
 /** Public enum describing the algorithmic parameters of int type that
  * LagrangianDualRelaxationSolverML has in addition to these of
  * LagrangianDualRelaxationSolver:
  *
  * - intLearnNodes [20]: the number of nodes branched by strong branching
  *   before the model is fitted;
  *
  * - intTrainEpochs [300]: the epochs of the fitting;
  *
  * - intHidden [16]: the units of the hidden layer of the model;
  *
  * - intLearnSeed [0]: the seed of the random initialization of the
  *   model. */

 enum int_par_type_LDRSML {
  intLearnNodes = intLastLDRSPar ,  ///< nodes of strong branching
  intTrainEpochs ,                  ///< epochs of the fitting
  intHidden ,                       ///< hidden units of the model
  intLearnSeed ,                    ///< seed of the model
  intLastLDRSMLPar  ///< first allowed new int parameter for derived classes
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the algorithmic parameters of double type that
  * LagrangianDualRelaxationSolverML has in addition to these of
  * LagrangianDualRelaxationSolver:
  *
  * - dblLearnRate [1e-2]: the step of Adam in the fitting. */

 enum dbl_par_type_LDRSML {
  dblLearnRate = dblLastPPHPar ,    ///< step of the fitting
  dblLastLDRSMLPar  ///< first allowed new double parameter
  };

 /// public enum for the string algorithmic parameters
 /** Public enum describing the algorithmic parameters of string type that
  * LagrangianDualRelaxationSolverML has in addition to these of
  * LagrangianDualRelaxationSolver:
  *
  * - strOnlineModel [""]: the file to which the model is written once
  *   fitted, in the format of strBranchModel; with the empty string nothing
  *   is written. */

 enum str_par_type_LDRSML {
  strOnlineModel = strLastLDRSPar ,  ///< the file of the fitted model
  strLastLDRSMLPar  ///< first allowed new str parameter for derived classes
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor

 LagrangianDualRelaxationSolverML( void ) : LagrangianDualRelaxationSolver()
 {}

 /// destructor

 ~LagrangianDualRelaxationSolverML() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 using LagrangianDualRelaxationSolver::set_par;  // keep the others visible

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 void set_par( idx_type par , std::string && value ) override;

 [[nodiscard]] idx_type int_par_first_is( void ) const override {
  return( intLastLDRSMLPar );
  }

 [[nodiscard]] idx_type dbl_par_first_is( void ) const override {
  return( dblLastLDRSMLPar );
  }

 [[nodiscard]] idx_type str_par_first_is( void ) const override {
  return( strLastLDRSMLPar );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type par )
  const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type par )
  const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// give the Solver the GlobalInformation of the search
 /** Besides what LagrangianDualRelaxationSolver does, makes sure that the
  * GlobalInformation has the Collection of the data and of the model of
  * eOnline [see str_OnlineBranching]. */

 void set_global_information( GlobalInformation * gi ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
 /// the candidate to branch upon, by eOnline or as the base class does

 Index choose( std::vector< Cand > & cands ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// call f on the OnlineData, under the lock of their Collection (if any)

 template< class F >
 void with_data( F && f );

/*--------------------------------------------------------------------------*/
 /// fit the model on the data of the nodes of strong branching

 torch::nn::Sequential fit( const std::vector< std::vector< Features > > & x ,
			    const std::vector< std::vector< double > > & y );

/*--------------------------------------------------------------------------*/
 /// write the model to strOnlineModel, in the format of strBranchModel

 void write_model( torch::nn::Sequential & net ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int f_learn_nodes = 20;                ///< intLearnNodes
 int f_train_epochs = 300;              ///< intTrainEpochs
 int f_hidden = 16;                     ///< intHidden
 int f_learn_seed = 0;                  ///< intLearnSeed
 double f_learn_rate = 1e-2;            ///< dblLearnRate
 std::string f_online_model;            ///< strOnlineModel

 /// the Collection of the GlobalInformation with the data, nullptr if none

 std::shared_ptr< Collection< OnlineData > > f_online;

 /// the data when there is no GlobalInformation

 OnlineData f_local;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( LagrangianDualRelaxationSolverML ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* LagrangianDualRelaxationSolverML.h included */

/*--------------------------------------------------------------------------*/
/*--------------- End File LagrangianDualRelaxationSolverML.h --------------*/
/*--------------------------------------------------------------------------*/
