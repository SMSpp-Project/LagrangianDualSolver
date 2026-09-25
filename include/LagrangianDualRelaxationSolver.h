#pragma once

#include "PrimalProximalHeur.h"
#include "AbstractPath.h"
#include "ChangeSolver.h"

namespace SMSpp_di_unipi_it
{

    /// LagrangianChange is a change that accept as a block only a LagrangianDualBlock

    class LagrangianChange : public AbstractChange
    {
    public:
        enum LagrangianChangeType
        {
            eDeleteUB = eLastACTtype, ///< rimuove il vincolo di UB da f_Block
            eDeleteLB,                ///< rimuove il vincolo di LB da f_Block
            eLastLagrangianChangeType ///< primo valore libero per le classi derivate
        };

        LagrangianChange() : AbstractChange() {}

        LagrangianChange(int type, std::vector<double> value,
                         std::vector<AbstractPath> paths)
            : AbstractChange(static_cast<AbstractChange::AbstractChangeType>(type), std::move(value), std::move(paths)) {}

        ~LagrangianChange() = default;

        Change *apply(Block *block, bool doUndo = false,
                      ModParam issueMod = eNoBlck,
                      ModParam issueAMod = eNoBlck) override
        {
            if (!dynamic_cast<AbstractBlock *>(block))
                throw std::invalid_argument(
                    "LagrangianChange::apply: block is not an AbstractBlock");
            switch (f_type)
            {
            case eDeleteUB:
            case eDeleteLB:
                throw std::invalid_argument(
                    "LagrangianChange::apply: eDeleteUB/eDeleteLB richiedono lo stato "
                    "interno del solver e vanno applicati tramite "
                    "LagrangianDualRelaxationSolver::apply(), non genericamente");
            default:
                return AbstractChange::apply(block, doUndo, issueMod, issueAMod);
            }
        }

    private:
        SMSpp_insert_in_factory_h;
    };
    /*---------------------------------------------------------------------------------*/
    /*----------------------LagrangianDualRelaxationSolver-----------------------------*/
    /*---------------------------------------------------------------------------------*/
    class LagrangianDualRelaxationSolver : public RelaxationSolver,
                                           public PrimalProximalHeur
    {
    public:
        // "import" basic types from Block
        using Index = Block::Index;

        /*--------------------------------------------------------------------------*/
        /*------------------------- RESERVED NAMES / KEYS --------------------------*/
        /*--------------------------------------------------------------------------*/
        /// name of the Collection std::map<ColVariable *, gpool_el> that maps the variables to their corresponding solutions
        static constexpr const char *str_VarToSol = "map_varToSol";

        /// key (in map_varToSol) for the Collection of purged columns and their corresponding solutions
        static constexpr const char *str_PurgedColumns = "purgedColumns";

        /// name of the group of dynamic constraint
        static constexpr const char *str_BranchBounds = "BranchBounds";

        /// type T of the Collection of purged columns
        using PurgedColumn = std::map<ColVariable *, std::vector<LagBFunction::gpool_el>>;

        /*-------------------------------------------------------------------------------------*/
        // Adding parameters for branching and applying changes to the master or subproblem
        //  TODO expand with new branching strategies, e.g. most fractional, pseudo-costs, strong branching, etc.
        enum BranchStrategy
        {
            mostFractional = 0, ///< branch on the variable with the most fractional value
        };
        enum ApplyStrategy
        {
            Master = 0,     ///< apply the change to the master problem
            Subproblem = 1, ///< apply the change to the subproblem
        };

        enum int_par_type_LDRS
        {
            ApplyStrategy = intLastPPHPar, ///< strategy for applying the change to the master or subproblem
            BranchStrategy,                ///< branching strategy
            intLastLDRSPar                 ///< first allowed new int parameter for derived classes
        };

        [[nodiscard]] idx_type int_par_first_is(void) const override
        {
            return (intLastLDRSPar);
        }

        void set_par(idx_type par, int val) override;

        [[nodiscard]] int get_dflt_int_par(idx_type par) const override;

        [[nodiscard]] int get_int_par(idx_type par) const override;

        [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override;

        [[nodiscard]] const std::string &int_par_idx2str(idx_type par) const override;

        /*-------------------------------------------------------------------------------------*/

        LagrangianDualRelaxationSolver();

        ~LagrangianDualRelaxationSolver() override;

        int compute(bool changedvars = true) override;

        // TODO understand if it's correct, since we obtain the value from the primal solver
        OFValue get_lb() override;
        OFValue get_ub() override;

        OFValue get_true_lb() override;
        OFValue get_true_ub() override;

        bool has_true_var_solution() override;
        bool new_true_var_solution() override;

        void set_Block(Block *block) override
        {
            PrimalProximalHeur::set_Block(block);
            std::list<FRowConstraint> emptyList = {};
            if (LagrDual && !LagrDual->get_dynamic_constraint<FRowConstraint>(str_BranchBounds))
                LagrDual->add_dynamic_constraint(emptyList, std::string(str_BranchBounds));
        }

        void get_true_var_solution(Configuration *solc = nullptr) override;

        Solution *get_Solution(Configuration *solc = nullptr) override;

        Solution *get_true_solution(Configuration *solc = nullptr) override;

        std::vector<Change *> branch() override;

        Change *apply(Change *change, bool doUndo = false) override;

        // TODO capire se ha senso
        void set_global_information(GlobalInformation *gi) override;

    private:
        // bool PPHdone = false;
        int branchingStrategy = mostFractional;
        int applyStrategy = Master;
        // map the (p)variable to a pair, the first one containig the UB (+1) linear function and the second one containing the LB (-1) linear function
        std::unordered_map<Variable *, std::pair<LinearFunction *, LinearFunction *>> map_varToLF;
        std::shared_ptr<Collection<PurgedColumn>> map_varToSol = nullptr;
        std::map<ColVariable *, AbstractPath> map_varToPath;

        // Remove bound from the LagrangianDualBlock, can return doUndo
        Change *removeBound(ColVariable *pv, bool isLB, bool doUndo,
                            const std::vector<AbstractPath> &paths);

        SMSpp_insert_in_factory_h; // insert LagrangianDualRelaxationSolver in the factory
    };
}