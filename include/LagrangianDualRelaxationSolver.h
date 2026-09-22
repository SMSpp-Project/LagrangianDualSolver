#pragma once

#include "PrimalProximalHeur.h"
#include "AbstractPath.h"
#include "ChangeSolver.h"

namespace SMSpp_di_unipi_it
{
    class LagrangianChange : public Change
    {
    public:
        enum LagrangianChangeType
        {
            eEmpty = 0,      ///< empty change, used for initialization
            eChgObj,         ///< change objective coefficient of a variable
            eChgSense,       ///< change sense of the objective
            eChgIntegrality, ///< change integrality of a variable
            eFixX,           ///< fix a variable to a value
            eUnfixX,         ///< unfix a variable
            eChgLB,          ///< change lower bound of a variable
            eChgUB           ///< change upper bound of a variable
        };
        /*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

        // constructor
        LagrangianChange() : f_type(eEmpty), v_data(), v_paths() {}

        LagrangianChange(int type, std::vector<double> value, std::vector<AbstractPath> paths)
            : f_type(type), v_data(std::move(value)), v_paths(std::move(paths)) {}

        // decostructor
        ~LagrangianChange() = default;

        /*-------------------- PUBLIC METHODS OF THE CLASS -------------------------*/
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        void deserialize(const netCDF::NcGroup &group) override;

        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        void serialize(netCDF::NcGroup &group) const override;

        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        Change *apply(Block *block, bool doUndo = false,
                      ModParam issueMod = eNoBlck,
                      ModParam issueAMod = eNoBlck) override;

        // getter
        [[nodiscard]] int get_type() const { return f_type; }
        [[nodiscard]] const std::vector<double> &get_data() const { return v_data; }
        [[nodiscard]] const std::vector<AbstractPath> &get_paths() const { return v_paths; }

    protected:
        int f_type;                        ///< type of the change
        std::vector<double> v_data;        ///< value of the change
        std::vector<AbstractPath> v_paths; ///< vector of abstract path (variables and constraints) involved in the change

    private:
        SMSpp_insert_in_factory_h;
    }; // end of class LagrangianChange

    /*---------------------------------------------------------------------------------*/
    /*----------------------LagrangianDualRelaxationSolver-----------------------------*/
    /*---------------------------------------------------------------------------------*/
    class LagrangianDualRelaxationSolver : public RelaxationSolver, public PrimalProximalHeur
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
        std::unordered_map<Variable *, LinearFunction *> map_varToLF;
        std::shared_ptr<Collection<PurgedColumn>> map_varToSol = nullptr;

        SMSpp_insert_in_factory_h; // insert LagrangianDualRelaxationSolver in the factory
    };
}