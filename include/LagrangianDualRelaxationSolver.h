#pragma once

#include "PrimalProximalHeur.h"
#include "RelaxationSolver.h"

namespace SMSpp_di_unipi_it
{

    // by claude, da controllare
    static ColVariable *get_static_variable_by_index(Block *block, Index idx)
    {
        ColVariable *result = nullptr;
        Index cur = 0;

        for (const auto &el : block->get_static_variables())
        {
            if (result)
                break;

            if (un_any_thing_0(ColVariable, el, {
                    if (cur == idx)
                        result = &var;
                    ++cur;
                }))
                continue;

            if (un_any_thing_1(ColVariable, el, {
                    for (Index j = 0; j < var.size(); ++j, ++cur)
                        if (cur == idx)
                        {
                            result = var.data() + j;
                            break;
                        }
                }))
                continue;

            if (un_any_thing_K(ColVariable, el, {
                    for (Index j = 0; j < var.num_elements(); ++j, ++cur)
                        if (cur == idx)
                        {
                            result = var.data() + j;
                            break;
                        }
                }))
                continue;
        }
        if (!result)
            throw std::invalid_argument(
                "LagrangianChange: variable index out of range");

        return (result);
    }

    class LagrangianChange : public Change
    {
        enum LagrangianChangeType
        {
            eEmpty = 0       ///< empty change, used for initialization
            eChgObj,         ///< change objective coefficient of a variable
            eChgSense,       ///< change sense of the objective
            eChgIntegrality, ///< change integrality of a variable
            eFixX,           ///< fix a variable to a value
            eUnfixX,         ///< unfix a variable
            // eChgLB,          ///< change lower bound of a variable
            // eChgUB,          ///< change upper bound of a variable
        };
        /*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

        // constructor
        LagrangianChange(int type, std::vector<double> value) : f_type(type), f_value(value) {}
        // decostructor
        ~LagrangianChange() = default;
        /*-------------------- PUBLIC METHODS OF THE CLASS -------------------------*/
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        void deserialize(const netCDF::NcGroup &group) override
        {
            auto ftype = group.getAtt("LagrangianChange_type");
            if (ftype.isNull())
                throw std::invalid_argument("LagrangianChange_type attribute not found in netCDF group");
            ftype.getValues(&f_type);
            // read data
            netCDF::NcDim ni = group.getDim("dim");
            netCDF::NcVar data = group.getVar("Data");
            if (data.isNull())
                v_data.clear();
            else
            {
                v_data.resize(ni.getSize());
                data.getVar(v_data.data());
            }
        }
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        void serialize(netCDF::NcGroup &group) const override
        {

            // always call the method of the base class first
            Change::serialize(group);

            group.putAtt("LagrangianChange_type", netCDF::NcInt(), f_type);

            netCDF::NcDim ni = group.addDim("dim", v_data.size());
            (group.addVar("Data",
                          netCDF::NcDouble(), ni))
                .putVar(v_data.data());
        }
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        Change *apply(Block *block, bool doUndo = false,
                      ModParam issueMod = eNoBlck,
                      ModParam issueAMod = eNoBlck) override
        {
            Change *returnChange = nullptr;
            switch (f_type)
            {
            case eEmpty:
                throw std::invalid_argument("LagrangianChange: empty change cannot be applied");
            case eChgObj:
            {
                const auto block = Index(v_data[0]);
                const auto idx = Index(v_data[1]);
                auto obj = static_cast<FRealObjective *>(block->get_objective());
                auto fobj = obj->get_function();
                auto pv = get_static_variable_by_index(block, idx);

                // --- Quadratic Case -----------------------------------------------
                if (auto qf = dynamic_cast<DQuadFunction *>(fobj))
                {
                    if (v_data.size() != 4)
                        throw std::invalid_argument(
                            "LagrangianChange: eChgObj on quadratic objective needs 4 values");

                    const auto new_c1 = v_data[2];
                    const auto new_c2 = v_data[3];

                    auto pos = qf->is_active(pv);

                    double old_c1 = 0.0, old_c2 = 0.0;
                    if (pos < qf->get_num_active_var())
                    {
                        old_c1 = qf->get_linear_coefficient(pos);
                        old_c2 = qf->get_quadratic_coefficient(pos);
                        qf->modify_term(pos, new_c1, new_c2, issueMod);
                    }
                    else
                        qf->add_variable(pv, new_c1, new_c2, issueMod);

                    if (doUndo)
                        returnChange = (new LagrangianChange(eChgObj,
                                                             {double(block), double(idx), old_c1, old_c2}));
                }

                // --- Linear Case -----------------------------------------------
                if (auto lf = dynamic_cast<LinearFunction *>(fobj))
                {
                    if (v_data.size() != 3)
                        throw std::invalid_argument(
                            "LagrangianChange: eChgObj on linear objective needs 3 values");

                    const auto new_c1 = v_data[2];
                    auto pos = lf->is_active(pv);

                    const double old_c1 = (pos < lf->get_num_active_var())
                                              ? lf->get_coefficient(pos)
                                              : 0.0;

                    if (pos < lf->get_num_active_var())
                        lf->modify_coefficient(pos, new_c1, issueMod);
                    else
                        lf->add_variable(pv, new_c1, issueMod);

                    if (doUndo)
                        returnChange = (new LagrangianChange(eChgObj, {double(block), double(idx), old_c1}));
                }

                throw std::invalid_argument(
                    "LagrangianChange: objective Function type not supported");
            }

            case eChgSense:
            {
                // apply change to the block
                if (doUndo)
                    returnChange = new LagrangianChange(eChgSense, std::vector<double>{v_data[0], static_cast<double>(block->get_objective()->get_sense())});
                block->get_objective()->set_sense(static_cast<int>(v_data[1]));
            }
            case eChgIntegrality:
            {
                const auto block = Index(v_data[0]);
                const auto idx = Index(v_data[1]);
                const bool new_integer = (v_data[2] != 0.0);

                auto pv = get_static_variable_by_index(block, idx);
                const bool old_integer = pv->is_integer();

                pv->is_integer(new_integer, issueMod);
                if (doUndo)
                    returnChange = new LagrangianChange(eChgIntegrality,
                                                        {double(block), double(idx), old_integer ? 1.0 : 0.0});
            }
            case eFixX:
            {
                const auto block = Index(v_data[0]);
                const auto idx = Index(v_data[1]);
                const auto fix_value = v_data[2];

                auto pv = get_static_variable_by_index(block, idx);
                const bool was_fixed = pv->is_fixed();
                const double old_value = pv->get_value();

                pv->set_value(fix_value);
                pv->is_fixed(true, issueMod);

                if (doUndo)
                    returnChange = was_fixed ? new LagrangianChange(eFixX, {double(block), double(idx), old_value}) : new LagrangianChange(eUnfixX, {double(block), double(idx)});
            }
            case eUnfixX:
            {
                const auto block = Index(v_data[0]);
                const auto idx = Index(v_data[1]);

                auto pv = get_static_variable_by_index(block, idx);
                const double old_value = pv->get_value();

                pv->is_fixed(false, issueMod);

                if (doUndo)
                {
                    returnChange = new LagrangianChange(eFixX, {double(block), double(idx), old_value});
                }
            }
            default:
                throw std::invalid_argument("LagrangianChange: unknown change type");
            }
            return returnChange;

        protected:
            int f_type;                 ///< type of the change
            std::vector<double> v_data; ///< value of the change
        private:
        }

        SMSpp_insert_in_factory_h;
    }; // end of class LagrangianChange

    /*---------------------------------------------------------------------------------*/
    /*----------------------LagrangianDualRelaxationSolver-----------------------------*/
    /*---------------------------------------------------------------------------------*/
    class LagrangianDualRelaxationSolver : public RelaxationSolver, PrimalProximalHeur
    {
    public:
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

        void set_par override(idx_type par, int val)
        {
            switch (par)
            {
            case ApplyStrategy:
                applyStrategy = val;
                break;
            case BranchStrategy:
                branchingStrategy = val;
                break;
            default:
                PrimalProximalHeur::set_par(par, val);
                break;
            }
        }
        [[nodiscard]] int get_dflt_int_par(idx_type par) const override
        {
            switch (par)
            {
            case ApplyStrategy:
                return Master;
            case BranchStrategy:
                return mostFractional;
            default:
                return PrimalProximalHeur::get_dflt_int_par(par);
            }
        }

        [[nodiscard]] int get_int_par(idx_type par) const override
        {
            switch (par)
            {
            case ApplyStrategy:
                return applyStrategy;
            case BranchStrategy:
                return branchingStrategy;
            default:
                return PrimalProximalHeur::get_int_par(par);
            }
        }

        [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override
        {
            if (name == "ApplyStrategy")
                return ApplyStrategy;
            if (name == "BranchStrategy")
                return BranchStrategy;
            return PrimalProximalHeur::int_par_str2idx(name);
        }

        [[nodiscard]] std::string int_par_idx2str(idx_type par) const override
        {
            switch (par)
            {
            case ApplyStrategy:
                return "ApplyStrategy";
            case BranchStrategy:
                return "BranchStrategy";
            default:
                return PrimalProximalHeur::int_par_idx2str(par);
            }

            /*-------------------------------------------------------------------------------------*/

            LagrangianDualRelaxationSolver() : PrimalProximalHeur(),
                                               // RelaxationSolver(),
                                               PPHdone(false),
                                               branchingStrategy(mostFractional),
                                               applyStrategy(Master)
            {
            }

            ~LagrangianDualRelaxationSolver() override
            {
            }

            int compute(bool changedvars = true) override
            {
                PPHdone = false;
                int status = this->PrimalProximalHeur::LagrangianDualSolver::compute(changedvars);

                return status;
            }

            OFValue get_lb() override { return this->PrimalProximalHeur::LagrangianDualSolver::get_lb(); }
            OFValue get_ub() override { return this->PrimalProximalHeur::LagrangianDualSolver::get_ub(); }

            OFValue get_true_lb() override
            {
                if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                }
                return this->PrimalProximalHeur::get_lb();
            }
            OFValue get_true_ub() override
            {
                if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                }
                return this->PrimalProximalHeur::get_ub();
            }

            bool has_true_var_solution() override
            {
                if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                }
                return this->PrimalProximalHeur::has_var_solution();
            }
            bool new_true_var_solution() override
            {
                if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                }
                return this->PrimalProximalHeur::new_var_solution();
            }

            void get_true_var_solution(Configuration *solc = nullptr) override
            {
                if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                }
                this->PrimalProximalHeur::get_var_solution(solc);
            }

            std::vector<Change *> branch() override
            {
                std::vector<Change *> changes;
                switch (branchingStrategy)
                {
                case mostFractional:
                {
                    ColVariable *mostFracVar = nullptr;
                    double maxFractionality = -1.0;
                    Index sb = 0, selectedSB = 0;
                    for (const auto &sbd : idx_to_var_sbi1) // for each subblock
                    {
                        for (const auto &dv : sbd) // for each variable in the subblock
                        {
                            const auto pv = dv.second;
                            double value = pv->get_value();
                            double fractionality = std::abs(value - std::round(value));
                            if (fractionality > maxFractionality)
                            {
                                maxFractionality = fractionality;
                                mostFracVar = pv;
                                selectedSB = sb;
                            }
                        }
                        ++sb;
                    }
                    break;
                }
                default:
                    throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: branching strategy not implemented"));
                }

                changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {double(selectedSB), double(mostFracVar->get_index()), std::ceil(mostFracVar->get_value())}));
                changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {double(selectedSB), double(mostFracVar->get_index()), std::floor(mostFracVar->get_value())}));
                return changes;
            }

            Change *apply(Change * change, bool doUndo = false) override
            {
                Change *undoChange = nullptr;
                if (!change)
                    throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: null change"));
                auto c = dynamic_cast<LagrangianChange *>(change);
                if (!c)
                    throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: change is not a LagrangianChange"));
                if (c->get_type() == LagrangianChange::eFixX)
                {
                    switch (applyStrategy)
                    {
                    case Master:
                    {
                        auto index = c->get_data()[1];
                        auto value = c->get_data()[2];
                        auto lbf = v_LBF[Index(c->get_data()[0])];
                        auto inner = lbf->get_inner_block();
                        auto pv = get_static_variable_by_index(inner, Index(index));
                        // build the contraints
                        auto con = new BoxConstraint(pv, value, value);
                        lbf->add_dynamic_constraint(std::list<Constraint *>{con});
                        if (doUndo)
                        {
                            undoChange = new LagrangianChange(LagrangianChange::eUnfixX, {c->get_data()[0], c->get_data()[1]});
                        }
                    }
                    case Subproblem:
                    {
                        auto index = c->get_data()[1];
                        auto value = c->get_data()[2];
                        auto lbf = v_LBF[Index(c->get_data()[0])];
                        auto inner = lbf->get_inner_block();
                        auto pv = get_static_variable_by_index(inner, Index(index));

                        pv->fix_value(value);
                        if (doUndo)
                        {
                            undoChange = new LagrangianChange(LagrangianChange::eUnfixX, {c->get_data()[0], c->get_data()[1]});
                        }
                    }
                    default:
                        throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
                    }
                    else if (c->get_type() == LagrangianChange::eUnfixX)
                    {
                        switch (applyStrategy)
                        {
                        case Master:
                        {
                            auto value = c->get_data()[2];
                            auto lbf = v_LBF[Index(c->get_data()[0])];
                            auto pv = get_static_variable_by_index(lbf->get_inner_block(), Index(c->get_data()[1]));
                            // remove the constraint
                            bool found = false;
                            BoxConstraint *removeCon = nullptr;
                            for (Index i = 0; i < lbf->get_dynamic_constraints().size(); ++i)
                            {
                                auto con = lbf->get_dynamic_constraints()[i];
                                if (auto box = dynamic_cast<BoxConstraint *>(con))
                                {
                                    if (box->is_active(pv) && box->get_rhs() == value && box->get_lhs() == value)
                                    {
                                        found = true;
                                        removeCon = box;
                                        if (doUndo)
                                        {
                                            undoChange = new LagrangianChange(LagrangianChange::eFixX, {c->get_data()[0], c->get_data()[1], value});
                                        }
                                        break;
                                    }
                                }
                            }
                            if (!found)
                                throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: constraint to unfix not found"));
                            lbf->remove_dynamic_constraint(new std::list<Constraint *>{removeCon});
                            delete removeCon;
                        }
                        case Subproblem:
                        {
                            auto index = c->get_data()[1];
                            auto lbf = v_LBF[Index(c->get_data()[0])];
                            auto inner = lbf->get_inner_block();
                            auto pv = get_static_variable_by_index(inner, Index(index));
                            if (doUndo)
                            {
                                undoChange = new LagrangianChange(LagrangianChange::eFixX, {c->get_data()[0], c->get_data()[1], pv->get_value()});
                            };
                        }
                            pv->unfix_value();
                        }
                    default:
                        throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
                    }
                    else undoChange = c->apply(this->f_Block, doUndo);
                }
                return undoChange;
            }

        private:
            bool PPHdone = false;
            int branchingStrategy = mostFractional;
            int applyStrategy = Master;
        };
    }
}
