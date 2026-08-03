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
            eChgLB,          ///< change lower bound of a variable
            eChgUB,          ///< change upper bound of a variable
        };
        /*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

        // constructor
        LagrangianChange(int type, std::vector<double> value, std::vector<AbstractPath> paths) : f_type(type), f_value(std::move(value)), v_paths(std::move(paths)) {}
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
            // read AbstractPath
            auto pg = group.getGroup("VariablesPath");
            if (pg.isNull())
                v_paths.clear();
            else
                AbstractPath::deserialize(v_paths, pg);
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
            if (!v_paths.empty())
            {
                auto pg = group.addGroup("VariablesPath");
                AbstractPath::serialize(v_paths, pg);
            }
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
                // Change Obj require in data[0] index, data[1] new coefficient, data[2] new quadratic coefficient (if quadratic)
                // in this case constains AbstractPath contains only 2 element, the variable and the objective function
            case eChgObj:
            {
                Variable *pv = nullptr;
                Function *fobj = nullptr;
                if (v_paths.size() != 2)
                    throw std::invalid_argument(
                        "LagrangianChange: eChgObj requires 2 AbstractPath elements (variable and objective function)");
                for (const auto &path : v_paths)
                {
                    auto node_type = path.get_last_node(block).type;

                    if (node_type == 'V' || node_type == 'v')
                        pv = path.get_element<Variable>(block);
                    else if (node_type == 'O')
                    {
                        auto *obj = path.get_element<FRealObjective>(block);
                        fobj = obj->get_function();
                    }
                }
                if (!pv || !fobj)
                    throw std::invalid_argument(
                        "LagrangianChange: eChgObj requires a variable and an objective function in AbstractPath");

                // --- Quadratic Case -----------------------------------------------
                if (auto qf = dynamic_cast<DQuadFunction *>(fobj))
                {
                    if (v_data.size() != 2)
                        throw std::invalid_argument(
                            "LagrangianChange: eChgObj on quadratic objective needs 2 values");

                    auto pos = qf->is_active(pv);
                    const auto new_c1 = v_data[0];
                    const auto new_c2 = v_data[1];

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
                                                             {old_c1, old_c2}, std::vector<AbstractPath>{v_paths}));
                }

                // --- Linear Case -----------------------------------------------
                if (auto lf = dynamic_cast<LinearFunction *>(fobj))
                {
                    if (v_data.size() != 1)
                        throw std::invalid_argument(
                            "LagrangianChange: eChgObj on linear objective needs 1 value");

                    const auto new_c1 = v_data[0];
                    auto pos = lf->is_active(pv);

                    const double old_c1 = (pos < lf->get_num_active_var())
                                              ? lf->get_coefficient(pos)
                                              : 0.0;

                    if (pos < lf->get_num_active_var())
                        lf->modify_coefficient(pos, new_c1, issueMod);
                    else
                        lf->add_variable(pv, new_c1, issueMod);

                    if (doUndo)
                        returnChange = (new LagrangianChange(eChgObj, {old_c1}, std::vector<AbstractPath>{v_paths}));
                }

                throw std::invalid_argument(
                    "LagrangianChange: objective Function type not supported");
                break;
            }

            // v_paths contains only the objective function, and v_data[0] contains the new sense
            case eChgSense:
            {
                // apply change to the block
                auto obj = v_paths[0].get_element<Objective>(block);
                if (doUndo)
                    returnChange = new LagrangianChange(eChgSense, std::vector<double>{static_cast<double>(block->get_objective()->get_sense())}, std::vector<AbstractPath>{v_paths});
                obj->set_sense(static_cast<int>(v_data[0]));
                break;
            }
            // v_paths constains the variable, v_data[0] contains the the new integrality
            case eChgIntegrality:
            {
                const bool new_integer = (v_data[0] != 0.0);
                auto pv = v_paths[0].get_element<Variable>(block);
                const bool old_integer = pv->is_integer();

                pv->is_integer(new_integer, issueMod);
                if (doUndo)
                    returnChange = new LagrangianChange(eChgIntegrality,
                                                        {old_integer ? 1.0 : 0.0}, std::vector<AbstractPath>{v_paths});
                break;
            }
            // v_paths constains the variable, v_data[0] contains the new fixed value
            case eFixX:
            {
                const auto fix_value = v_data[0];

                auto pv = v_paths[0].get_element<Variable>(block);
                const bool was_fixed = pv->is_fixed();

                if (doUndo)
                    returnChange = was_fixed ? new LagrangianChange(eFixX, {pv->get_value()}, std::vector<AbstractPath>{v_paths}) : new LagrangianChange(eUnfixX, {}, std::vector<AbstractPath>{v_paths});
                pv->set_value(fix_value);
                pv->is_fixed(true, issueMod);
                break;
            }
            // v_paths constains the variable
            case eUnfixX:
            {
                auto pv = v_paths[0].get_element<Variable>(block);

                pv->is_fixed(false, issueMod);

                if (doUndo)
                    returnChange = new LagrangianChange(eFixX, {pv->get_value()}, std::vector<AbstractPath>{v_paths});
                break;
            }
            // v_paths constains the variable, v_data[0] contains the new lower bound
            case eChgLB:
            {
                auto pv = v_paths[0].get_element<Variable>(block);
                const auto new_lb = v_data[0];
                bool found = false;
                for (Index i = 0; i < pv->get_num_active(); ++i)
                {
                    auto *dep = var->get_active(i);
                    if (auto *c = dynamic_cast<BoxConstraint *>(dep))
                    {
                        if (doUndo)
                            returnChange = new LagrangianChange(eChgLB, {c->get_lhs()}, std::vector<AbstractPath>{v_paths});
                        c->set_lhs(new_lb, issueMod);
                        found = true;
                        break;
                    }
                    else if (auto *c = dynamic_cast<LBConstraint *>(dep))
                    {
                        if (doUndo)
                            returnChange = new LagrangianChange(eChgLB, {c->get_lhs()}, std::vector<AbstractPath>{v_paths});
                        c->set_lhs(new_lb, issueMod);
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    if (doUndo)
                        returnChange = new LagrangianChange(eChgLB, {pv->get_lb()}, std::vector<AbstractPath>{v_paths});
                    auto con = new LBConstraint(pv->get_block(), pv, new_lb);
                    Block *blk = pv->get_block();
                    Index idx = Inf<Index>();

                    auto &d_constraints = blk->get_dynamic_constraints(); // c_Vec_any &

                    for (Index i = 0; i < d_constraints.size(); ++i)
                    {
                        if (d_constraints[i].type() == typeid(std::list<LBConstraint>))
                        {
                            idx = i;
                            break;
                        }
                    }

                    if (idx == Inf<Index>())
                    {
                        // nessun gruppo di LBConstraint dinamici: lo registriamo ora
                        std::list<LBConstraint> empty_list;
                        blk->set_dynamic_constraint(std::move(empty_list), "LBConstraint");
                        idx = d_constraints.size() - 1;
                    }

                    auto *d_list = std::any_cast<std::list<LBConstraint>>(&d_constraints[idx]);
                    // d_list è garantito non-nullptr qui, perché idx è stato appena
                    // verificato/creato per contenere esattamente std::list<LBConstraint>

                    std::list<LBConstraint> newlist;
                    newlist.emplace_back(blk, pv, new_lb);

                    blk->add_dynamic_constraints(*d_list, newlist, issueMod);
                }
                break;
            }
            case eChgUB:
            {
                auto pv = v_paths[0].get_element<Variable>(block);
                const auto new_ub = v_data[0];
                bool found = false;
                for (Index i = 0; i < pv->get_num_active(); ++i)
                {
                    auto *dep = var->get_active(i);
                    if (auto *c = dynamic_cast<BoxConstraint *>(dep))
                    {
                        if (doUndo)
                            returnChange = new LagrangianChange(eChgUB, {c->get_rhs()}, std::vector<AbstractPath>{v_paths});
                        c->set_rhs(new_ub, issueMod);
                        found = true;
                        break;
                    }
                    else if (auto *c = dynamic_cast<UBConstraint *>(dep))
                    {
                        if (doUndo)
                            returnChange = new LagrangianChange(eChgUB, {c->get_rhs()}, std::vector<AbstractPath>{v_paths});
                        c->set_rhs(new_ub, issueMod);
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    if (doUndo)
                        returnChange = new LagrangianChange(eChgUB, {pv->get_ub()}, std::vector<AbstractPath>{v_paths});
                    auto con = new UBConstraint(pv->get_block(), pv, new_ub);
                    auto blk = pv->get_block();
                    Index idx = Inf<Index>();
                    for (Index i = 0; i < d_constraints.size(); ++i)
                    {
                        if (d_constraints[i].type() == typeid(std::list<UBConstraint>))
                        {
                            idx = i;
                            break;
                        }
                    }
                    if (idx == Inf<Index>())
                    {
                        // nessun gruppo di UBConstraint dinamici: lo registriamo ora
                        std::list<UBConstraint> empty_list;
                        blk->set_dynamic_constraint(std::move(empty_list), "UBConstraint");
                        idx = d_constraints.size() - 1;
                    }
                    auto *d_list = std::any_cast<std::list<UBConstraint>>(&d_constraints[idx]);
                    std::list<UBConstraint> newlist;
                    newlist.emplace_back(blk, pv, new_ub);
                    blk->add_dynamic_constraints(*d_list, std::move(newlist), issueMod);
                }
                break;
            }
            default:
                throw std::invalid_argument("LagrangianChange: unknown change type");
            }
            return returnChange;

        protected:
            int f_type;                        ///< type of the change
            std::vector<double> v_data;        ///< value of the change
            std::vector<AbstractPath> v_paths; ///< vector of abstract path (variables and constraints) involved in the change

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
