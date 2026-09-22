#include "LagrangianDualRelaxationSolver.h"

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using OFValue = Solver::OFValue;

/*---------------------------------------------------------------------------------*/
/*-------------------------------- LagrangianChange --------------------------------*/
/*---------------------------------------------------------------------------------*/

void LagrangianChange::deserialize(const netCDF::NcGroup &group)
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
        v_paths = AbstractPath::vector_deserialize(pg);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void LagrangianChange::serialize(netCDF::NcGroup &group) const
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

Change *LagrangianChange::apply(Block *block, bool doUndo,
                                ModParam issueMod,
                                ModParam issueAMod)
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
        ColVariable *pv = nullptr;
        Function *fobj = nullptr;
        if (v_paths.size() != 2)
            throw std::invalid_argument(
                "LagrangianChange: eChgObj requires 2 AbstractPath elements (variable and objective function)");
        for (const auto &path : v_paths)
        {
            auto node_type = path.get_last_node(block).type;

            if (node_type == 'V' || node_type == 'v')
                pv = path.get_element<ColVariable>(block);
            else if (node_type == 'O')
            {
                auto *obj = dynamic_cast<FRealObjective *>(
                    path.get_element<Objective>(block));
                if (obj)
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
        else if (auto lf = dynamic_cast<LinearFunction *>(fobj))
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
        else
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
        auto pv = v_paths[0].get_element<ColVariable>(block);
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

        auto pv = v_paths[0].get_element<ColVariable>(block);
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
        auto pv = v_paths[0].get_element<ColVariable>(block);

        pv->is_fixed(false, issueMod);

        if (doUndo)
            returnChange = new LagrangianChange(eFixX, {pv->get_value()}, std::vector<AbstractPath>{v_paths});
        break;
    }
    // v_paths constains the variable, v_data[0] contains the new lower bound
    case eChgLB:
    {
        auto pv = v_paths[0].get_element<ColVariable>(block);
        const auto new_lb = v_data[0];
        bool found = false;
        for (Index i = 0; i < pv->get_num_active(); ++i)
        {
            auto *dep = pv->get_active(i);
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
            throw std::invalid_argument("LagrangianChange: eChgLB requires a variable with an existing BoxConstraint or LBConstraint");
            /*                     if (doUndo)
                                    returnChange = new LagrangianChange(eChgLB, {pv->get_lb()}, std::vector<AbstractPath>{v_paths});
                                // auto con = new LBConstraint(pv->get_Block(), pv, new_lb);
                                Block *blk = pv->get_Block();
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
                                        throw std::invalid_argument("No dynamic LBConstraint group found in block. Please register a dynamic LBConstraint group before applying LagrangianChange eChgLB.");
                                }

                                auto *d_list = boost::any_cast<std::list<LBConstraint>>(&d_constraints[idx]);
                                // d_list è garantito non-nullptr qui, perché idx è stato appena
                                // verificato/creato per contenere esattamente std::list<LBConstraint>

                                std::list<LBConstraint> newlist;
                                newlist.emplace_back(blk, pv, new_lb);

                                blk->add_dynamic_constraints(*d_list, newlist, issueMod);
             */
        }
        break;
    }
    case eChgUB:
    {
        auto pv = v_paths[0].get_element<ColVariable>(block);
        const auto new_ub = v_data[0];
        bool found = false;
        for (Index i = 0; i < pv->get_num_active(); ++i)
        {
            auto *dep = pv->get_active(i);
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
            throw std::invalid_argument("LagrangianChange: eChgUB requires a variable with an existing BoxConstraint or UBConstraint");
            /*                     if (doUndo)
                                    returnChange = new LagrangianChange(eChgUB, {pv->get_ub()}, std::vector<AbstractPath>{v_paths});
                                // auto con = new UBConstraint(pv->get_Block(), pv, new_ub);
                                auto blk = pv->get_Block();
                                Index idx = Inf<Index>();
                                auto &d_constraints = blk->get_dynamic_constraints();
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
                                    throw std::invalid_argument("No dynamic UBConstraint group found in block. Please register a dynamic UBConstraint group before applying LagrangianChange eChgUB.");
                                }
                                auto *d_list = boost::any_cast<std::list<UBConstraint>>(&d_constraints[idx]);
                                std::list<UBConstraint> newlist;
                                newlist.emplace_back(blk, pv, new_ub);
                                blk->add_dynamic_constraints(*d_list, std::move(newlist), issueMod);
             */
        }
        break;
    }
    default:
    {
        throw std::invalid_argument("LagrangianChange: unknown change type");
    }
    }
    return returnChange;
}

/*---------------------------------------------------------------------------------*/
/*----------------------LagrangianDualRelaxationSolver-----------------------------*/
/*---------------------------------------------------------------------------------*/

void LagrangianDualRelaxationSolver::set_par(idx_type par, int val)
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

int LagrangianDualRelaxationSolver::get_dflt_int_par(idx_type par) const
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

int LagrangianDualRelaxationSolver::get_int_par(idx_type par) const
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

LagrangianDualRelaxationSolver::idx_type LagrangianDualRelaxationSolver::int_par_str2idx(const std::string &name) const
{
    if (name == "ApplyStrategy")
        return ApplyStrategy;
    if (name == "BranchStrategy")
        return BranchStrategy;
    return PrimalProximalHeur::int_par_str2idx(name);
}

const std::string &LagrangianDualRelaxationSolver::int_par_idx2str(idx_type par) const
{
    static const std::string apply_strategy_name = "ApplyStrategy";
    static const std::string branch_strategy_name = "BranchStrategy";
    switch (par)
    {
    case ApplyStrategy:
        return apply_strategy_name;
    case BranchStrategy:
        return branch_strategy_name;
    default:
        return PrimalProximalHeur::int_par_idx2str(par);
    }
}

/*-------------------------------------------------------------------------------------*/

// TODO fathom node through cutoff (dblupcutoff/dbldowncutoff)

LagrangianDualRelaxationSolver::LagrangianDualRelaxationSolver()
    : PrimalProximalHeur(),
      // RelaxationSolver(),
      // PPHdone(false),
      branchingStrategy(mostFractional),
      applyStrategy(Master),
      map_varToLF()
{
}

LagrangianDualRelaxationSolver::~LagrangianDualRelaxationSolver()
{
}

int LagrangianDualRelaxationSolver::compute(bool changedvars)
{
    // PPHdone = false;
    int status = this->PrimalProximalHeur::compute(changedvars);
    if (status >= ThinComputeInterface::kOK && (status == Solver::kLowPrecision || status < ThinComputeInterface::kError) && get_Lagrangian_initial_solution().size() == NumStatVar)
        return ThinComputeInterface::kOK;
    return status;
}

// TODO understand if it's correct, since we obtain the value from the primal solver
OFValue LagrangianDualRelaxationSolver::get_lb() { return this->PrimalProximalHeur::LagrangianDualSolver::get_lb(); }
OFValue LagrangianDualRelaxationSolver::get_ub() { return this->PrimalProximalHeur::LagrangianDualSolver::get_ub(); }

OFValue LagrangianDualRelaxationSolver::get_true_lb()
{
    /*             if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                } */
    return this->PrimalProximalHeur::get_lb();
}
OFValue LagrangianDualRelaxationSolver::get_true_ub()
{
    /*             if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                } */
    return this->PrimalProximalHeur::get_ub();
}

bool LagrangianDualRelaxationSolver::has_true_var_solution()
{
    /*             if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                } */
    return this->PrimalProximalHeur::has_var_solution();
}
bool LagrangianDualRelaxationSolver::new_true_var_solution()
{
    /*             if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                } */
    return this->PrimalProximalHeur::new_var_solution();
}

void LagrangianDualRelaxationSolver::get_true_var_solution(Configuration *solc)
{
    /*             if (!PPHdone)
                {
                    this->PrimalProximalHeur::compute();
                    PPHdone = true;
                } */
    this->PrimalProximalHeur::get_var_solution(solc);
}

Solution *LagrangianDualRelaxationSolver::get_Solution(Configuration *solc)
{
    f_Block->lock(this);
    this->LagrangianDualSolver::get_var_solution(solc);
    auto solution = f_Block->get_Solution(solc);
    f_Block->unlock(this);
    return solution;
}

Solution *LagrangianDualRelaxationSolver::get_true_solution(Configuration *solc)
{
    f_Block->lock(this);
    get_true_var_solution(solc);
    auto solution = f_Block->get_Solution(solc);
    f_Block->unlock(this);
    return solution;
}

std::vector<Change *> LagrangianDualRelaxationSolver::branch()
{
    std::vector<Change *> changes;
    // obtain the solution of the lagrangian relaxation
    const auto &init_sol = get_Lagrangian_initial_solution();
    if (init_sol.empty())
        throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: no initial solution available"));
    if (init_sol.size() != NumStatVar)
        throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: initial solution size does not match number of variables"));

    // obtain the variable for choose the branching variable
    ColVariable *mostFracVar = nullptr;
    double bestCriterionValue = -1.0;
    double branchValue = 0.0;

    Index kVar = 0;
    switch (branchingStrategy)
    {
    case mostFractional:
    {
        // Definition of sol (which is used for init_sol)
        //  for( Index index = 0 ; index < n_sub ; ++index )
        //      for( Index ivar = 0 ; ivar < pos_id_sbi[ index ] ; ++ivar )
        //      {
        //          sol[ kvar ] = idx_to_var_sbi1[ index ][ ivar ].second->get_value();
        //          kvar++;
        //      }

        for (const auto &sbd : idx_to_var_sbi1) // for each subblock
        {
            for (const auto &dv : sbd) // for each variable in the subblock
            {
                const double value = init_sol[kVar++];
                double fractionality = std::abs(value - std::round(value));
                if (fractionality > bestCriterionValue)
                {
                    bestCriterionValue = fractionality;
                    mostFracVar = dv.second;
                    branchValue = value;
                }
            }
        }
        break;
    }
    default:
        throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: branching strategy not implemented"));
    }
    if (!mostFracVar)
        throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: no fractional variable found"));

    switch (applyStrategy)
    {
    case Master:
        changes.push_back(new LagrangianChange(LagrangianChange::eChgUB, {std::floor(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        changes.push_back(new LagrangianChange(LagrangianChange::eChgLB, {std::ceil(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        return changes;
    case Subproblem:
        changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {std::floor(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {std::ceil(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        // changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {1}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        // changes.push_back(new LagrangianChange(LagrangianChange::eFixX, {0}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));

        return changes;
    default:
        throw(std::runtime_error("LagrangianDualRelaxationSolver::branch: apply strategy not implemented"));
    }
}

Change *LagrangianDualRelaxationSolver::apply(Change *change, bool doUndo)
{
    Change *undoChange = nullptr;
    if (!change)
        throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: null change"));
    auto c = dynamic_cast<LagrangianChange *>(change);
    if (!c)
        throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: change is not a LagrangianChange"));
    auto pv = c->get_paths()[0].get_element<ColVariable>(this->f_Block);
    if (!pv)
        throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: variable not found in block"));
    auto lbf = v_LBF[Block2Index(pv->get_Block())];
    if (c->get_type() == LagrangianChange::eChgLB)
    {
        auto value = c->get_data()[0];

        switch (applyStrategy)
        {
        case Master:
        {
            // Index found_pos = Inf<Index>();
            double oldLB = pv->get_lb();
            if (map_varToLF.contains(pv) && map_varToLF[pv]->get_coefficient(0) == -1.0)
            {
                // trovata: aggiorna solo il termine costante (cioè "value")
                auto *gi = map_varToLF[pv];
                oldLB = gi->get_constant_term();
                gi->set_constant_term(value);
            }
            else
            {
                // g(x) = value -x = -1*x + value
                auto *g = new LinearFunction(LinearFunction::v_coeff_pair{{pv, -1.0}}, //-1 coefficent
                                             value);                                   // constant term
                ColVariable *y = new ColVariable();
                y->is_positive(true);
                lbf->add_dual_pairs(LagBFunction::v_dual_pair{{y, g}});
                map_varToLF[pv] = g;
            }
            if (doUndo)
            {
                undoChange = new LagrangianChange(LagrangianChange::eChgLB, {oldLB}, std::vector<AbstractPath>{c->get_paths()});
            }
            break;
        }
        case Subproblem:
        {
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented yet for eChgLB"));
        }
        default:
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
        }
    }
    else if (c->get_type() == LagrangianChange::eChgUB)
    {
        auto value = c->get_data()[0];

        switch (applyStrategy)
        {
        case Master:
        {
            double oldUB = pv->get_ub();
            if (map_varToLF.contains(pv) && map_varToLF[pv]->get_coefficient(0) == 1.0)
            {
                // trovata: aggiorna solo il termine costante (cioè "value")
                auto *gi = map_varToLF[pv];
                oldUB = gi->get_constant_term() * -1.0; // store the old upper bound
                gi->set_constant_term(-value);
            }
            else
            {
                // g(x) = x - value = 1*x - value

                auto *g = new LinearFunction(LinearFunction::v_coeff_pair{{pv, 1.0}}, // 1 coefficent
                                             -value);                                 // constant term
                ColVariable *y = new ColVariable();
                y->is_positive(true);
                lbf->add_dual_pairs(LagBFunction::v_dual_pair{{y, g}});
                map_varToLF[pv] = g;
            }
            if (doUndo)
            {
                undoChange = new LagrangianChange(LagrangianChange::eChgUB, {oldUB}, std::vector<AbstractPath>{c->get_paths()});
            }
            break;
        }
        case Subproblem:
        {
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented yet for eChgUB"));
        }
        default:
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
        }
    }
    else if (c->get_type() == LagrangianChange::eFixX)
    {
        auto value = c->get_data()[0];
        switch (applyStrategy)
        {
        case Master:
        {
            // possible to apply in the block, let's see
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented yet for eFixX"));
        }
        case Subproblem:
        {

            if (doUndo)
            {
                if (pv->is_fixed())
                {
                    undoChange = new LagrangianChange(LagrangianChange::eFixX, {pv->get_value()}, std::vector<AbstractPath>{c->get_paths()});
                }
                else
                {
                    undoChange = new LagrangianChange(LagrangianChange::eUnfixX, {}, std::vector<AbstractPath>{c->get_paths()});
                }
            }
            auto mvts = map_varToSol;
            const auto handler_id = lbf->set_event_handler(LagBFunction::eColumnPurged,
                                                           // TODO insert here lambda function
                                                           [lbf, pv, mvts]() -> int
                                                           {
                                                               auto el = lbf->release_current_purged_solution();
                                                               if (!el.sol)
                                                                   throw std::logic_error("eColumnPurged called without a Solution");
                                                               if (!mvts)
                                                                   throw std::runtime_error("LagrangianDualRelaxationSolver::apply: map_varToSol not set in global information");
                                                               // TODO controllare che sia giusta e funzioni
                                                               auto ok = mvts->write_with(str_PurgedColumns,
                                                                                          [pv, &el](PurgedColumn &pc)
                                                                                          {
                                                                                              pc[pv].push_back(std::move(el));
                                                                                          });
                                                               // temporarly, check if it works correctly
                                                               assert(ok);
                                                               return ThinComputeInterface::eContinue;
                                                           });
            pv->set_value(value);
            pv->is_fixed(true);
            lbf->reset_event_handler(LagBFunction::eColumnPurged, handler_id);
            break;
        }
        default:
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
        }
    }
    else if (c->get_type() == LagrangianChange::eUnfixX)
    {
        switch (applyStrategy)
        {
        case Master:
        {
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented yet for eUnfixX"));
        }
        case Subproblem:
        {
            if (doUndo)
            {
                undoChange = new LagrangianChange(LagrangianChange::eFixX, {pv->get_value()}, std::vector<AbstractPath>{c->get_paths()});
            }
            pv->is_fixed(false);
            if (!f_global_information)
                throw std::runtime_error("LagrangianDualRelaxationSolver::apply: global information not set");

            if (!map_varToSol)
            {
                // TODO capire se è sbagliato in quanto non ho rimosso colonne
                throw std::runtime_error("LagrangianDualRelaxationSolver::apply: map_varToSol not found in global information");
            }
            std::vector<LagBFunction::gpool_el> pruned_solutions;
            const bool found = map_varToSol->read_with(str_PurgedColumns,
                                                       [pv, &pruned_solutions](const PurgedColumn &pc)
                                                       {
                                                           auto it = pc.find(pv);
                                                           if (it != pc.end())
                                                               pruned_solutions = it->second;
                                                           else
                                                               pruned_solutions.clear();
                                                       });

            if (!found)
                pruned_solutions.clear();
            lbf->restore_purged_solutions(std::move(pruned_solutions));

            break;
        }
        default:
            throw(std::runtime_error("LagrangianDualRelaxationSolver::apply: apply strategy not implemented"));
        }
    }
    else
        undoChange = c->apply(this->f_Block, doUndo);

    return undoChange;
}

// TODO capire se ha senso
void LagrangianDualRelaxationSolver::set_global_information(GlobalInformation *gi)
{
    ChangeSolver::set_global_information(gi);
    map_varToSol = nullptr;
    if (!gi)
        return;
    if (!gi->exists(str_VarToSol))
    {
        gi->add_to_Universe<PurgedColumn>(str_VarToSol);
    }
    map_varToSol = gi->get_from_Universe<PurgedColumn>(str_VarToSol);
    if (map_varToSol && !map_varToSol->contains(str_PurgedColumns))
    {
        map_varToSol->write(str_PurgedColumns, PurgedColumn{});
    }
}
// register LagrangianDualRelaxationSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0(LagrangianDualRelaxationSolver);
SMSpp_insert_in_factory_cpp_0(LagrangianChange);
