#include "LagrangianDualRelaxationSolver.h"

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using OFValue = Solver::OFValue;

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
      map_varToLF(),
      map_varToPath()
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
        changes.push_back(new LagrangianChange(AbstractChange::eChgUB, {std::floor(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, LagrDual)}));
        changes.push_back(new LagrangianChange(AbstractChange::eChgLB, {std::ceil(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, LagrDual)}));
        return changes;
    case Subproblem:
        changes.push_back(new LagrangianChange(AbstractChange::eFixX, {std::floor(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, LagrDual)}));
        changes.push_back(new LagrangianChange(AbstractChange::eFixX, {std::ceil(branchValue)}, std::vector<AbstractPath>{AbstractPath(mostFracVar, LagrDual)}));
        // changes.push_back(new LagrangianChange(AbstractChange::eFixX, {1}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));
        // changes.push_back(new LagrangianChange(AbstractChange::eFixX, {0}, std::vector<AbstractPath>{AbstractPath(mostFracVar, this->f_Block)}));

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
    auto pv = c->get_paths()[0].get_element<ColVariable>(LagrDual);
    if (!pv)
        throw(std::invalid_argument("LagrangianDualRelaxationSolver::apply: variable not found in block"));
    auto lbf = v_LBF[Block2Index(pv->get_Block())];
    if (c->get_type() == AbstractChange::eChgLB)
    {
        auto value = c->get_data()[0];

        switch (applyStrategy)
        {
        case Master:
        {
            // Index found_pos = Inf<Index>();
            double oldLB = pv->get_lb();
            if (map_varToLF.contains(pv) && map_varToLF[pv].second != nullptr /*&& map_varToLF[pv].second->get_coefficient(0) == -1.0*/)
            {
                // trovata: aggiorna solo il termine costante (cioè "value")
                auto *gi = map_varToLF[pv].second;
                oldLB = gi->get_constant_term();
                gi->set_constant_term(value);
                if (doUndo)
                {
                    undoChange = new LagrangianChange(AbstractChange::eChgLB, {oldLB}, std::vector<AbstractPath>{c->get_paths()});
                }
            }
            else
            {
                // g(x) = value -x = -1*x + value
                auto *g = new LinearFunction(LinearFunction::v_coeff_pair{{pv, -1.0}}, //-1 coefficent
                                             value);                                   // constant term

                std::list<FRowConstraint> cons_list;
                cons_list.emplace_back(LagrDual, -Inf<RowConstraint::RHSValue>(), 0, g);

                auto *group = LagrDual->get_dynamic_constraint<FRowConstraint>(str_BranchBounds);
                if (!group)
                    throw std::runtime_error("group BranchBounds not found in block");

                LagrDual->add_dynamic_constraints(*group, cons_list); // plurale: aggiunge al gruppo esistente

                if (group->empty())
                    throw std::runtime_error("Invalid BranchBounds group");

                map_varToLF[pv].second = static_cast<LinearFunction *>(group->back().get_function());

                if (doUndo)
                {
                    undoChange = new LagrangianChange(LagrangianChange::eDeleteLB, {}, std::vector<AbstractPath>{c->get_paths()});
                }
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
    else if (c->get_type() == AbstractChange::eChgUB)
    {
        auto value = c->get_data()[0];

        switch (applyStrategy)
        {
        case Master:
        {
            double oldUB = pv->get_ub();
            if (map_varToLF.contains(pv) && map_varToLF[pv].first != nullptr /* && map_varToLF[pv].first->get_coefficient(0) == 1.0 */)
            {
                // trovata: aggiorna solo il termine costante (cioè "value")
                auto *gi = map_varToLF[pv].first;
                oldUB = gi->get_constant_term() * -1.0; // store the old upper bound
                gi->set_constant_term(-value);
                if (doUndo)
                {
                    undoChange = new LagrangianChange(AbstractChange::eChgUB, {oldUB}, std::vector<AbstractPath>{c->get_paths()});
                }
            }
            else
            {
                // g(x) = x - value = 1*x - value

                auto *g = new LinearFunction(LinearFunction::v_coeff_pair{{pv, 1.0}}, // 1 coefficent
                                             -value);
                // constant term
                std::list<FRowConstraint> cons_list;
                cons_list.emplace_back(LagrDual, -Inf<RowConstraint::RHSValue>(), 0, g);

                auto *group = LagrDual->get_dynamic_constraint<FRowConstraint>(str_BranchBounds);
                if (!group)
                    throw std::runtime_error("group BranchBounds not found in block");

                LagrDual->add_dynamic_constraints(*group, cons_list); // plurale: aggiunge al gruppo esistente

                if (group->empty())
                    throw std::runtime_error("Invalid BranchBounds group");

                map_varToLF[pv].first = static_cast<LinearFunction *>(group->back().get_function());
                if (doUndo)
                {
                    undoChange = new LagrangianChange(LagrangianChange::eDeleteUB, {}, std::vector<AbstractPath>{c->get_paths()});
                }
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
    else if (c->get_type() == AbstractChange::eFixX)
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
                    undoChange = new LagrangianChange(AbstractChange::eFixX, {pv->get_value()}, std::vector<AbstractPath>{c->get_paths()});
                }
                else
                {
                    undoChange = new LagrangianChange(AbstractChange::eUnfixX, {}, std::vector<AbstractPath>{c->get_paths()});
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
    else if (c->get_type() == AbstractChange::eUnfixX)
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
                undoChange = new LagrangianChange(AbstractChange::eFixX, {pv->get_value()}, std::vector<AbstractPath>{c->get_paths()});
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
    else if (c->get_type() == LagrangianChange::eDeleteUB)
    {
    }
    else if (c->get_type() == LagrangianChange::eDeleteLB)
    {
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

// nel .cpp
Change *LagrangianDualRelaxationSolver::removeBound(
    ColVariable *pv, bool isLB, bool doUndo,
    const std::vector<AbstractPath> &paths)
{
    auto &slot = isLB ? map_varToLF[pv].second : map_varToLF[pv].first;

    if (!slot)
        throw std::logic_error(
            "LagrangianDualRelaxationSolver::removeBound: no activate constraint");

    auto *con = dynamic_cast<FRowConstraint *>(slot->get_Observer());
    if (!con)
        throw std::logic_error(
            "LagrangianDualRelaxationSolver::removeBound: Observer is not a FRowConstraint");

    Change *undoChange = nullptr;
    if (doUndo)
    {
        double oldValue = isLB ? slot->get_constant_term()
                               : -slot->get_constant_term();
        undoChange = new LagrangianChange(
            isLB ? AbstractChange::eChgLB : AbstractChange::eChgUB,
            {oldValue}, paths);
    }

    auto *group = LagrDual->get_dynamic_constraint<FRowConstraint>(str_BranchBounds);
    auto it = std::find_if(group->begin(), group->end(),
                           [con](const FRowConstraint &c2)
                           { return &c2 == con; });
    if (it == group->end())
        throw std::logic_error(
            "LagrangianDualRelaxationSolver::removeBound: constraint not found in the group");

    LagrDual->remove_dynamic_constraint(*group, it);

    slot = nullptr;
    return undoChange;
}

// register LagrangianDualRelaxationSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0(LagrangianDualRelaxationSolver);
SMSpp_insert_in_factory_cpp_0(LagrangianChange);
