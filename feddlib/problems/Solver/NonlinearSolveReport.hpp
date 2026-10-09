#ifndef FEDD_NONLINEAR_SOLVE_REPORT_HPP
#define FEDD_NONLINEAR_SOLVE_REPORT_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_decl.hpp>
#include <cmath>

namespace FEDD {

/** @brief Outcome of one nonlinear solve, independently of cancellation policy. */
struct NonlinearSolveReport {
    bool converged = false;
    bool linearSolvesConverged = true;
    bool finite = true;
    int nonlinearIterations = 0;
    int linearSolveCount = 0;
    int totalLinearIterations = 0;
    int maximumLinearIterations = 0; ///< Zero for direct solvers or unspecified limits.
    double finalCriterion = 0.;

    double averageLinearIterations() const {
        return linearSolveCount ? double(totalLinearIterations) / linearSolveCount : 0.;
    }
    bool reliable() const { return converged && linearSolvesConverged && finite; }
};

/** @brief Read the iteration limit of the selected Belos method, not another solver's limit. */
inline int configuredLinearIterationLimit(const Teuchos::ParameterList& parameters) {
    if (!parameters.isSublist("ThyraSolver")) return 0;
    const auto& thyra = parameters.sublist("ThyraSolver");
    if ((thyra.isParameter("Linear Solver Type") && thyra.get<std::string>("Linear Solver Type") != "Belos") ||
        !thyra.isSublist("Linear Solver Types")) return 0;
    const auto& types = thyra.sublist("Linear Solver Types");
    if (!types.isSublist("Belos")) return 0;
    const auto& belos = types.sublist("Belos");
    if (!belos.isSublist("Solver Types")) return 0;
    const auto& solvers = belos.sublist("Solver Types");
    std::string method = "Block GMRES";
    if (belos.isParameter("Solver Type")) {
        if (belos.isType<std::string>("Solver Type")) method = belos.get<std::string>("Solver Type");
        else method = Teuchos::getStringValue<Thyra::EBelosSolverType>(belos, "Solver Type");
    }
    if (!solvers.isSublist(method)) return 0;
    const auto& selected = solvers.sublist(method);
    return selected.isParameter("Maximum Iterations") ? selected.get<int>("Maximum Iterations") : 1000;
}
} // namespace FEDD
#endif
