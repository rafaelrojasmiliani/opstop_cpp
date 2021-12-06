#include <opstop/excursion_cost.hpp>
namespace opstop {

ExcursionCost::ExcursionCost() : CostTerm("excursion_cost") {}
double ExcursionCost::GetCost() const {
  return GetVariables()
      ->GetComponent("parametrization_variables")
      ->GetValues()(1);
}

void ExcursionCost::FillJacobianBlock(std::string _var_set,
                                      Jacobian &_jac) const {

  _jac.coeffRef(0, 1) = 1;
  _jac.coeffRef(0, 0) = 0;
}
} // namespace opstop
