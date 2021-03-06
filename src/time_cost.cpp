#include <opstop/time_cost.hpp>
namespace opstop {

TimeCost::TimeCost() : CostTerm("time_cost") {}
double TimeCost::GetCost() const {
  return GetVariables()
      ->GetComponent("parametrization_variables")
      ->GetValues()(0);
}

void TimeCost::FillJacobianBlock(std::string _var_set, Jacobian &_jac) const {

  _jac.coeffRef(0, 0) = 1;
  _jac.coeffRef(0, 1) = 0;
}
} // namespace opstop
