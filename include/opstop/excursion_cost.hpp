#ifndef EXCURSION_COST_H
#define EXCURSION_COST_H

#include <eigen3/Eigen/Core>
#include <ifopt/cost_term.h>
namespace opstop {

class ExcursionCost : public ifopt::CostTerm {
private:
public:
  ExcursionCost();
  virtual ~ExcursionCost() = default;
  double GetCost() const override;
  void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
};
} // namespace opstop

#endif /* EXCURSION_COST_H */
