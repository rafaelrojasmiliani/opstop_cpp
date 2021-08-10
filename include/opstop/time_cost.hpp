#ifndef TIME_COST_H
#define TIME_COST_H

#include <eigen3/Eigen/Core>
#include <ifopt/cost_term.h>
namespace opstop {

class TimeCost : public ifopt::CostTerm {
private:
public:
  TimeCost();
  virtual ~TimeCost() = default;
  double GetCost() const override;
  void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
};
} // namespace opstop

#endif /* TIME_COST_H */
