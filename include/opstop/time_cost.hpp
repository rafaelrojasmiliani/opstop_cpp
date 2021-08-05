#ifndef TIME_COST_H
#define TIME_COST_H

#include <eigen3/Eigen/Core>
#include <ifopt/cost_term.h>
class TimeCost : public ifopt::CostTerm {
private:
public:
  TimeCost();
  virtual ~TimeCost() {}
  double GetCost() const override;
  void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
};

#endif /* TIME_COST_H */
