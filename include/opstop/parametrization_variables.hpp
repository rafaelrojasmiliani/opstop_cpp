#ifndef PARAMETRIZATION_VARIABLES_H
#define PARAMETRIZATION_VARIABLES_H

#include <gsplines/Functions/FunctionBase.hpp>
#include <ifopt/variable_set.h>
namespace opstop {

class ParametrizationVariables : public ifopt::VariableSet {
private:
  Eigen::Vector2d values_;
  const double ti_;
  const double exec_time_;
  ifopt::Component::Component::VecBound bounds_;

public:
  ParametrizationVariables(double _ti, double _exec_time);
  ParametrizationVariables(const ParametrizationVariables &_that);
  virtual ~ParametrizationVariables() {}
  void SetVariables(const Eigen::VectorXd &_vec) override;
  Eigen::VectorXd GetValues() const override;
  ifopt::Component::VecBound GetBounds() const override;
  bool
  set_for_acceleration_bounds(const gsplines::functions::FunctionBase &_trj,
                              double _ti, const Eigen::VectorXd &_acc_bounds);
};

} // namespace opstop
#endif /* PARAMETRIZATION_VARIABLES_H */
