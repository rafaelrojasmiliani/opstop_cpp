#include <opstop/parametrization_variables.hpp>

ParametrizationVariables::ParametrizationVariables(double _ti,
                                                   double _exec_time)
    : VariableSet(2, "ParametrizationVariables"), ti_(_ti),
      exec_time_(_exec_time), bounds_({ifopt::Bounds(0.0, ifopt::inf),
                                       ifopt::Bounds(ti_, exec_time_)}) {}

void ParametrizationVariables::SetVariables(const Eigen::VectorXd &_vec) {
  values_ = _vec;
}

Eigen::VectorXd ParametrizationVariables::GetValues() const { return values_; }
ifopt::Component::VecBound ParametrizationVariables::GetBounds() const {
  return bounds_;
}
ParametrizationVariables::ParametrizationVariables(
    const ParametrizationVariables &_that)
    : VariableSet(_that), ti_(_that.ti_), exec_time_(_that.exec_time_),
      bounds_(_that.bounds_), values_(_that.values_) {}
