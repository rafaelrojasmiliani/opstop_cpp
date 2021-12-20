#include <iostream>
#include <opstop/parametrization_variables.hpp>
namespace opstop {
ParametrizationVariables::ParametrizationVariables(double _ti,
                                                   double _exec_time)
    : VariableSet(2, "parametrization_variables"), ti_(_ti),
      exec_time_(_exec_time), bounds_({ifopt::Bounds(0.0, ifopt::inf),
                                       ifopt::Bounds(ti_, exec_time_)}) {

  double xi = 2;
  double eta = 3.0 / 5.0 * xi;
  values_(0) = (_exec_time - _ti) * xi;
  values_(1) = (_exec_time - _ti) * eta + _ti;
}

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

bool ParametrizationVariables::set_for_acceleration_bounds(
    const gsplines::functions::FunctionBase &_trj, double _ti,
    const Eigen::VectorXd &_acc_bounds) {

  std::unique_ptr<gsplines::functions::FunctionBase> vel = _trj.deriv(1);
  std::unique_ptr<gsplines::functions::FunctionBase> acc = _trj.deriv(2);

  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
      500, _trj.get_domain().first, _trj.get_domain().second);

  Eigen::MatrixXd vel_val = vel->value(time_spam);
  Eigen::MatrixXd acc_val = acc->value(time_spam);

  Eigen::VectorXd xi_vec(Eigen::VectorXd::Zero(_trj.get_codom_dim()));
  for (std::size_t i = 0; i < _trj.get_codom_dim(); i++) {

    double vel_max = vel_val.col(i).array().abs().maxCoeff();
    double acc_max = acc_val.col(i).array().abs().maxCoeff();

    if (acc_max > _acc_bounds(i))
      return false;

    xi_vec(i) = 16.0 / 9.0 * vel_max / (_acc_bounds(i) - acc_max) /
                (_trj.get_domain().second - _ti);

    if (xi_vec(i) > 2.0) {
      return false;
    }
  }

  double xi = 2.0;
  double eta = 3.0 / 5.0 * xi;
  values_(0) = (_trj.get_domain().second - _ti) * xi;
  values_(1) = (_trj.get_domain().second - _ti) * eta + _ti;
  return true;
}
} // namespace opstop
