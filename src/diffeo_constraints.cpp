#include <opstop/diffeo_constraints.hpp>

namespace opstop {

DiffeoConstraints::DiffeoConstraints(double _ti, double _exec_time)
    : ConstraintSet(2, "diffeo_constraints"), ti_(_ti), exec_time_(_exec_time),
      value_(2) {

  bounds_vector_.push_back(ifopt::Bounds(0.0, ifopt::inf));
  bounds_vector_.push_back(ifopt::Bounds(0.0, ifopt::inf));
}
Eigen::VectorXd DiffeoConstraints::GetValues() const {
  // helper_.set_diffeo(_x(0), _x(1));
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();

  double Ts = x(0);
  double sf = x(1);

  value_(0) = 3.0 * Ts - 5.0 * sf + 5.0 * ti_;
  value_(1) = -2.0 * Ts + 5.0 * sf - 5.0 * ti_;

  return value_;
}

ifopt::Component::VecBound DiffeoConstraints::GetBounds() const {
  return bounds_vector_;
}

void DiffeoConstraints::FillJacobianBlock(std::string _set_name,
                                          Jacobian &_jac_block) const {

  _jac_block.coeffRef(0, 0) = 3.0;
  _jac_block.coeffRef(0, 1) = -5.0;

  _jac_block.coeffRef(1, 0) = -2.0;
  _jac_block.coeffRef(1, 1) = 5.0;
}
} // namespace opstop
