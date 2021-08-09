#include <opstop/diffeo_constraints.hpp>

DiffeoConstraints::DiffeoConstraints(double _ti, double _exec_time)
    : ConstraintSet(1, "diffeo_constraints"), ti_(_ti), exec_time_(_exec_time),
      value_(1) {

  bounds_vector_.push_back(ifopt::Bounds(ifopt::inf, 0.0));
}
Eigen::VectorXd DiffeoConstraints::GetValues() const {
  // helper_.set_diffeo(_x(0), _x(1));
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();

  double Ts = x(0);
  double sf = x(1);
  value_(0) = 8.0 * Ts - 15.0 * sf + 15.0 * ti_;

  return value_;
}

ifopt::Component::VecBound DiffeoConstraints::GetBounds() const {
  return bounds_vector_;
}

void DiffeoConstraints::FillJacobianBlock(std::string _set_name,
                                          Jacobian &_jac_block) const {

  _jac_block.coeffRef(0, 0) = 8.0;
  _jac_block.coeffRef(0, 1) = -15.0;
}
