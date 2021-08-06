#include <opstop/diffeo_constraints.hpp>

DiffeoConstraints::DiffeoConstraints(double _ti, double _exec_time)
    : ConstraintSet(4, "diffeo_constraints"), ti_(_ti), exec_time_(_exec_time),
      value_(4) {}
Eigen::VectorXd DiffeoConstraints::GetValues() const {
  // helper_.set_diffeo(_x(0), _x(1));
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();

  double Ts = x(0);
  double sf = x(1);
  value_(0) = sf - ti_;
  value_(1) = sf - exec_time_;
  value_(2) = Ts;
  value_(3) = 8.0 * Ts - 15.0 * sf + 15.0 * ti_;

  return value_;
}

ifopt::Component::VecBound DiffeoConstraints::GetBounds() const {
  return bounds_vector_;
}

void DiffeoConstraints::FillJacobianBlock(std::string _set_name,
                                          Jacobian &_jac_block) const {

  _jac_block.coeffRef(0, 1) = 1.0;
  _jac_block.coeffRef(1, 1) = 1.0;
  _jac_block.coeffRef(2, 0) = 1.0;
  _jac_block.coeffRef(3, 0) = 8.0;
  _jac_block.coeffRef(3, 1) = -15.0;
}
