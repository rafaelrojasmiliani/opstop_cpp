#include <opstop/velocity_constraint.hpp>
namespace opstop {

VelocityConstraints::VelocityConstraints(
    const gsplines::functions::FunctionExpression &_curve, std::size_t _nglp,
    double _ti, std::vector<double> &_bound)
    : ConstraintSet(_curve.get_codom_dim() * _nglp, "velocity_constraints"),
      helper_(_curve, _nglp, _ti), value_buff_(_nglp * _curve.get_codom_dim()) {

  for (std::size_t i = 0; i < _curve.get_codom_dim(); i++) {
    ifopt::Bounds b(-_bound[i], _bound[i]);
    for (std::size_t j = 0; j < _nglp; j++) {
      bounds_vector_.push_back(b);
    }
  }
}

Eigen::VectorXd VelocityConstraints::GetValues() const {
  // helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_Ts();
  return Eigen::Map<const Eigen::VectorXd>(helper_.q_diff_1_wrt_t_buff_.data(),
                                           helper_.q_diff_1_wrt_t_buff_.size());
}

ifopt::Component::VecBound VelocityConstraints::GetBounds() const {
  return bounds_vector_;
}

void VelocityConstraints::FillJacobianBlock(std::string _set_name,
                                            Jacobian &_jac_block) const {}

Eigen::VectorXd VelocityConstraints::__GetValues(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  return Eigen::Map<const Eigen::VectorXd>(helper_.q_diff_1_wrt_t_buff_.data(),
                                           helper_.q_diff_1_wrt_t_buff_.size());
}

Eigen::MatrixXd VelocityConstraints::__GetJacobian(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_sf();

  Eigen::MatrixXd result(GetRows(), 2);

  result.col(0) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_1_wrt_t_diff_wrt_Ts_buff_.data(), GetRows());
  result.col(1) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_1_wrt_t_diff_wrt_sf_buff_.data(), GetRows());

  return result;
}
} // namespace opstop
