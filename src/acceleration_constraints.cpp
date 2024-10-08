#include <exception>
#include <opstop/acceleration_constraints.hpp>

namespace opstop {

AccelerationConstraints::AccelerationConstraints(
    const gsplines::functions::FunctionBase &_curve, std::size_t _nglp,
    double _ti, const std::vector<double> &_bound)
    : ConstraintSet(_curve.get_codom_dim() * _nglp, "acceleration_constraints"),
      helper_(_curve, _nglp, _ti), value_buff_(_nglp * _curve.get_codom_dim()) {

  if (_curve.get_codom_dim() != _bound.size()) {
    throw std::logic_error("Incorrect number of constraints");
  }

  for (std::size_t i = 0; i < _curve.get_codom_dim(); i++) {
    ifopt::Bounds b(-_bound[i], _bound[i]);
    for (std::size_t j = 0; j < _nglp; j++) {
      bounds_vector_.push_back(b);
    }
  }
}

Eigen::VectorXd AccelerationConstraints::GetValues() const {
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  return Eigen::Map<const Eigen::VectorXd>(helper_.q_diff_2_wrt_t_buff_.data(),
                                           helper_.q_diff_2_wrt_t_buff_.size());
}

ifopt::Component::VecBound AccelerationConstraints::GetBounds() const {
  return bounds_vector_;
}

void AccelerationConstraints::FillJacobianBlock(std::string _set_name,
                                                Jacobian &_jac_block) const {

  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_sf();

  Eigen::Map<const Eigen::VectorXd> col_0(
      helper_.q_diff_2_wrt_t_diff_wrt_Ts_buff_.data(), GetRows());
  Eigen::Map<const Eigen::VectorXd> col_1(
      helper_.q_diff_2_wrt_t_diff_wrt_sf_buff_.data(), GetRows());

  for (std::size_t uici = 0; uici < GetRows(); uici++) {
    _jac_block.coeffRef(uici, 0) = col_0(uici);
    _jac_block.coeffRef(uici, 1) = col_1(uici);
  }
}

Eigen::VectorXd
AccelerationConstraints::__GetValues(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  return Eigen::Map<const Eigen::VectorXd>(helper_.q_diff_2_wrt_t_buff_.data(),
                                           helper_.q_diff_2_wrt_t_buff_.size());
}

Eigen::MatrixXd
AccelerationConstraints::__GetJacobian(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_sf();

  Eigen::MatrixXd result(GetRows(), 2);

  result.col(0) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_2_wrt_t_diff_wrt_Ts_buff_.data(), GetRows());
  result.col(1) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_2_wrt_t_diff_wrt_sf_buff_.data(), GetRows());

  return result;
}
} // namespace opstop
