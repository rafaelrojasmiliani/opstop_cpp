
#include <opstop/torque_constraint.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
namespace opstop {

TorqueConstraint::TorqueConstraint(
    const gsplines::functions::FunctionExpression &_curve, std::size_t _nglp,
    double _ti, const std::vector<double> &_bound,
    const pinocchio::Model &_model)
    : ConstraintSet(_curve.get_codom_dim() * _nglp, "torque_constraint"),
      helper_(_curve, _nglp, _ti), torque_buff_(_nglp * _curve.get_codom_dim()),
      model_(_model), data_(model_) {

  assert(_bound.size() == _curve.get_codom_dim());

  for (std::size_t i = 0; i < _curve.get_codom_dim(); i++) {
    ifopt::Bounds b(-_bound[i], _bound[i]);
    for (std::size_t j = 0; j < _nglp; j++) {
      bounds_vector_.push_back(b);
    }
  }
}

Eigen::VectorXd TorqueConstraint::GetValues() const {
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::rnea(model_, data_, helper_.q_val_buff_.row(uici).transpose(),
                    helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
                    helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());
  }

  return torque_buff_;
}

ifopt::Component::VecBound TorqueConstraint::GetBounds() const {
  return bounds_vector_;
}

void TorqueConstraint::FillJacobianBlock(std::string _set_name,
                                         Jacobian &_jac_block) const {
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  std::size_t codom_dim = helper_.position_->get_codom_dim();
  Eigen::MatrixXd result(helper_.nglp_ * codom_dim, 2);
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();

  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_partial_diff_wrt_Ts();
  helper_.compute_q_partial_diff_wrt_sf();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_sf();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_sf();

  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::computeRNEADerivatives(
        model_, data_, helper_.q_val_buff_.row(uici).transpose(),
        helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
        helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());
    result.col(0).segment(uici * codom_dim, codom_dim) =
        data_.dtau_dq * helper_.q_diff_wrt_Ts_buff_.row(uici).transpose() +
        data_.dtau_dv *
            helper_.q_diff_1_wrt_t_diff_wrt_Ts_buff_.row(uici).transpose() +
        data_.M.selfadjointView<Eigen::Upper>() *
            helper_.q_diff_2_wrt_t_diff_wrt_Ts_buff_.row(uici).transpose();

    result.col(1).segment(uici * codom_dim, codom_dim) =
        data_.dtau_dq * helper_.q_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.dtau_dv *
            helper_.q_diff_1_wrt_t_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.M.selfadjointView<Eigen::Upper>() *
            helper_.q_diff_2_wrt_t_diff_wrt_sf_buff_.row(uici).transpose();
  }
  //_jac_block = result;
}

Eigen::VectorXd TorqueConstraint::__GetValues(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();

  std::size_t codom_dim = helper_.position_->get_codom_dim();
  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::rnea(model_, data_, helper_.q_val_buff_.row(uici).transpose(),
                    helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
                    helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());

    torque_buff_.segment(uici * codom_dim, codom_dim) = data_.tau;
  }

  return torque_buff_;
}
Eigen::MatrixXd TorqueConstraint::__GetJacobian(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();

  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_partial_diff_wrt_Ts();
  helper_.compute_q_partial_diff_wrt_sf();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_1_wrt_t_partial_diff_wrt_sf();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_2_wrt_t_partial_diff_wrt_sf();

  std::size_t codom_dim = helper_.position_->get_codom_dim();
  Eigen::MatrixXd result(helper_.nglp_ * codom_dim, 2);

  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::computeRNEADerivatives(
        model_, data_, helper_.q_val_buff_.row(uici).transpose(),
        helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
        helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());

    result.col(0).segment(uici * codom_dim, codom_dim) =
        data_.dtau_dq * helper_.q_diff_wrt_Ts_buff_.row(uici).transpose() +
        data_.dtau_dv *
            helper_.q_diff_1_wrt_t_diff_wrt_Ts_buff_.row(uici).transpose() +
        data_.M.selfadjointView<Eigen::Upper>() *
            helper_.q_diff_2_wrt_t_diff_wrt_Ts_buff_.row(uici).transpose();

    result.col(1).segment(uici * codom_dim, codom_dim) =
        data_.dtau_dq * helper_.q_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.dtau_dv *
            helper_.q_diff_1_wrt_t_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.M.selfadjointView<Eigen::Upper>() *
            helper_.q_diff_2_wrt_t_diff_wrt_sf_buff_.row(uici).transpose();
  }
  return result;
}
} // namespace opstop
