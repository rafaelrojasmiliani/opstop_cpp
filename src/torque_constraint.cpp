
#include <opstop/torque_constraint.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>

TorqueConstraint::TorqueConstraint(
    const gsplines::functions::FunctionExpression &_curve, std::size_t _nglp,
    double _ti, std::vector<double> &_bound)
    : ConstraintSet(_curve.get_codom_dim() * _nglp, "torque_constraint"),
      helper_(_curve, _nglp, _ti),
      torque_buff_(_nglp * _curve.get_codom_dim()) {

  assert(_bound.size() == _curve.get_codom_dim());

  for (std::size_t i = 0; i < _curve.get_codom_dim(); i++) {
    ifopt::Bounds b(-_bound[i], _bound[i]);
    bounds_vector_.push_back(b);
  }
}

Eigen::VectorXd TorqueConstraint::GetValues() const {
  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::rnea(model_, data_, helper_.q_val_buff_.row(uici).transpose(),
                    helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
                    helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());
    torque_buff_.segment(uici * helper_.position_->get_codom_dim(),
                         helper_.nglp_) = data_.tau;
  }

  return torque_buff_;
}

ifopt::Component::VecBound TorqueConstraint::GetBounds() const {
  return bounds_vector_;
}

void TorqueConstraint::FillJacobianBlock(std::string _set_name,
                                         Jacobian &_jac_block) const {
  for (std::size_t uici = 0; uici < helper_.nglp_; uici++) {
    pinocchio::computeRNEADerivatives(
        model_, data_, helper_.q_val_buff_.row(uici).transpose(),
        helper_.q_diff_1_wrt_t_buff_.row(uici).transpose(),
        helper_.q_diff_2_wrt_t_buff_.row(uici).transpose());
    Eigen::MatrixXd result =
        data_.dtau_dq * helper_.q_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.dtau_dv *
            helper_.q_diff_1_wrt_t_diff_wrt_sf_buff_.row(uici).transpose() +
        data_.M *
            helper_.q_diff_2_wrt_t_diff_wrt_sf_buff_.row(uici).transpose();
  }
}
