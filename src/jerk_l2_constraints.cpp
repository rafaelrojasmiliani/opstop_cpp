
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <opstop/jerk_l2_constraints.hpp>
namespace opstop {

JerkL2Constraints::JerkL2Constraints(
    const gsplines::functions::FunctionBase &_curve, std::size_t _nglp,
    double _ti, double _bound)
    : ConstraintSet(1, "jerk_l2_constraints"), helper_(_curve, _nglp, _ti),
      value_buff_(1), bounds_vector_({ifopt::Bounds(-ifopt::inf, _bound)}),
      alpha_(0.0), glw_(_nglp) {

  std::tie(std::ignore, glw_) =
      gsplines::collocation::legendre_gauss_lobatto_points_and_weights(_nglp);
}

Eigen::VectorXd JerkL2Constraints::GetValues() const {
  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  double result =
      (helper_.q_diff_3_wrt_t_buff_.array().pow(2.0).rowwise().sum().array() *
       glw_.array())
          .array()
          .sum();

  value_buff_(0) = result;
  return value_buff_;
}

ifopt::Component::VecBound JerkL2Constraints::GetBounds() const {
  return bounds_vector_;
}

void JerkL2Constraints::FillJacobianBlock(std::string _set_name,
                                          Jacobian &_jac_block) const {

  Eigen::VectorXd x =
      GetVariables()->GetComponent("parametrization_variables")->GetValues();
  helper_.set_diffeo(x(0), x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_diff_3_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_3_wrt_t_partial_diff_wrt_sf();

  Eigen::Map<const Eigen::VectorXd> col_0(
      helper_.q_diff_3_wrt_t_diff_wrt_Ts_buff_.data(), GetRows());
  Eigen::Map<const Eigen::VectorXd> col_1(
      helper_.q_diff_3_wrt_t_diff_wrt_sf_buff_.data(), GetRows());

  for (std::size_t uici = 0; uici < GetRows(); uici++) {
    _jac_block.coeffRef(uici, 0) = col_0(uici);
    _jac_block.coeffRef(uici, 1) = col_1(uici);
  }
}

Eigen::VectorXd JerkL2Constraints::__GetValues(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  double result =
      (helper_.q_diff_3_wrt_t_buff_.array().pow(2.0).rowwise().sum().array() *
       glw_.array())
          .array()
          .sum();

  value_buff_(0) = result * 0.5;
  return value_buff_;
}

Eigen::MatrixXd JerkL2Constraints::__GetJacobian(Eigen::Vector2d &_x) const {

  helper_.set_diffeo(_x(0), _x(1));
  helper_.compute_s_and_its_derivatives_wrt_tau();
  helper_.compute_q_and_its_derivatives_wrt_s();
  helper_.compute_q_and_its_derivatives_wrt_t();
  helper_.compute_s_and_its_derivatives_wrt_Ts_sf();
  helper_.compute_q_diff_3_wrt_t_partial_diff_wrt_Ts();
  helper_.compute_q_diff_3_wrt_t_partial_diff_wrt_sf();

  Eigen::MatrixXd result(GetRows(), 2);

  result.col(0) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_3_wrt_t_diff_wrt_Ts_buff_.data(), GetRows());
  result.col(1) = Eigen::Map<const Eigen::VectorXd>(
      helper_.q_diff_3_wrt_t_diff_wrt_sf_buff_.data(), GetRows());

  result(0, 0) = ((helper_.q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() *
                   helper_.q_diff_3_wrt_t_buff_.array())
                      .rowwise()
                      .sum()
                      .array() *
                  glw_.array())
                     .array()
                     .sum();

  result(0, 1) = ((helper_.q_diff_3_wrt_t_diff_wrt_sf_buff_.array() *
                   helper_.q_diff_3_wrt_t_buff_.array())
                      .rowwise()
                      .sum()
                      .array() *
                  glw_.array())
                     .array()
                     .sum();

  return result;
}

} // namespace opstop
