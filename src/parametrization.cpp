#include <opstop/parametrization.hpp>

ParametrizedCurveHelper::ParametrizedCurveHelper(
    std::unique_ptr<const gsplines::functions::FunctionExpression> _curve,
    std::size_t _nglp, double _ti)
    : position_(_curve->clone()), velocity_(_curve->deriv()),
      acceleration_(velocity_->deriv()), jerk_(acceleration_->deriv()),
      snap_(jerk_->deriv()) {

  for (std::size_t coeff = 0; coeff < 5; coeff++) {

    for (std::size_t deriv = 0; deriv < 4; deriv++) {

      s_bar_deriv_coeff[deriv][coeff] = 0.0;

      s_bar_deriv_coeff_diff_wrt_si[deriv][coeff] = 0.0;

      s_bar_deriv_coeff_diff_wrt_Ts[deriv][coeff] = 0.0;
    }
  }

  std::tie(glp_, std::ignore) =
      gsplines::collocation::legendre_gauss_lobatto_points_and_weights(_nglp);

  s_bar_deriv_coeff_diff_wrt_si[0][0] = 0.0;
  s_bar_deriv_coeff_diff_wrt_si[0][1] = 0.0;
  s_bar_deriv_coeff_diff_wrt_si[0][3] = 10.0;
  s_bar_deriv_coeff_diff_wrt_si[0][4] = -15.0;
  s_bar_deriv_coeff_diff_wrt_si[0][5] = 6.0;

  s_bar_deriv_coeff_diff_wrt_Ts[0][0] = 0.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][1] = 1.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][3] = -6.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][4] = 8.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][5] = -3.0;

  for (std::size_t deriv = 1; deriv < 4; deriv++) {
    for (std::size_t coeff = 0; coeff < 6 - deriv; coeff++) {
      s_bar_deriv_coeff_diff_wrt_si[deriv][coeff] =
          s_bar_deriv_coeff_diff_wrt_si[deriv - 1][coeff + 1] * (coeff + 1);

      s_bar_deriv_coeff_diff_wrt_Ts[deriv][coeff] =
          s_bar_deriv_coeff_diff_wrt_Ts[deriv - 1][coeff + 1] * (coeff + 1);
    }
  }
}
void ParametrizedCurveHelper::set_diffeo(double _Ts, double _sf) {

  s_bar_deriv_coeff[0][0] = ti_;
  s_bar_deriv_coeff[0][1] = _Ts;
  s_bar_deriv_coeff[0][3] = -6.0 * _Ts + 10.0 * _sf - 10.0 * ti_;
  s_bar_deriv_coeff[0][4] = 8.0 * _Ts - 15.0 * _sf + 15 * ti_;
  s_bar_deriv_coeff[0][5] = -3.0 * _Ts + 6.0 * _sf - 6.0 * ti_;

  for (std::size_t deriv = 1; deriv < 4; deriv++) {
    for (std::size_t coeff = 0; coeff < 6 - deriv; coeff++) {
      s_bar_deriv_coeff[deriv][coeff] =
          s_bar_deriv_coeff[deriv - 1][coeff + 1] * (coeff + 1);
    }
  }

  sf_ = _sf;
  Ts_ = _Ts;
}

void ParametrizedCurveHelper::compute_acceleration() {}
void ParametrizedCurveHelper::compute_jerk() {}
void ParametrizedCurveHelper::compute_snap() {}

void ParametrizedCurveHelper::compute_q_and_its_derivatives_wrt_s() {

  position_->value(s_val_buff_, q_val_buff_);
  velocity_->value(s_val_buff_, q_diff_1_wrt_s_buff_);
  acceleration_->value(s_val_buff_, q_diff_2_wrt_s_buff_);
  jerk_->value(s_val_buff_, q_diff_3_wrt_s_buff_);
  snap_->value(s_val_buff_, q_diff_4_wrt_s_buff_);
}
void ParametrizedCurveHelper::compute_s_and_its_derivatives_wrt_tau() {

  eval_pol_deg_5(glp_, s_bar_deriv_coeff[0], s_val_buff_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff[1], s_diff_1_wrt_tau_buff_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff[2], s_diff_2_wrt_tau_buff_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff[3], s_diff_3_wrt_tau_buff_);
}

void ParametrizedCurveHelper::compute_q_and_its_derivatives_wrt_t() {

  q_diff_1_wrt_t_buff_.array() = q_diff_1_wrt_s_buff_.array().rowwise() *
                                 s_diff_1_wrt_tau_buff_.transpose().array() /
                                 Ts_;
  // ---- Second Derivative -----
  q_diff_2_wrt_s_buff_.array() =
      q_diff_2_wrt_s_buff_.array().rowwise() *
          Eigen::pow(s_diff_1_wrt_tau_buff_.transpose().array(), 2) +
      q_diff_1_wrt_s_buff_.array().rowwise() *
          s_diff_2_wrt_tau_buff_.transpose().array();

  q_diff_2_wrt_t_buff_ /= std::pow(Ts_, 2);

  // ---- Third Derivative -----
  q_diff_3_wrt_t_buff_ =
      q_diff_3_wrt_s_buff_.array().rowwise() *
          Eigen::pow(s_diff_3_wrt_tau_buff_.transpose().array(), 3) +
      q_diff_2_wrt_s_buff_.array().rowwise() *
          s_diff_1_wrt_tau_buff_.transpose().array() *
          s_diff_2_wrt_tau_buff_.transpose().array() * 3 +
      q_diff_1_wrt_s_buff_.array().rowwise() *
          s_diff_3_wrt_tau_buff_.transpose().array();
  q_diff_3_wrt_t_buff_ /= std::pow(Ts_, 3);
}
