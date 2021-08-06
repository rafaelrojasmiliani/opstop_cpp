#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <opstop/parametrization.hpp>

void eval_pol_deg_5(const Eigen::VectorXd &_tau, const double _coeff[6],
                    Eigen::VectorXd &_result) {

  _result.setConstant(_coeff[5]);

  for (int i = 4; i >= 0; i--) {
    _result.array() *= _tau.array();

    _result.array() += _coeff[i];
  }
}
ParametrizedCurveHelper::ParametrizedCurveHelper(
    const gsplines::functions::FunctionExpression &_curve, std::size_t _nglp,
    double _ti)
    : position_(_curve.clone()), velocity_(_curve.deriv()),
      acceleration_(velocity_->deriv()), jerk_(acceleration_->deriv()),
      snap_(jerk_->deriv()),
      /* Parametrization and its derivatives wrt tau */
      s_val_buff_(_nglp), s_diff_1_wrt_tau_buff_(_nglp),
      s_diff_2_wrt_tau_buff_(_nglp), s_diff_3_wrt_tau_buff_(_nglp),
      /* Derivatives of the Parametrization and its derivatives wrt tau wrt Ts
       */
      s_val_diff_wrt_Ts_(_nglp), s_diff_1_wrt_tau_diff_wrt_Ts_(_nglp),
      s_diff_2_wrt_tau_diff_wrt_Ts_(_nglp),
      s_diff_3_wrt_tau_diff_wrt_Ts_(_nglp),
      /* Derivatives of the Parametrization and its derivatives wrt tau wrt sf
       */
      s_val_diff_wrt_sf_(_nglp), s_diff_1_wrt_tau_diff_wrt_sf_(_nglp),
      s_diff_2_wrt_tau_diff_wrt_sf_(_nglp),
      s_diff_3_wrt_tau_diff_wrt_sf_(_nglp),
      /* Values of the trajectory and its derivatives wrt the nominal time*/
      q_val_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_1_wrt_s_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_2_wrt_s_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_3_wrt_s_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_4_wrt_s_buff_(_nglp, _curve.get_codom_dim()),
      /* Values of the trajectory and its derivatives wrt time */
      q_diff_1_wrt_t_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_2_wrt_t_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_3_wrt_t_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_4_wrt_t_buff_(_nglp, _curve.get_codom_dim()),
      /* Derivatives of the trajectory and its derivatives wrt time wrt Ts*/
      q_diff_wrt_Ts_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_1_wrt_t_diff_wrt_Ts_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_2_wrt_t_diff_wrt_Ts_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_3_wrt_t_diff_wrt_Ts_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_4_wrt_t_diff_wrt_Ts_buff_(_nglp, _curve.get_codom_dim()),
      /* Derivatives of the trajectory and its derivatives wrt time wrt sf*/
      q_diff_wrt_sf_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_1_wrt_t_diff_wrt_sf_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_2_wrt_t_diff_wrt_sf_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_3_wrt_t_diff_wrt_sf_buff_(_nglp, _curve.get_codom_dim()),
      q_diff_4_wrt_t_diff_wrt_sf_buff_(_nglp, _curve.get_codom_dim()), ti_(_ti),
      nglp_(_nglp)

{

  for (std::size_t coeff = 0; coeff < 5; coeff++) {

    for (std::size_t deriv = 0; deriv < 4; deriv++) {

      s_bar_deriv_coeff[deriv][coeff] = 0.0;

      s_bar_deriv_coeff_diff_wrt_sf[deriv][coeff] = 0.0;

      s_bar_deriv_coeff_diff_wrt_Ts[deriv][coeff] = 0.0;
    }
  }

  std::tie(glp_, std::ignore) =
      gsplines::collocation::legendre_gauss_lobatto_points_and_weights(_nglp);

  glp_.array() += 1.0;
  glp_ /= 2.0;

  s_bar_deriv_coeff_diff_wrt_sf[0][0] = 0.0;
  s_bar_deriv_coeff_diff_wrt_sf[0][1] = 0.0;
  s_bar_deriv_coeff_diff_wrt_sf[0][3] = 10.0;
  s_bar_deriv_coeff_diff_wrt_sf[0][4] = -15.0;
  s_bar_deriv_coeff_diff_wrt_sf[0][5] = 6.0;

  s_bar_deriv_coeff_diff_wrt_Ts[0][0] = 0.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][1] = 1.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][3] = -6.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][4] = 8.0;
  s_bar_deriv_coeff_diff_wrt_Ts[0][5] = -3.0;

  for (std::size_t deriv = 1; deriv < 4; deriv++) {
    for (std::size_t coeff = 0; coeff < 6 - deriv; coeff++) {
      s_bar_deriv_coeff_diff_wrt_sf[deriv][coeff] =
          s_bar_deriv_coeff_diff_wrt_sf[deriv - 1][coeff + 1] * (coeff + 1);

      s_bar_deriv_coeff_diff_wrt_Ts[deriv][coeff] =
          s_bar_deriv_coeff_diff_wrt_Ts[deriv - 1][coeff + 1] * (coeff + 1);
    }
  }
}

void ParametrizedCurveHelper::set_diffeo(double _Ts, double _sf) {

  s_bar_deriv_coeff[0][0] = ti_;
  s_bar_deriv_coeff[0][1] = _Ts;
  s_bar_deriv_coeff[0][3] = -6.0 * _Ts + 10.0 * _sf - 10.0 * ti_;
  s_bar_deriv_coeff[0][4] = 8.0 * _Ts - 15.0 * _sf + 15.0 * ti_;
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

void ParametrizedCurveHelper::compute_s_and_its_derivatives_wrt_Ts_sf() {

  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_Ts[0], s_val_diff_wrt_Ts_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_Ts[1],
                 s_diff_1_wrt_tau_diff_wrt_Ts_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_Ts[2],
                 s_diff_2_wrt_tau_diff_wrt_Ts_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_Ts[3],
                 s_diff_3_wrt_tau_diff_wrt_Ts_);

  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_sf[0], s_val_diff_wrt_sf_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_sf[1],
                 s_diff_1_wrt_tau_diff_wrt_sf_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_sf[2],
                 s_diff_2_wrt_tau_diff_wrt_sf_);
  eval_pol_deg_5(glp_, s_bar_deriv_coeff_diff_wrt_sf[3],
                 s_diff_3_wrt_tau_diff_wrt_sf_);
}

void ParametrizedCurveHelper::compute_q_and_its_derivatives_wrt_t() {

  q_diff_1_wrt_t_buff_ = q_diff_1_wrt_s_buff_.array().colwise() *
                         s_diff_1_wrt_tau_buff_.array() / Ts_;
  // ---- Second Derivative -----
  q_diff_2_wrt_t_buff_ =
      q_diff_2_wrt_s_buff_.array().colwise() *
          Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 2) +
      q_diff_1_wrt_s_buff_.array().colwise() * s_diff_2_wrt_tau_buff_.array();

  q_diff_2_wrt_t_buff_ /= std::pow(Ts_, 2);

  // ---- Third Derivative -----
  q_diff_3_wrt_t_buff_ =
      q_diff_3_wrt_s_buff_.array().colwise() *
          Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 3) +
      q_diff_2_wrt_s_buff_.array().colwise() *
          (s_diff_1_wrt_tau_buff_.array() * s_diff_2_wrt_tau_buff_.array() *
           3.0) +
      q_diff_1_wrt_s_buff_.array().colwise() * s_diff_3_wrt_tau_buff_.array();
  q_diff_3_wrt_t_buff_ /= std::pow(Ts_, 3);
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_partial_diff_wrt_Ts() {
  q_diff_wrt_Ts_buff_ =
      q_diff_1_wrt_s_buff_.array().colwise() * s_val_diff_wrt_Ts_.array();
  return q_diff_wrt_Ts_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_1_wrt_t_partial_diff_wrt_Ts() {

  q_diff_1_wrt_t_diff_wrt_Ts_buff_.array() =
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() * s_diff_1_wrt_tau_buff_.array());

  q_diff_1_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_1_wrt_tau_diff_wrt_Ts_.array();

  q_diff_1_wrt_t_diff_wrt_Ts_buff_ /= Ts_;

  q_diff_1_wrt_t_diff_wrt_Ts_buff_.array() +=
      -q_diff_1_wrt_t_buff_.array() / Ts_;
  return q_diff_wrt_Ts_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_2_wrt_t_partial_diff_wrt_Ts() {

  q_diff_2_wrt_t_diff_wrt_Ts_buff_ =
      q_diff_3_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() *
       Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 2));

  q_diff_2_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_buff_.array() * s_diff_1_wrt_tau_diff_wrt_Ts_.array() *
       2);

  q_diff_2_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() * s_diff_2_wrt_tau_buff_.array());

  q_diff_2_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_2_wrt_tau_diff_wrt_Ts_.array();

  q_diff_2_wrt_t_diff_wrt_Ts_buff_ /= std::pow(Ts_, 2);
  q_diff_2_wrt_t_diff_wrt_Ts_buff_ += -2 * q_diff_2_wrt_t_buff_ / Ts_;

  return q_diff_2_wrt_t_diff_wrt_Ts_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_3_wrt_t_partial_diff_wrt_Ts() {

  q_diff_3_wrt_t_diff_wrt_Ts_buff_ =
      q_diff_4_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() *
       Eigen::pow(s_diff_1_wrt_tau_diff_wrt_Ts_.array(), 3));

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_3_wrt_s_buff_.array().colwise() *
      (Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 2) *
       s_diff_1_wrt_tau_diff_wrt_Ts_.array() * 3);

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_3_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() * s_diff_1_wrt_tau_buff_.array() *
       s_diff_2_wrt_tau_buff_.array() * 3);

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_diff_wrt_Ts_.array() * s_diff_2_wrt_tau_buff_.array() *
       3);

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_buff_.array() * s_diff_2_wrt_tau_diff_wrt_Ts_.array() *
       3);

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_Ts_.array() * s_diff_3_wrt_tau_buff_.array());

  q_diff_3_wrt_t_diff_wrt_Ts_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_3_wrt_tau_diff_wrt_Ts_.array();

  q_diff_3_wrt_t_diff_wrt_Ts_buff_ /= std::pow(Ts_, 3);

  q_diff_3_wrt_t_diff_wrt_Ts_buff_ += -3 * q_diff_3_wrt_t_buff_ / Ts_;

  return q_diff_3_wrt_t_diff_wrt_Ts_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_partial_diff_wrt_sf() {
  q_diff_wrt_sf_buff_ =
      q_diff_1_wrt_s_buff_.array().colwise() * s_val_diff_wrt_sf_.array();
  return q_diff_wrt_sf_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_1_wrt_t_partial_diff_wrt_sf() {

  q_diff_1_wrt_t_diff_wrt_sf_buff_.array() =
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() * s_diff_1_wrt_tau_buff_.array());

  q_diff_1_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_1_wrt_tau_diff_wrt_sf_.array();

  q_diff_1_wrt_t_diff_wrt_sf_buff_ /= Ts_;

  return q_diff_wrt_sf_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_2_wrt_t_partial_diff_wrt_sf() {

  q_diff_2_wrt_t_diff_wrt_sf_buff_ =
      q_diff_3_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() *
       Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 2));

  q_diff_2_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_buff_.array() * s_diff_1_wrt_tau_diff_wrt_sf_.array() *
       2);

  q_diff_2_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() * s_diff_2_wrt_tau_buff_.array());

  q_diff_2_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_2_wrt_tau_diff_wrt_sf_.array();

  q_diff_2_wrt_t_diff_wrt_sf_buff_ /= std::pow(Ts_, 2);

  return q_diff_2_wrt_t_diff_wrt_sf_buff_;
}

const Eigen::MatrixXd &
ParametrizedCurveHelper::compute_q_diff_3_wrt_t_partial_diff_wrt_sf() {

  q_diff_3_wrt_t_diff_wrt_sf_buff_ =
      q_diff_4_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() *
       Eigen::pow(s_diff_1_wrt_tau_diff_wrt_sf_.array(), 3));

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_3_wrt_s_buff_.array().colwise() *
      (Eigen::pow(s_diff_1_wrt_tau_buff_.array(), 2) *
       s_diff_1_wrt_tau_diff_wrt_sf_.array() * 3);

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_3_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() * s_diff_1_wrt_tau_buff_.array() *
       s_diff_2_wrt_tau_buff_.array() * 3);

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_diff_wrt_sf_.array() * s_diff_2_wrt_tau_buff_.array() *
       3);

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_diff_1_wrt_tau_buff_.array() * s_diff_2_wrt_tau_diff_wrt_sf_.array() *
       3);

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_2_wrt_s_buff_.array().colwise() *
      (s_val_diff_wrt_sf_.array() * s_diff_3_wrt_tau_buff_.array());

  q_diff_3_wrt_t_diff_wrt_sf_buff_.array() +=
      q_diff_1_wrt_s_buff_.array().colwise() *
      s_diff_3_wrt_tau_diff_wrt_sf_.array();

  q_diff_3_wrt_t_diff_wrt_sf_buff_ /= std::pow(Ts_, 3);

  return q_diff_3_wrt_t_diff_wrt_sf_buff_;
}

gsplines::functions::FunctionExpression get_diffeo(double _ti, double _Ts,
                                                   double _sf) {

  Eigen::VectorXd pol_coeff(6);
  pol_coeff << _ti, _Ts, 0.0, -6.0 * _Ts + 10.0 * _sf - 10.0 * _ti,
      8.0 * _Ts - 15.0 * _sf + 15.0 * _ti, -3.0 * _Ts + 6.0 * _sf - 6.0 * _ti;

  gsplines::functions::CanonicPolynomial pol({0, 1}, pol_coeff);

  gsplines::functions::FunctionExpression tau_exp =
      (1.0 / _Ts) *
      (gsplines::functions::Identity({_ti, _ti + _Ts}) +
       gsplines::functions::ConstFunction({_ti, _ti + _Ts}, 1, -_ti));

  return gsplines::functions::Identity({0, _ti}).concat(pol.compose(tau_exp));
}
