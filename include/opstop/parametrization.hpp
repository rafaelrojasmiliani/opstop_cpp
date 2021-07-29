#include <Eigen/Core>
#include <gsplines++/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
#ifndef PARAMETRIZATION
#define PARAMETRIZATION

void eval_pol_deg_5(const Eigen::VectorXd &_tau, const double _coeff[6],
                    Eigen::VectorXd &_result) {

  _result.setConstant(_coeff[5]);

  for (int i = 4; i >= 0; i--) {
    _result.array() *= _tau.array();

    _result.array() += _coeff[i];
  }
}

class ParametrizedCurveHelper {
private:
  std::size_t nglp_;

  double s_bar_deriv_coeff[4][6];

  double s_bar_deriv_coeff_diff_wrt_si[4][6];

  double s_bar_deriv_coeff_diff_wrt_Ts[4][6];

  Eigen::VectorXd glp_;

  double ti_;

  std::unique_ptr<gsplines::functions::FunctionExpression> position_;
  std::unique_ptr<gsplines::functions::FunctionExpression> velocity_;
  std::unique_ptr<gsplines::functions::FunctionExpression> acceleration_;
  std::unique_ptr<gsplines::functions::FunctionExpression> jerk_;
  std::unique_ptr<gsplines::functions::FunctionExpression> snap_;

  mutable Eigen::VectorXd vector_buff_1_;
  mutable Eigen::VectorXd vector_buff_2_;

  mutable Eigen::VectorXd s_val_buff_;
  mutable Eigen::VectorXd s_diff_1_wrt_tau_buff_;
  mutable Eigen::VectorXd s_diff_2_wrt_tau_buff_;
  mutable Eigen::VectorXd s_diff_3_wrt_tau_buff_;

  mutable Eigen::VectorXd q_val_buff_;
  mutable Eigen::VectorXd q_diff_1_wrt_s_buff_;
  mutable Eigen::VectorXd q_diff_2_wrt_s_buff_;
  mutable Eigen::VectorXd q_diff_3_wrt_s_buff_;
  mutable Eigen::VectorXd q_diff_4_wrt_s_buff_;

  mutable Eigen::VectorXd q_diff_1_wrt_t_buff_;
  mutable Eigen::VectorXd q_diff_2_wrt_t_buff_;
  mutable Eigen::VectorXd q_diff_3_wrt_t_buff_;
  mutable Eigen::VectorXd q_diff_4_wrt_t_buff_;

  mutable double Ts_;
  mutable double sf_;

public:
  /*
  const Eigen::VectorXd &get_domain_buff() const { return domain_buff_; }
  const Eigen::MatrixXd &get_position_buff() const { return position_buff_; }
  const Eigen::MatrixXd &get_velocity_buff() const { return velocity_buff_; }
  const Eigen::MatrixXd &get_acceleration_buff() const {
    return acceleration_buff_;
  }
  const Eigen::MatrixXd &get_jerk_buff() const { return jerk_buff_; }
  const Eigen::MatrixXd &get_snap_buff() const { return snap_buff_; }
*/
  ParametrizedCurveHelper(
      std::unique_ptr<const gsplines::functions::FunctionExpression> _curve,
      std::size_t _nglp, double _ti);

  void compute_speed_at_gl_points(Eigen::VectorXd &_result);

  void set_diffeo(double _Ts, double _si);
  void compute_domain_points();
  void compute_position();
  void compute_velocity();
  void compute_acceleration();
  void compute_jerk();
  void compute_snap();
  void compute_q_and_its_derivatives_wrt_s();
  void compute_s_and_its_derivatives_wrt_tau();
  void compute_q_and_its_derivatives_wrt_t();
};
#endif /* ifndef PARAMETRIZATION*/
