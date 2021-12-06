#include <Eigen/Core>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionBase.hpp>
#ifndef PARAMETRIZATION
#define PARAMETRIZATION

namespace opstop {

gsplines::functions::FunctionExpression get_diffeo(double _ti, double _Ts,
                                                   double _sf);

gsplines::functions::CanonicPolynomial
get_diffeo_wrt_tau(double _ti, double _Ts, double _sf);

gsplines::functions::FunctionExpression get_tau_inv(double _ti, double _Ts);
gsplines::functions::FunctionExpression get_tau(double _ti, double _Ts);
void eval_pol_deg_5(const Eigen::VectorXd &_tau, const double _coeff[6],
                    Eigen::VectorXd &_result);

class ParametrizedCurveHelper {
public:
  std::size_t nglp_;

  double s_bar_deriv_coeff[4][6];

  double s_bar_deriv_coeff_diff_wrt_sf[4][6];

  double s_bar_deriv_coeff_diff_wrt_Ts[4][6];

  Eigen::VectorXd glp_;

  double ti_;

  std::unique_ptr<gsplines::functions::FunctionBase> position_;
  std::unique_ptr<gsplines::functions::FunctionBase> velocity_;
  std::unique_ptr<gsplines::functions::FunctionBase> acceleration_;
  std::unique_ptr<gsplines::functions::FunctionBase> jerk_;
  std::unique_ptr<gsplines::functions::FunctionBase> snap_;

  mutable Eigen::VectorXd s_val_buff_;
  mutable Eigen::VectorXd s_diff_1_wrt_tau_buff_;
  mutable Eigen::VectorXd s_diff_2_wrt_tau_buff_;
  mutable Eigen::VectorXd s_diff_3_wrt_tau_buff_;

  mutable Eigen::VectorXd s_val_diff_wrt_Ts_;
  mutable Eigen::VectorXd s_diff_1_wrt_tau_diff_wrt_Ts_;
  mutable Eigen::VectorXd s_diff_2_wrt_tau_diff_wrt_Ts_;
  mutable Eigen::VectorXd s_diff_3_wrt_tau_diff_wrt_Ts_;

  mutable Eigen::VectorXd s_val_diff_wrt_sf_;
  mutable Eigen::VectorXd s_diff_1_wrt_tau_diff_wrt_sf_;
  mutable Eigen::VectorXd s_diff_2_wrt_tau_diff_wrt_sf_;
  mutable Eigen::VectorXd s_diff_3_wrt_tau_diff_wrt_sf_;

  mutable Eigen::MatrixXd q_val_buff_;
  mutable Eigen::MatrixXd q_diff_1_wrt_s_buff_;
  mutable Eigen::MatrixXd q_diff_2_wrt_s_buff_;
  mutable Eigen::MatrixXd q_diff_3_wrt_s_buff_;
  mutable Eigen::MatrixXd q_diff_4_wrt_s_buff_;

  mutable Eigen::MatrixXd q_diff_1_wrt_t_buff_;
  mutable Eigen::MatrixXd q_diff_2_wrt_t_buff_;
  mutable Eigen::MatrixXd q_diff_3_wrt_t_buff_;
  mutable Eigen::MatrixXd q_diff_4_wrt_t_buff_;

  mutable Eigen::MatrixXd q_diff_wrt_Ts_buff_;
  mutable Eigen::MatrixXd q_diff_1_wrt_t_diff_wrt_Ts_buff_;
  mutable Eigen::MatrixXd q_diff_2_wrt_t_diff_wrt_Ts_buff_;
  mutable Eigen::MatrixXd q_diff_3_wrt_t_diff_wrt_Ts_buff_;
  mutable Eigen::MatrixXd q_diff_4_wrt_t_diff_wrt_Ts_buff_;

  mutable Eigen::MatrixXd q_diff_wrt_sf_buff_;
  mutable Eigen::MatrixXd q_diff_1_wrt_t_diff_wrt_sf_buff_;
  mutable Eigen::MatrixXd q_diff_2_wrt_t_diff_wrt_sf_buff_;
  mutable Eigen::MatrixXd q_diff_3_wrt_t_diff_wrt_sf_buff_;
  mutable Eigen::MatrixXd q_diff_4_wrt_t_diff_wrt_sf_buff_;

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
  ParametrizedCurveHelper(const gsplines::functions::FunctionBase &_curve,
                          std::size_t _nglp, double _ti);

  void compute_speed_at_gl_points(Eigen::VectorXd &_result);

  void set_diffeo(double _Ts, double _si);

  void compute_q_and_its_derivatives_wrt_s();
  void compute_s_and_its_derivatives_wrt_tau();
  void compute_s_and_its_derivatives_wrt_Ts_sf();
  void compute_q_and_its_derivatives_wrt_t();

  const Eigen::MatrixXd &compute_q_partial_diff_wrt_Ts();
  const Eigen::MatrixXd &compute_q_diff_1_wrt_t_partial_diff_wrt_Ts();
  const Eigen::MatrixXd &compute_q_diff_2_wrt_t_partial_diff_wrt_Ts();
  const Eigen::MatrixXd &compute_q_diff_3_wrt_t_partial_diff_wrt_Ts();

  const Eigen::MatrixXd &compute_q_partial_diff_wrt_sf();
  const Eigen::MatrixXd &compute_q_diff_1_wrt_t_partial_diff_wrt_sf();
  const Eigen::MatrixXd &compute_q_diff_2_wrt_t_partial_diff_wrt_sf();
  const Eigen::MatrixXd &compute_q_diff_3_wrt_t_partial_diff_wrt_sf();
};
} // namespace opstop
#endif /* ifndef PARAMETRIZATION*/
