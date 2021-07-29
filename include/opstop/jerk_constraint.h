#ifndef JERK_CONSTRAINT
#define JERK_CONSTRAINT
#include <gsplines++/Functions/Function.hpp>
#include <gsplines++/PiecewiseFunction.hpp>
#include <ifopt/constraint_set.h>

class JerkConstraint : public ifopt::ConstraintSet {
public:
  JerkConstraint(gsplines::PiecewiseFunction &_motion,
                 double _start_to_stop_time, double _bound,
                 std::size_t _number_gl_points);

  Eigen::VectorXd GetValues() const override;
  ifopt::Component::VecBound GetBounds() const override;

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override;

  virtual ~JerkConstraint();

private:
  Eigen::VectorXd values_;
  double exec_time_;
  double start_to_stop_time_;
  ifopt::Component::VecBound bounds_;
  const std::unique_ptr<gsplines::functions::FunctionExpression> motion_;
  const std::unique_ptr<gsplines::functions::FunctionExpression>
      motion_velocity_;
  const std::unique_ptr<gsplines::functions::FunctionExpression>
      motion_acceleration_;
  const std::unique_ptr<gsplines::functions::FunctionExpression> motion_jerk_;
  const std::unique_ptr<gsplines::functions::FunctionExpression> motion_snap_;
  std::size_t number_of_gl_points_;
  const std::pair<const Eigen::VectorXd, const Eigen::VectorXd>
      gauss_lobatt_points_weights_;
};
#endif /* ifndef JERK_CONSTRAINT */
