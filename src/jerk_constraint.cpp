
#include <opstop/jerk_constraint.h>

JerkConstraint::JerkConstraint(gsplines::PiecewiseFunction &_motion,
                               double _start_to_stop_time, double _bound,
                               std::size_t _number_gl_points)
    : ConstraintSet(_number_gl_points, "JerkConstraint"),
      motion_(_motion.clone()), motion_velocity_(_motion.deriv()),
      motion_acceleration_(_motion.deriv(2)), motion_snap_(_motion.deriv(3)),
      start_to_stop_time_(_start_to_stop_time),
      number_of_gl_points_(_number_gl_points) {

  ifopt::Bounds default_bound(-_bound, _bound);
  bounds_ = ifopt::Component::VecBound(_number_gl_points, default_bound);
}

Eigen::VectorXd JerkConstraint::GetValues() const {

  // static Eigen::VectorXd tauv(GetRows());
  static Eigen::VectorXd result(1);

  // tauv = GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();
  // result(0) = tauv.sum();
  return result;
}
ifopt::Component::VecBound JerkConstraint::GetBounds() const { return bounds_; }

void JerkConstraint::FillJacobianBlock(std::string _set_name,
                                       Jacobian &_jac_block) const {
  // for (unsigned int i = 0;
  //    i < GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows();
  //  i++)
  // _jac_block.coeffRef(0, i) = 1.0;
}

JerkConstraint::~JerkConstraint() {}
