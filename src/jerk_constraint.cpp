
#include <gsplines++/Collocation/GaussLobattoPointsWeights.hpp>
#include <opstop/jerk_constraint.h>
namespace opstop {

JerkConstraint::JerkConstraint(gsplines::PiecewiseFunction &_motion,
                               double _start_to_stop_time, double _bound,
                               std::size_t _number_gl_points)
    : ConstraintSet(_number_gl_points, "JerkConstraint"),
      motion_(_motion.clone()), motion_velocity_(_motion.deriv()),
      motion_acceleration_(_motion.deriv(2)), motion_jerk_(_motion.deriv(3)),
      motion_snap_(_motion.deriv(4)), start_to_stop_time_(_start_to_stop_time),
      number_of_gl_points_(_number_gl_points),
      gauss_lobatt_points_weights_(
          gsplines::collocation::legendre_gauss_lobatto_points_and_weights(
              _number_gl_points)) {

  ifopt::Bounds default_bound(-_bound, _bound);
  bounds_ = ifopt::Component::VecBound(_number_gl_points, default_bound);
}

Eigen::VectorXd JerkConstraint::GetValues() const {

  tauv = GetVariables()->GetComponent("Tssf")->GetValues();
  result(0) = tauv.sum();

  return std::move(result);
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
} // namespace opstop
