#include <gsplines/Functions/FunctionExpression.hpp>
#include <ifopt/constraint_set.h>
#include <opstop/parametrization.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#ifndef TORQUE_CONSTRAINT_H
#define TORQUE_CONSTRAINT_H

namespace opstop {

class TorqueConstraint : public ifopt::ConstraintSet {
private:
  mutable pinocchio::Model model_;
  mutable pinocchio::Data data_;
  mutable ParametrizedCurveHelper helper_;
  mutable Eigen::VectorXd torque_buff_;

  ifopt::Component::VecBound bounds_vector_;

public:
  TorqueConstraint(const gsplines::functions::FunctionExpression &_curve,
                   std::size_t _nglp, double _ti,
                   const std::vector<double> &_bound,
                   const pinocchio::Model &_model);

  Eigen::VectorXd GetValues() const override;

  ifopt::Component::VecBound GetBounds() const override;

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override;

  virtual ~TorqueConstraint() {}

  Eigen::VectorXd __GetValues(Eigen::Vector2d &_x) const;

  Eigen::MatrixXd __GetJacobian(Eigen::Vector2d &_x) const;
};

} // namespace opstop
#endif /* TORQUE_CONSTRAINT_H */
