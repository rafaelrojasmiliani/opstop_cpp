
#ifndef VELOCITY_CONSTRAINT_H
#define VELOCITY_CONSTRAINT_H
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <ifopt/constraint_set.h>
#include <opstop/parametrization.hpp>

class VelocityConstraints : public ifopt::ConstraintSet {
private:
  mutable ParametrizedCurveHelper helper_;
  mutable Eigen::VectorXd value_buff_;
  ifopt::Component::VecBound bounds_vector_;

public:
  VelocityConstraints(const gsplines::functions::FunctionExpression &_curve,
                      std::size_t _nglp, double _ti,
                      std::vector<double> &_bound);

  virtual ~VelocityConstraints() = default;

  Eigen::VectorXd GetValues() const override;

  ifopt::Component::VecBound GetBounds() const override;

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override;

  Eigen::VectorXd __GetValues(Eigen::Vector2d &_x) const;

  Eigen::MatrixXd __GetJacobian(Eigen::Vector2d &_x) const;
};

#endif /* VELOCITY_CONSTRAINT_H */
