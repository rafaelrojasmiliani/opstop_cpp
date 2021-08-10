#ifndef DIFFEO_CONSTRAINTS_H
#define DIFFEO_CONSTRAINTS_H
#include <ifopt/constraint_set.h>
#include <opstop/parametrization.hpp>

namespace opstop {
class DiffeoConstraints : public ifopt::ConstraintSet {
private:
  ifopt::Component::VecBound bounds_vector_;

  const double ti_;
  const double exec_time_;

  mutable Eigen::VectorXd value_;

public:
  DiffeoConstraints(double _ti, double _exec_time);
  Eigen::VectorXd GetValues() const override;

  ifopt::Component::VecBound GetBounds() const override;

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override;

  virtual ~DiffeoConstraints() = default;
};

} // namespace opstop
#endif /* DIFFEO_CONSTRAINTS_H */
