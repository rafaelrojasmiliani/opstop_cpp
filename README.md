# Minimum-Time Path-Consistent Stopping Library for Collaborative Robots
##Overview
Smooth stopping of collaborative robots is crucial in environments where safety is essential.
However, abrupt stops can result in potential collisions; therefore, a collision-free stop must lie within the original collision-free path.
This computes a minimum-time stopping trajectory along a collision-free path. It is based on the [gsplines library](https://github.com/rafaelrojasmiliani/gsplines_cpp).

The core feature of this library is the formulation of a time minimization problem for a specific-kind of stopping parametrization of a curve (a [gspline](https://github.com/rafaelrojasmiliani/gsplines_cpp) curve) with smoothness constraints.

This means that when an emergency stop is initiated, the robot will decelerate smoothly along its original path rather than deviating.
In addition, smoothness and torque constrains ensure a gentle and a dynamically feasible stop.

This library is particularly useful in scenarios where robots work in close proximity to humans or other robots, and an immediate yet safe stop is essential due to an unexpected event or emergency.

# Publications

This library was used to publish

- Rojas, Rafael A., Andrea Giusti, and Renato Vidoni. "Online Computation of Time-Optimization-Based, Smooth and Path-Consistent Stop Trajectories for Robots." Robotics 11.4 (2022): 70. (**We published here primarily due to an impending deadline, though a more prestigious journal could have been an option under different circumstances**)
```
@article{rojas2022online,
  title={Online Computation of Time-Optimization-Based, Smooth and Path-Consistent Stop Trajectories for Robots},
  author={Rojas, Rafael A and Giusti, Andrea and Vidoni, Renato},
  journal={Robotics},
  volume={11},
  number={4},
  pages={70},
  year={2022},
  publisher={MDPI}
}
```
