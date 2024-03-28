# Minimum-Time Path-Consistent Stopping Library for Collaborative Robots

Smooth stopping of collaborative robots is crucial in environments where safety is essential.
However, abrupt stops can result in potential collisions; therefore, a collision-free stop must lie within the original collision-free path.
This computes a minimum-time stopping trajectory along a collision-free path. It is based on the [gsplines library](https://github.com/rafaelrojasmiliani/gsplines_cpp).

The core feature of this library is the formulation of a time minimization problem for a specific-kind of stopping parametrization of a curve (a [gspline](https://github.com/rafaelrojasmiliani/gsplines_cpp) curve) with smoothness constraints.

This means that when an emergency stop is initiated, the robot will minimize the stopping time by decelerating smoothly along its original path rather than deviating.
In addition, smoothness and torque constrains ensure a gentle and a dynamically feasible stop.

This library is particularly useful in scenarios where robots work in proximity to humans or other robots, and an immediate yet safe stop is essential due to an unexpected event or emergency.

- ROS implementation [here](https://github.com/rafaelrojasmiliani/opstop_ros)
- Compatibility with moveit [here](https://github.com/rafaelrojasmiliani/gsplines_moveit)
- Contact: Rafael A. Rojas rafaelrojasmiliani@gmail.com.
- Docker containers with this library already installed
    - *vim awesome plugins for development and moveit* rafa606/moveit-opstop-vim-dev:noetic
    - *vim awesome plugins for development and awesome ros packages*: rafa606/ros-opstop-vim-dev:noetic
- **Remark on real-time stop** The time to compute the stopping trajectory depends strongly in the optimizer. Currently this repo uses ipopt. To have good performances (<10ms with a good pc) use the `ma27` linear solver.

# Example

```python
import gsplines
import gsplines.plot as gplot
import numpy as np
import opstop
from gsplines.optimization import minimum_jerk_path


model_file = 'path_to_urdf_robot_description'

dim = 7   # number of joints of the robot
# Generate a random numpy array of wayponts
# (each row is a waypoint in R^n)
number_of_waypoints = 5
waypoints = np.random.rand(number_of_waypoints, dim)
# The a minimum jerk trajectory with execution time of 5s
trj = minimum_jerk_path(waypoints).linear_scaling_new_execution_time(5.0)
# Select the time to stop as the 60% of the time.
stop_ti = trj.get_domain()[1]*0.6
# Get a parametrization that minimizes the time and  does not
# allow an increment in the acceleration larger than 50%
optimal_parametrization = opstop.minimum_time_bounded_acceleration(
    trj, stop_ti, 1.5, str(model_file), 8)
# Obtain the stopping trajectory
stop_trj = trj.compose(optimal_parametrization)
# Plot the nominal and the stopping trajectory
gplot.plot_compare([stop_trj, trj], ['red', 'blue'], [
                   'Emergency Stop Trajectory',
                   'Original Trajectory'], _show=True, _up_to_deriv=2)
```
This code will plot two trajectories. The blue is the original trajectory of the robot in the joint space. The red is the emergency stop trajectory that minimizes the time, avoid accelerations larger than 50% of the original trajectory and respect the torque constraints (defined in the urdf).
![alt text](img/plot.png)

# Installation

## Installation with ROS
First get the dependencies
```bash
sudo apt-get install libeigen3-dev ros-noetic-hpp-fcl robotpkg-pinocchio coinor-libipopt-dev ros-noetic-ifopt
```
Don't forget to export the `robotpkg` paths
```bash
export PATH=/opt/openrobots/bin:$PATH
export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/opt/openrobots/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=/opt/openrobots:$CMAKE_PREFIX_PATH
```

Install the gsplines

```bash
wget https://github.com/rafaelrojasmiliani/gsplines_cpp/releases/download/master/gsplines-0.0.1-amd64.deb
sudo dpkg -i gsplines-0.0.1-amd64.deb
```

Install this package

```bash
wget https://github.com/rafaelrojasmiliani/opstop_cpp/releases/download/master/opstop-0.0.1-gcc-11-amd64.deb
sudo dpkg -i opstop-0.0.1-amd64.deb
```

If you are working with pythond bindings and gcc-11 you might need the gsplines and optstop compiled with the same version

```bash
wget https://github.com/rafaelrojasmiliani/gsplines_cpp/releases/download/master/gsplines-0.0.1-gcc-11-amd64.deb
wget https://github.com/rafaelrojasmiliani/opstop_cpp/releases/download/master/opstop-0.0.1-gcc-11-amd64.deb
sudo dpkg -i gsplines-0.0.1-gcc-11-amd64.deb  opstop-0.0.1-gcc-11-amd64.deb
```

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
