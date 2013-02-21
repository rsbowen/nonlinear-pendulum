nonlinear-pendulum
==================
This is, for practice, a non-linear pendulum PDE solver.  We are discretizing

\theta''(t)= - sin(\theta(t))

to get

(1/h^2) (\theta\_(i-1) - 2 \theta\_i + \theta\_(i+1)) + sin(\theta\_i) = 0

which is nonlinear.  I'll use a Newton update step.
