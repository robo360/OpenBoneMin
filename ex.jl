using DifferentialEquations
using Plots
A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
tspan = (0.0,1.0)

function f(u,p,t)
      u0 = rand(4,2)
      u = A*u
 end
prob = ODEProblem(f, 0,tspan)
sol = solve(prob)
plot(sol)
