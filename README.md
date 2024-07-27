# Numerical-Differentiation
This repo contains two Matlab functions to compute the numerical differentiation of a vector.

Numerical differentiation of N uniformly spaced points implemented based on high-accuracy numerical differentiation formulas from [1]. 
>  `dy = numDiff(y)`
> 
where:
- `y` is a column vector of uniformly spaced points of length N
- `dy` is the numerical derivative of y


Numerical differentiation of unevenly spaced points (non-equispaced data) implemented based on [1]. 
> `[dyNewton,dyLagrange] = numDiff2(xn,xm,ym)`
> 
where:
- `xn` is the query vector x of length N
- `xm` is the breakpoint (or knot) vector x of length M
- `ym` is the breakpoint (or knot) vector y of length M
- `dyNewton` is the derivative vector of length N computed using the piecewise constant approximation method (Newton) 
- `dyLagrange` is the derivative vector of length N computed using the Lagrange interpolating polynomials method

## Reference: 

[1] Numerical Methods for Engineers Book by Raymond P. Canale and Steven C. Chapra. 
