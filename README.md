# Numerical-Differentiation

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
- `xn` is input query vector x points of length N
- `xm` is the table vector of x points of length M
- `ym` is the table vector of y points of length M
- `dyNewton`  is the derivative of length N based on the Newton method (piecewise constant) 
- `dyLagrange` is the derivative of length N based on the Lagrange interpolating polynomials method 

## Reference: 

[1] Numerical Methods for Engineers Book by Raymond P. Canale and Steven C. Chapra. 
