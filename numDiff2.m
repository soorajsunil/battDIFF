function [dyNewton,dyLagrange] = numDiff2(xn,xm,ym)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical differentiation of unevenly spaced points (nonequispaced data)
% writer: S. Sunil (sunil11@uwindsor.ca) 
% last commit on 02 June 2024 by S. Sunil 
% Implemented based on: 
% [1] Numerical Methods for Engineers Book
% by Raymond P. Canale and Steven C. Chapra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% dynewton   - Netwon method 
% dyLagrange - Lagrange method
% xm  - table x points of length M
% ym  - table y points of length M
% xn  - input x points of length N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Newton method (piecewise constant) % (bad approach)
slope = diff(ym)./diff(xm);
dyNewton = nan(size(xn));
for m = 1:length(xm)-1
    dyNewton(xn>=xm(m) & xn<=xm(m+1)) = slope(m);
end

% % % Lagrange interpolating polynomial method
dyLagrange = nan(size(xn));
for i=1:length(xn)
    [x0,x1,x2,y0,y1,y2] = pickClosestPoints(xm,ym,xn(i));
    dyLagrange(i) = lagrangeDifferentiation(x0,x1,x2,y0,y1,y2,xn(i));
end

end

function [x0,x1,x2,y0,y1,y2] = pickClosestPoints(xm,ym,x)

M = length(xm);
[~,idx] = min(abs(xm-x));

if idx ==1
    x0 = xm(1);
    x1 = xm(2);
    x2 = xm(end);
    y0 = ym(1);
    y1 = ym(2);
    y2 = ym(end);
elseif idx == M
    x0 = xm(M-2);
    x1 = xm(M-1);
    x2 = xm(M);
    y0 = ym(M-2);
    y1 = ym(M-1);
    y2 = ym(M);
else
    x0 = xm(idx-1);
    x1 = xm(idx);
    x2 = xm(idx+1);
    y0 = ym(idx-1);
    y1 = ym(idx);
    y2 = ym(idx+1);
end
end

function f2 = lagrangeDifferentiation(x0,x1,x2,y0,y1,y2,x)
f2 = ((2*x-x1-x2)/((x0-x1)*(x0-x2)))*y0 + ...
    ((2*x-x0-x2)/((x1-x0)*(x1-x2)))*y1 + ...
    ((2*x-x0-x1)/((x2-x0)*(x2-x1)))*y2;
end