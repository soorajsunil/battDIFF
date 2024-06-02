function dy = numDiff(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dy - derivative
% y  - input points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical differentiation of N uniformly spaced points.
% writer: S. Sunil (sunil11@uwindsor.ca) 
% last commit on 02 June 2024 by S. Sunil 
% Implemented based on: 
% [1] Numerical Methods for Engineers Book
% by Raymond P. Canale and Steven C. Chapra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

N  = length(y);
dy = nan(size(y));
h  = 1/(N-1); % step size

% High-Accuracy Numerical Differentiation Formulas
dy(1)   = forwardDiff(y,h,1);
dy(N)   = backwardDiff(y,h,N);
dy(2)   = centeredDiff(y,h,2);
dy(N-1) = centeredDiff(y,h,N-1);

for n=3:N-2
    dy(n) = centeredDiff2(y,h,n);
end

end

function fdash = forwardDiff(yn,h,n)
% The forward difference of accuracy O(h^2)
y   = yn(n);
yp1 = yn(n+1);
yp2 = yn(n+2);
fdash = (-yp2 + 4*yp1 - 3*y)/(2*h);
end

function fdash = backwardDiff(yn,h,n)
% The backward difference of accuracy O(h^2)
y   = yn(n);
yp1 = yn(n-1);
yp2 = yn(n-2);
fdash = (3*y - 4*yp1 + yp2)/(2*h);
end

function fdash = centeredDiff(yn,h,n)
% The centered difference of accuracy O(h^2)
ym1 = yn(n-1);
yp1 = yn(n+1);
fdash = (yp1-ym1)/(2*h);
end

function fdash = centeredDiff2(yn,h,n)
% The centered difference of accuracy O(h^4)
ym1 = yn(n-1);
ym2 = yn(n-2);
yp1 = yn(n+1);
yp2 = yn(n+2);
fdash = (-yp2 + 8*yp1 - 8*ym1 + ym2)/(12*h);
end