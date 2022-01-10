function [ rel, x00 ] = computerelerror(x,x0)
%COMPUTERELERROR compute the relative error to trace the convergence
%   Detailed explanation goes here

% vectorize the input
x = x(:);
x0 = x0(:);
innertemp = dot(x0,x);
c = innertemp/abs(innertemp);
rel = norm(c*x0-x)/norm(x0);
x00 = conj(c)*x;
end