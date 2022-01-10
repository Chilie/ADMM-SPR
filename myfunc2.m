function [f,varargout] = myfunc2(x,A,b)
[m,n] = size(A);
 f = 0.5*norm(A*x-b)^2;
if nargout > 1
    varargout{1} = A'*(A*x-b)/m;
end
end