function [grad] = mygrad2(x,A,b)
% m = length(y);
% temp = fft(y);
% tmpgrad = ifft(conj(temp).*(temp.*fft(x)-b));
grad = A'*(A*x-b);
end