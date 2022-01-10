function [x,err,k] = StormSpr(prob,opt)
maxiter = opt.maxiter;
xk = opt.x0;
sk = opt.sk;
s0 = opt.s0;
gamma = opt.gamma;
[m,n] = size(prob.A);
delta = opt.delta;
err = Inf;
for k = 1:maxiter
    % each iteration, extract a submatrix
    ind = sort(randsample(m,floor(gamma*m)));
    subA = prob.A(ind,:);
    subdata = prob.data(ind);
    suby = sign(subA*xk).*subdata;
    
    % solve sparsity-constrained least-squares
    myfunc = @(x) myfunc2(x,subA,suby);
    pars.tol = 1e-6;
    pars.iteron = 0;
    pars.maxit = 500;
    if err < 0.2
        sk = s0;
        out = IIHT(prob.d1,sk,myfunc,pars);
    else
        sk = 2*s0;
        out = IIHT(prob.d1,sk,myfunc,pars);
    end
    xkt = out.x;
%     xkt = BIHT(subA,sk,suby,0.01);
    diffx = norm(xkt-xk);
%     err = min(norm(prob.x0-xk),norm(prob.x0+xk))/norm(prob.x0);
    [err, x00] = computerelerror(xk, prob.x0);
    f_value = abs(prob.A*x00);
    if opt.verbosity
        fprintf('Iter %d error--%.3f \n',k, err);
    end
    if diffx < delta || err < delta || norm(f_value-prob.data)/norm(prob.data) < 1e-3
        break;
    end
    xk = xkt;   
end
% x = truncate(xk,s0);
x = truncate(x00,s0);
end