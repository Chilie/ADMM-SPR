function [x,err,k] = ADM2(prob,opt)
maxiter = opt.maxiter;
xk = opt.x0;
sk = opt.sk;
s0 = opt.s0;
[m,n] = size(prob.A);
lambda = zeros(m,1);
yk = prob.A*xk;
zk = yk;
rho = opt.rho;
err = Inf;
delta = opt.delta;
% nom2 = 1/(norm(prob.A).^2);
% tau = 0.01;%sqrt(nom2);
% sigma = nom2/tau;
% f = @(x) real(prob.data.*exp(1i*angle(x)));
% begin iteration
for k = 1:maxiter
%     for jj = 1:m
%         aa = prob.A(jj,:);
%         temp = aa*xk;
%         temp2 = (prob.data(jj)-abs(temp))*exp(1i*angle(temp));
%         xk = xk + (real(temp2)/norm(aa)^2)*aa';
%         xk = truncate(xk,2*sk);
%         fprintf('error--%.3f\n',min(norm(prob.x0-xk),norm(prob.x0+xk)));
%     end

% %% Linearized ADMM 
%       tempz = yk + tau*(prob.A*xk);
%     ykp = tempz - tau*f(tempz/tau);
%     tmp = xk - sigma*prob.A'*(2*ykp-yk);
% %     xk = truncate(xk - sigma*prob.A'*(2*ykp-yk),sk);
%     xk = sign(tmp).*max(abs(tmp)-sigma,0);
% yk = ykp;

%%

%     xk = BIHT(prob.A,sk,yk-lambda,0.01);
    tmp = prob.A*xk - lambda;
%     yk = real(prob.data.*exp(1i*angle(tmp)));
    yk = prob.data.*exp(1i*angle(tmp));
    zk = (yk+rho*tmp)/(1+rho);
    
    myfunc = @(x) myfunc2(x,prob.A,zk+lambda);
    pars.tol = 1e-6;
    pars.iteron = 0;
    if err < 0.5
        rho = 2.1;
        out = IIHT(prob.d1,sk,myfunc,pars);
    else % 3*sk
        out = IIHT(prob.d1,sk,myfunc,pars);
    end
    xk = out.x;
    lambdat = lambda+ zk-prob.A*xk; %zk-tmp;%lambda + (prob.A*xk-yk);
%     err = min(norm(prob.x0-xk),norm(prob.x0+xk))/norm(prob.x0);
    [err, x00] = computerelerror(xk,prob.x0);
    if opt.verbosity
        fprintf('error--%.3f || gap --%.3f\n',err,norm(lambdat-lambda));
    end
    lambda = lambdat;
    if err < delta
        break;
    end
end
% x = xk;
x = x00;
end