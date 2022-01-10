function [x,err,k] = ADM(prob,opt)
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

% [m,n] = size(prob.A);
% gamma = opt.gamma;
%     % each iteration, extract a submatrix
%     ind = sort(randsample(m,floor(gamma*m)));
%     subA = prob.A(ind,:);
%     subdata = prob.data(ind);
%     suby = sign(subA*xk).*subdata;
    
%     % solve sparsity-constrained least-squares
%     myfunc = @(x) myfunc2(x,subA,suby);
%     pars.tol = 1e-6;
%     pars.iteron = 0;
%     pars.maxit = 500;
%     out = IIHT(prob.d1,sk,myfunc,pars);


myfunc = @(x) myfunc2(x,prob.A,zk+lambda);
    pars.tol = 1e-6;
    pars.iteron = 0;
    pars.maxit = 300;
    if err < 0.2
        rho = 0.1;%2.1;
        sk = s0;
        out = IIHT(prob.d1,sk,myfunc,pars);
        % sol = BIHT(A,K,y,mu);
    elseif err < 0.05 %0.05
        rho = 4.01;%2.1;
        sk = s0;
        out = IIHT(prob.d1,sk,myfunc,pars);
        % sol = BIHT(A,K,y,mu);    
    else % 3*sk
        out = IIHT(prob.d1,2*s0,myfunc,pars);
        % sol = BIHT(A,K,y,mu);
    end
    xk = out.x;
    if strcmp(prob.type, 'real')
        xk = real(xk);
    end

%     xk = BIHT(prob.A,sk,yk-lambda,0.01);
    tmp = prob.A*xk - lambda;
    if strcmp(prob.Atype, 'comp')
        yk = prob.data.*exp(1i*angle(tmp));
    else
        yk = real(prob.data.*exp(1i*angle(tmp)));
    end
    % yk = prob.data.*exp(1i*angle(tmp));
    
    zk = (yk+rho*tmp)/(1+rho);
    lambdat = zk-tmp;%lambda + (prob.A*xk-yk);
%     err = min(norm(prob.x0-xk),norm(prob.x0+xk))/norm(prob.x0);
    [err, xk_t] = computerelerror(xk,prob.x0);
    f_value = abs(prob.A*xk_t);
    err_lambda = norm(lambdat-lambda)/norm(lambda);
    if opt.verbosity
        fprintf('Iter--%d  error--%.3f || gap --%.3f || loss --%.3f\n',k, err,err_lambda, norm(f_value-prob.data)/norm(prob.data));
    end

    % the below is one criterion to stop the interation
    % err_lambda = norm(lambdat-lambdat)/norm(lambda);
    if err_lambda < 0.001
        break;
    end

    lambda = lambdat;
    if err < delta || norm(f_value-prob.data)/norm(prob.data) < 1e-3
        break;
    end
end
% x = xk;
x = xk_t;
end