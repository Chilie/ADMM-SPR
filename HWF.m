function [x,err,k] = HWF(prob,opt)
maxiter = opt.maxiter;
xk = opt.x0;
sk = opt.sk;
s0 = opt.s0;
gamma = opt.gamma;
[m,n] = size(prob.A);
delta = 1e-12;%opt.delta;
delta_e = opt.delta;
y = prob.data.^2; % squared data
% define the related loss and gradient
wf_loss = @(x) sum(abs(prob.A*x).^2-y)/(4*m);
wf_grad = @(x) prob.A'*((abs(prob.A*x).^2-y).*(prob.A*x))/m ;
beta = 1e-8; % default parameter

% set the initialization
% estimate the signal norm
theta = sqrt(mean(y));
[~,ind0] = max(sum(y.*(abs(prob.A).^2)));
u_cur = sqrt(beta/2)*ones(n,1);
v_cur = sqrt(beta/2)*ones(n,1);
u_cur(ind0) = sqrt(sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
v_cur(ind0) = sqrt(-sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
x_cur0 = u_cur.^2 - v_cur.^2;
step = 0.1/mean(y).^(3/2);
for k = 1:maxiter
    % each iteration
    r = 2*step * wf_grad(x_cur0);
    if k > 1 && err < 0.5
%         step = step/5;
        x_cur = truncate(x_cur0 - r,sk);
    else
        u_cur = u_cur.*(1-r);
        v_cur = v_cur.*(1+r);
        x_cur = u_cur.^2 - v_cur.^2;
    end
    
    if k == 1000
       step = step /10;
    end
    % compute the metrics
%     x_tmp = truncate(x_cur,sk);
    diffx = norm(x_cur-x_cur0);
    xk = truncate(x_cur,s0);
%     err = min(norm(prob.x0-xk),norm(prob.x0+xk))/norm(prob.x0);
    [err, ~] = computerelerror(xk,prob.x0);
%     err_breg = min(breg_dist(prob.x0, xk,beta), breg_dist(-prob.x0, xk, beta));
%     err_breg = breg_dist(prob.x0, xk_t,beta);
    if opt.verbosity
%         fprintf('iter %d error--%.3f  error_breg--%.3f norm_r -- %.3e\n',k, err, err_breg,norm(r));
        fprintf('iter %d error--%.3f  norm_r -- %.3e\n',k, err,norm(r));
    end
    if diffx < delta || err < delta_e %diffx < delta || 
        break;
    end
    if err < 0.5
        x_cur0 = truncate(x_cur,sk);
    else
        x_cur0 = x_cur;   
    end
end
x = xk;
end