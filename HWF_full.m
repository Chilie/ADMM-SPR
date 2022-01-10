function [x,err,k] = HWF_full(prob,opt)
maxiter = opt.maxiter;
xk = opt.x0;
s0 = opt.s0;
gamma = opt.gamma;
[m,n] = size(prob.A);
delta = 1e-12;%opt.delta;
y = prob.data.^2; % squared data
% define the related loss and gradient
wf_loss = @(x) sum(abs(prob.A*x).^2-y)/(4*m);
wf_grad = @(x) prob.A'*((abs(prob.A*x).^2-y).*(prob.A*x))/m ;
% wf_grad = @(x) prob.A'*(prob.A*x - prob.data.*exp(1i*angle(prob.A*x)))/m ;
beta = 1e-8; % default parameter

% set the initialization
% estimate the signal norm
theta = sqrt(mean(y));
[~,ind0] = max(sum(y.*(abs(prob.A).^2)));
u_cur = sqrt(beta/2)*ones(n,1);
v_cur = sqrt(beta/2)*ones(n,1);
u_cur(ind0) = sqrt(sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
v_cur(ind0) = sqrt(-sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
x_cur0 = u_cur.^2 - v_cur.^2;%xk; %u_cur.^2 - v_cur.^2;
step = 0.1/mean(y)^2;
% step = 0.1;
for k = 1:maxiter
    % each iteration
    r = 2*step * wf_grad(x_cur0);
    if k>1 && err<0.5
        beta = beta /5;
        x_cur = x_cur0 - sqrt(abs(x_cur0).^2+ beta).*r;
    else
%     u_cur = u_cur.*(1-r);
%     v_cur = v_cur.*(1+r);
%     x_cur = u_cur.^2 - v_cur.^2;
    x_cur = x_cur0 - sqrt(abs(x_cur0).^2+ beta).*r;
    end
    
    % compute the metrics
    diffx = norm(x_cur-x_cur0);
    % xk = truncate(x_cur,s0);
%     err = min(norm(prob.x0-x_cur),norm(prob.x0+x_cur));
    [err_0, ~] = computerelerror(x_cur,prob.x0);
    [err, xk_tmp] = computerelerror(x_cur,prob.x0);
%     err_breg = min(breg_dist(prob.x0, x_cur,beta), breg_dist(-prob.x0, x_cur, beta));
    if opt.verbosity
%         fprintf('error--%.3f  error_breg--%.3f \n',err, err_breg);
        fprintf('iter: %d  error--%.3f \n',k,err);
    end
    if diffx < delta || err < 1e-3
        break;
    end
    x_cur0 = x_cur;   
end
% x = xk;
x = truncate(x_cur,s0);;
end