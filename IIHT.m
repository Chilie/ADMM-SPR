function out = IIHT(n,s, func, pars,varargin)
% A solver for sparsity constraints models:
%
%                    min f(x),  s.t. ||x||_0<=s,
%
%  where f: R^n->R and s<<n.
%
% Written by 16/01/2016, Shenglong Zhou
%
%
% Inputs:
%     n       : Dimension of the solution x, (required)
%     func    : function handle defines the function f(x) and its gradient             
%     pars:     Parameters are all OPTIONAL
%               pars.iteron --  Results will  be shown for each iteration if pars.iteron=1 (default)
%                               Results won't be shown for each iteration if pars.iteron=0 
%               pars.maxit  --  Maximum nonumber of iteration.  pars.maxit=5000 (default) 
%               pars.tol    --  Tolerance of stopping criteria. pars.maxit=1e-5 (default) 
%
% Outputs:
%     out.x:             The sparse solution x 
%     out.obj:           f(x)
%     out.time           CPU time
%     out.iter:          Number of iterations
%     out.error:         error ||g(x~=0)||
%     out.normg:         ||g||
%
%  This solver was created based on the algorithm proposed by  
%  Pan, L., Zhou, S., Xiu, N. & Qi, H.D. (2017). A convergent iterative hard 
%  thresholding for nonnegative sparsity optimization. Pacific Journal of 
%  Optimization, 13(2), 325-353.
%
%%%%%%%    Send your comments and suggestions to                     %%%%%%
%
%%%%%%%    shenglong.zhou@soton.ac.uk                                %%%%%%
% 
%%%%%%%    Warning: Accuracy may not be guaranteed!!!!!              %%%%%%
warning off;

if nargin<3; error('Imputs are not enough!\n'); end
if nargin<4; pars=[]; end
if isfield(pars,'iteron');iteron = pars.iteron; else; iteron = 1;        end
if isfield(pars,'maxit'); maxit  = pars.maxit;  else; maxit  = 1e4;      end
if isfield(pars,'tol');   tol    = pars.tol;    else; tol = 1e-6*sqrt(n);end  

t0     = tic;
sigma0 = 1e-4;
if nargin > 4
    x      = varargin{1};
xo     = varargin{1};
else
x      = zeros(n,1);
xo     = zeros(n,1);
end

% main body
if iteron 
fprintf(' Start to run the sover...\n'); 
fprintf('\n Iter    Error        f(x)       Time \n'); 
fprintf('---------------------------------------\n');
end
[f,g]  = func(x);
scale  = 0; 
if max(f,norm(g))/n>1; scale = 1; end
sl     = n*(scale==1)+(scale==0); 

for iter=1:maxit      
    x_old  = x;
    [f,g]  = func(x);
    if scale     
    f      = f/n;  
    g      = g/n;     
    end
    alpha  = sqrt(iter);
    
    % find a proper or the best alpha0   
    [mx,T] = maxk(x_old-alpha*g,s,'ComparisonMethod','abs');
    x      = xo; 
    x(T)   = mx;
    
    % Line search for setp size alpha
    fx_old = f;
    fx     = func(x);
    if scale     
    fx     = fx/n;    
    end
    
    for j  = 1:10
        if (fx < fx_old-.5*sigma0*sum((x-x_old).^2)); break; end
        alpha  = alpha/2;
        [mx,T] = maxk(x_old-alpha*g,s,'ComparisonMethod','abs');
        x      = xo; 
        x(T)   = mx;
        fx     = func(x);
        if scale; fx  = fx/n; end
    end

    % Stop criteria 
	residual = sl*norm(g(T))/max(1,norm(mx)); 
    if iteron
       fprintf('%4d    %5.2e    %5.2e   %5.2fsec\n',iter,residual,fx*sl,toc(t0)); 
    end
 
	if residual<tol || abs(fx-fx_old)<1e-10*(1+abs(fx))  
       break; 
    end  

end

if iteron
fprintf('---------------------------------------\n');
end

out.x    = x;
out.obj  = fx*sl;
out.iter = iter;
out.time = toc(t0);
out.error= residual;
out.normg= norm(g)*sl;
if  out.normg<1e-5 && iteron
fprintf(' A global optimal solution might be found\n');
fprintf(' because of ||g(x)||=%5.2e!\n',out.normg);  
end

end
