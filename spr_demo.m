% define the problem
addpath(genpath('.'));
Alg =  {'ADM','StormSpr','COSMAP','SparTAF', 'HWF_full'};
TestSp = [30:10:80]; % tested sparsity levels
OS = [2:0.25:3];
relax_factor = 2;

Nprob = 10;

nalg = length(Alg);
nsp = length(TestSp);
nos = length(OS);
Recov = zeros(Nprob,nalg,nos,nsp);
Recovredius = zeros(Nprob,nalg,nos,nsp);
Itertrace = zeros(Nprob,nalg,nos,nsp);
TTimer = zeros(Nprob,nalg,nos,nsp);

prob.ndim = 1;
prob.type = 'real';
prob.Atype = 'real';
prob.d1 = 400;
prob.d2 = 1;
n = prob.d1;
for isp = 1:nsp
    % generate signal
    sk = TestSp(isp);
    SK = floor(TestSp(isp)*(log(prob.d1))); % +log(100)
    for ios = 1:nos
        prob.os = OS(ios);
        m = round(prob.os*SK);
        
        % generate the problem
        % prob.x0 = randn(n,1);
        for iprob = 1:Nprob
            fprintf('\n Solve the #%d problem with os %.2e and Sparisty %d\n',iprob,ios,isp);
            prob.x0 = zeros(n,1);
            omega = randsample(n,sk);
            prob.x0(omega) = randn(sk,1);
%             prob.x0 = prob.x0/norm(prob.x0);
            prob.A = randn(m,n);
            
            prob.data = abs(prob.A*prob.x0);
            
            SNR = 0;%30
            if SNR
                prob.data = max(awgn(prob.data,SNR,'measured'),realmin);
            end
            % solve this problem
            
            % parameter
            opt.maxiter = 10000;
%             x0=zeros(n,1);
%             sk2 = 2*sk;
%             omega = randsample(n,sk2);
%             x0(omega) = 10*randn(sk2,1);

            opt.rho = 0.01;%0.01;
            opt.sk = relax_factor*sk;%2*sk
            opt.verbosity = 1;
            opt.s0 = sk;
            opt.delta = 0.01;
            
            y = prob.data.^2; % squared data
            beta = 1e-8; % default parameter
            
            % set the initialization
            % estimate the signal norm
%             theta = sqrt(mean(y));
%             [~,ind0] = max(sum(y.*(prob.A.^2)));
%             x0 = sqrt(beta/2)*randn(n,1);
%             u_cur = sqrt(beta/2)*randn(n,1);
%             v_cur = sqrt(beta/2)*randn(n,1);
%             u_cur(ind0) = sqrt(sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
%             v_cur(ind0) = sqrt(-sqrt(theta^2/12) + sqrt(theta^2/12 + beta^2/4));
%             x0(ind0) = u_cur(ind0).^2 - v_cur(ind0).^2;
            % generate the initilaization
%             x0 = randn(n,1);
% %             omega = randsample(n,opt.sk);
% %             x0(omega) = randn(opt.sk,1)/100;
%             [~,idnx] = max(abs(prob.x0));
%             x0(idnx) = 3*randn;
                        x0 = randn(n,1);
%                         [x0] = x_initial(prob,opt);
            opt.x0 = x0;
            
            opt.gamma = 0.6;%min(0.6,opt.sk*log10(n/0.001)/m);
            opt.verbosity = 0;
            opt.sk = relax_factor*sk;
            opt.s0 = sk;
            opt.delta = 0.01;
            opt.verbosity = 1;
            for ialg = 1:nalg
                tic;
                [x,err,k] = solve_spr(prob,opt,Alg{ialg});
                TTimer(iprob,1,ios,isp) = toc;
                Recovredius(iprob,1,ios,isp) = err;
                if err< 0.01
                    Recov(iprob,ialg,ios,isp) = 1;
                end
                Itertrace(iprob,ialg,ios,isp) = k;
            end
%             save('DataT1_random_init.mat','Recov','Recovredius','TTimer','Itertrace');
        end
%         save('DataT1_random_init.mat','Recov','Recovredius','TTimer','Itertrace');
    end
%     save('DataT1_random_init.mat','Recov','Recovredius','TTimer','Itertrace');
end



% function [x] = truncate(x,s)
% y = sort(abs(x),'descend');
% x = zeros(size(x));
% x(1:s) = y(1:s);
% end