function sol = BIHT(A,K,y,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Author                 Zhang Cheng, Yang Hairong
%%Modified time          2010/07/23
%%function               Code for Backtracking Itetative Hard Thresholding
%% A  - Measurement matrix
%% K  - sparsity level
%% y  - measurement vector
%% mu - parameter as in IHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N]   = size(A);
s       = zeros(N,1);
Ps      = zeros(M,1);
Residual    = y;
MAXITER=2000;
Phi=A;

done=0;
iter=1;


while ~done
        theta               =   s+mu*Phi'*(y-Phi*s);
        [ssort sortind]     =   sort(abs(theta),'descend');
        theta(sortind(K+1:end)) =   0;
        
        Index_S=find(s~=0);
        Index_theta=find(theta~=0);
        activeset=union(Index_S,Index_theta);

        Phi_x = Phi(:,activeset);        
        beta=inv(Phi_x'*Phi_x)*Phi_x'*y;

        [bsort bsortind]     =   sort(abs(beta),'descend');
        beta(bsortind(K+1:end)) =   0;
        
        s(activeset) = beta;
     
     iter=iter+1;
     res=y-Phi_x*beta;
     err=norm(res);
     
     if iter >= MAXITER
         display('Stopping. Maximum number of iterations reached!')
         done = 1; 
     end
     
     if err < 1e-6
         display('Success!')
         done = 1; 
     end   
end

sol = zeros(N,1);
sol(activeset) = beta;