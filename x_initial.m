%signal marginals
function [x] = x_initial(prob,opt)
A = prob.A;
y_abs = prob.data;
y_abs2 = prob.data.^2;
[m,n] = size(prob.A);
s = opt.s0;

MShat = zeros(s);

phi_sq = sum(y_abs2)/m;
phi = sqrt(phi_sq); %signal power
Marg = ((y_abs2)'*(abs(A).^2))'/m; % n x 1
[Mg MgS] = sort(Marg,'descend');
S0 = MgS(1:s); %pick top s-marginals
Shat = sort(S0); %store indices in sorted order
%supp(Shat) = 1; figure; plot(supp); %support indicator
AShat = A(:,Shat); % m x s %sensing sub-matrix

%% Truncated measurements
card_Marg = ceil(m/6);
%large measurements - amplitude flow
for i=1:m
    M_eval(i) = y_abs(i)/norm(AShat(i,:));
end
[Mm MmS] = sort(M_eval,'descend');
Io = MmS(1:card_Marg); %indices between 1 to m

%% Initialize x
%compute top singular vector according to thresholded sensing vectors and large measurements
for i = 1:card_Marg
    ii = Io(i);
    MShat = MShat + (y_abs2(ii))*AShat(ii,:)'*AShat(ii,:); % (s x s)
end

svd_opt = 'svd'; %more accurate, but slower for larger dimensions
svd_opt = 'power'; %approximate, faster

switch svd_opt
    case 'svd'
        [u,sigma,v] = svd(MShat);
        v1 = u(:,1); %top singular vector of MShat, normalized - s x 1
    case 'power'
        v1 = svd_power(MShat);
end

v = zeros(n,1);
v(Shat,1) = v1;
x_init = phi*v; %ensures that the energy/norm of the initial estimate is close to actual
x = x_init;
end