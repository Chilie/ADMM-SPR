% solve problom
function [x,err,k] = solve_spr(prob,opt, method)
%SOLVE_SPR solve sparse phase retrieval problem using different methods
if strcmp(method, 'ADM')
    fprintf('\nRunning ADMM . . .\n');
    [x,err,k] = ADM(prob,opt); % ADM_split
elseif strcmp(method, 'ADM2')
    fprintf('\nRunning ADMM . . .\n');
    [x,err,k] = ADM2(prob,opt);
elseif strcmp(method, 'StormSpr')
    fprintf('\nRunning StormSpr . . .\n');
    [x,err,k] = StormSpr(prob,opt);
elseif strcmp(method, 'COSMAP')
    fprintf('\nRunning COSMAP . . .\n');
    [x,err,k] = myCoPRAM(prob,opt);  
elseif strcmp(method, 'SparTAF')
    fprintf('\nRunning SparTAF . . .\n');
    [x,err,k] = mySparTAF(prob,opt);   
elseif strcmp(method, 'HWF')
    fprintf('\nRunning HWF . . .\n');
    [x,err,k] = HWF(prob,opt);
elseif strcmp(method, 'HWF_r')
    fprintf('\nRunning HWF . . .\n');
    [x,err,k] = HWF_r(prob,opt);
elseif strcmp(method, 'HWF_full')
    fprintf('\nRunning HWF_full . . .\n');
    [x,err,k] = HWF_full(prob,opt);
elseif strcmp(method, 'ADM_HWF')
    fprintf('\nRunning ADM_HWF . . .\n');
    [x,err,k] = ADM_HWF(prob,opt);
elseif strcmp(method, 'ADM_Linear')
    fprintf('\nRunning ADM_Linear . . .\n');
    [x,err,k] = ADM_Linear(prob,opt);
else
    error('No related method is found!');
end
