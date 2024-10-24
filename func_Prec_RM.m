% This function is used to implement Rate Maximization Precoding at BS
% Inputs: - W: BS-BD-RIS Channel
%         - H: BD-RIS-Uk Channel
%         - Theta: BD-RIS Scattering matrix
%         - P_max: Max power budget at BS
%         - init: "1" for UP initialization
%                 "2" for Random initialization
%                 "3" for WF initialization
%         - N0: Noise PSD
% Output: Precoder matrix for RM
% Note that WF initialization is used
function [P,sR] = func_Prec_RM(W,H,Theta,P_max,N0,init)
E = H*Theta*W;
K = size(E,1);

% Fmincon proceedure
if(init==1) % UP
    x0 = diag(func_Prec_UP(P_max,K));
elseif(init==2) % Random
    x0 = randn(K,1);
    x0 = sign(x0).*x0;
elseif(init==3) % WF
    x0 = diag(func_Prec_WF(W,H,Theta,P_max,N0,4));
end

Aa = -(eye(K,K));
b = zeros(K,1);
Aeq = ones(1,K);
Beq = P_max;
lb = zeros(K,1);
ub = P_max*ones(K,1);
nonlcon = [];
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'UseParallel',false,'FunctionTolerance',1e-6,'StepTolerance',1e-6, ...
    'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8);

[P,sR] = fmincon(@(P_n)-func_sR(P_n,E,N0),x0,Aa,b,Aeq,Beq,lb,ub,nonlcon,options);
sR = abs(sR);

P = diag(sqrt(P));
