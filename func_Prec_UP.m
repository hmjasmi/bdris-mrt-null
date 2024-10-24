% This function is used to implement uniform power allocation at BS
% Inputs: - K: Number of users
%         - P_max: Max power budget at BS
% Output: Precoder matrix for UP
function P = func_Prec_UP(P_max,K)
P = sqrt(P_max./K).*diag(ones(1,K));

