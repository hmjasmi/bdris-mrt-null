% This function is used to implement ZF at BS
% Inputs: - W: BS-BD-RIS Channel
%         - H: BD-RIS-Uk Channel
%         - Theta: BD-RIS Scattering matrix
%         - P_max: Max power budget at BS
% Output: Precoder matrix for ZF
% Note that normalization is done based on the whole matrix
% and not the column vectors of the matrix
function P = func_Prec_ZF(W,H,Theta,P_max)
E = H*Theta*W;
P = inv(E'*E)*E';
P = sqrt(P_max).*P./sqrt(trace(P*P'));

