% This function is used to implement MRT at BS
% Inputs: - W: BS-BD-RIS Channel
%         - H: BD-RIS-Uk Channel
%         - Theta: BD-RIS Scattering matrix
%         - P_max: Max power budget at BS
% Output: Precoder matrix for MRT
% Note that normalization is done based on the whole matrix
% and not the column vectors of the matrix
function P = func_Prec_MRT(W,H,Theta,P_max)
E = H*Theta*W;
P = E';
P = sqrt(P_max).*P./sqrt(trace(P*P'));

