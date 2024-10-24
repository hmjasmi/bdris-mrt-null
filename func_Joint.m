function [Theta,P,sR] = func_Joint(W_norm,H_norm,NG,P_max,N0)
% Compute the scattering matrix Theta and the precoding matrix
% given the normalized channels
% W_norm and H_norm, and the group size NG
% Inputs:   H_norm: normalized channel H
%           W_norm: normalized channel W
%           NG: group size (NG=N means fully connected)
%                          (NG=1 means single connected)
%           P_max: Max power budget at BS
%           N0: Noise PSD
% Outputs:  Theta: scattering matrix
%           P: precoding matrix
%           R: scattering matrix relaxed 

N = size(W_norm,1);
K = size(W_norm,2);
global sigma
global SIM
global N_max
global eta
sigma = sqrt(N0);
SIM = 20;       % No. Iterations in Alg. 1 
N_max = 50;     % Maximum No. Iterations to find optimal lambda using bi-section search 
eta = 1e-3;     % Convergence tol

H.Hr = H_norm';
H.G = W_norm;

P = P_max;

if(NG == 1)
    [obj_r_f,sR] = FP_single_connected_tor(H.Hr,H.G,P);
    P = obj_r_f.w;
    Theta = obj_r_f.phy;
elseif(NG == N)
    [obj_r_f,sR] = FP_fully_connected_tor(H.Hr,H.G,P);
    P = obj_r_f.w;
    Theta = obj_r_f.phy;
else
    S = N/NG;
    [obj_r_f,sR] = FP_group_connected_tor(H.Hr,H.G,P,S);
    P = obj_r_f.w;
    Theta = obj_r_f.phy;

end

