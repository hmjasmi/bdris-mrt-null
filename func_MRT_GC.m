function [Theta,R] = func_MRT_GC(W_norm,H_norm,NG)
% Compute the scattering matrix Theta given the normalized channels
% W_norm and H_norm, and the group size NG
% Inputs:   H_norm: normalized channel H
%           W_norm: normalized channel W
%           NG: group size (NG=N means fully connected)
%                          (NG=1 means single connected)
% Outputs:  Theta: scattering matrix symuni
%           R: scattering matrix relaxed 
N = size(W_norm,1);
G_matrix = W_norm*H_norm;
Theta = zeros(N,N);
R = zeros(N,N);
if(NG == 1) % Single connected
    R = diag(G_matrix');
    Theta = diag(R./abs(R));
    R = diag(R);
else % Group or fully connected
    G = N/NG;
    for g = 1:G
        G_matrix_g = G_matrix(NG*(g-1)+1:NG*g,NG*(g-1)+1:NG*g);
        R_g = G_matrix_g';
        R_g = sqrt(NG)*R_g./sqrt(trace(R_g*R_g'));
        Theta_g = alg_symuni(R_g);
        Theta(NG*(g-1)+1:NG*g,NG*(g-1)+1:NG*g) = Theta_g;
        R(NG*(g-1)+1:NG*g,NG*(g-1)+1:NG*g) = R_g;
    end

end