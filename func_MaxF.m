function [Theta,R] = func_MaxF(W_norm,H_norm,NG)
% Compute the scattering matrix Theta given the normalized channels
% W_norm and H_norm, and the group size NG
% Inputs:   H_norm: normalized channel H
%           W_norm: normalized channel W
%           NG: group size (NG=N means fully connected)
%                          (NG=1 means single connected)
% Outputs:  Theta: scattering matrix symuni
%           R: Relaxed scattering matrix

N = size(W_norm,1);
K = size(W_norm,2);
if(NG == 1) % Single-connected
    A = kr(W_norm.',H_norm);
elseif(NG == N)
    A = kron(W_norm.',H_norm);
else % Group-connected
    A = [];
    G = N/NG;
    for i = 1:G
        W_g = W_norm(1+NG*(i-1):i*NG,1:K);
        H_g = H_norm(1:K,1+NG*(i-1):i*NG);
        A_g = kron(W_g.',H_g);
        A = [A A_g];
    end
end

[v,D] = eig(A'*A);
D = diag(D);
[D_max, idx_max] = max(D);
v_dom = v(:,idx_max);

if(NG == 1) % Single connected
    Theta = diag(v_dom./abs(v_dom));
    R = diag(v_dom);
else % Group or fully connected
    G = N/NG;
    for g = 1:G
        v_g = v_dom(NG^2*(g-1)+1:NG^2*g,1);
        V_g = reshape(sqrt(NG)*v_g,[NG,NG]);
        Theta_g = alg_symuni(V_g);
        Theta(NG*(g-1)+1:NG*g,NG*(g-1)+1:NG*g) = Theta_g;
        R(NG*(g-1)+1:NG*g,NG*(g-1)+1:NG*g) = V_g;
    end
    if(G==1)
        Theta = Theta_g;
        R = V_g;
    end

end





