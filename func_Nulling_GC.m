function [Theta,obj_f] = func_Nulling_GC(W_norm,H_norm,NG,init,num_iter,epsilon,varepsilon)
% Compute the scattering matrix Theta given the normalized channels
% W_norm and H_norm, and the group size NG
% Inputs:   H_norm: normalized channel H
%           W_norm: normalized channel W
%           NG: group size (NG=N means fully connected)
%                          (NG=1 means single connected)
%           init: "1" random initialization
%                 "0" MRT initialization
%           num_iter: number of AO iterations
%           epsilon: AO convergence tol1
%           varepsilon: AO convergence tol2
% Outputs:  Theta: scattering matrix symuni

% clc
% clear
% close all
% N = 4;
% K = 2;
% M = 2;
% NG = 2;
% init = 0;
% num_iter = 1000;
% epsilon = 1e-8;
% varepsilon = 1e-8;
% W_norm = 1./sqrt(2).*(randn(N,M) + 1j*randn(N,M));
% H_norm = 1./sqrt(2).*(randn(K,N) + 1j*randn(K,N));
%%
N = size(W_norm,1);
K = size(W_norm,2);
A = kron(W_norm.',H_norm);
if(init == 1) 
    Theta_init = randn(N,N)+1j*randn(N,N);
elseif(init == 0)
    Theta_init = func_MRT_GC(W_norm,H_norm,NG);
end

if(NG == 1) % Single-connected
    A_tilde = kr(W_norm.',H_norm);
    for i = 1:size(A_tilde,1)
        [row,col] = ind2sub([K K],i);
        A_cell{row,col} = A_tilde(i,:).';
    end
    A_d = cell(1,K);
    [A_d{:}] = A_cell{~~eye(size(A_cell))};
    A_i = cell(1,K*(K-1));
    [A_i{:}] = A_cell{~eye(size(A_cell))};
    A_d = cat(2,A_d{:});
    A_i = cat(2,A_i{:});

    % AO
    Theta_tilde = diag(Theta_init);
    Theta_tilde = Theta_tilde./abs(Theta_tilde);
    obj_f(1) = norm(A_i.'*Theta_tilde)^2; 
    B = A_i.';
    A_pre = B'/( B* B')*B; 
    for i_iter = 1:num_iter
        Theta_tilde = Theta_tilde-A_pre*Theta_tilde;
        Theta_tilde(abs(Theta_tilde)<1e-10) = 1;
        Theta_tilde = Theta_tilde./abs(Theta_tilde); 
        obj_f(i_iter) = norm(A_i.'*Theta_tilde)^2; 
        if i_iter>1 && ((obj_f(i_iter) < epsilon) || (abs((obj_f(i_iter)-obj_f(i_iter-1))/obj_f(i_iter-1))<varepsilon))
            Theta = diag(Theta_tilde);
            break;
        end
    end
    Theta = diag(Theta_tilde);

elseif(NG == N) % Fully-connected
    A_tilde = A;
    for i = 1:size(A_tilde,1)
        [row,col] = ind2sub([K K],i);
        A_cell{row,col} = A_tilde(i,:).';
    end
    A_d = cell(1,K);
    [A_d{:}] = A_cell{~~eye(size(A_cell))};
    A_i = cell(1,K*(K-1));
    [A_i{:}] = A_cell{~eye(size(A_cell))};
    A_d = cat(2,A_d{:});
    A_i = cat(2,A_i{:});
    
    % AO
    Theta_tilde = Theta_init;
    [U,S,V] = svd(Theta_tilde);
    Theta_tilde = U*V'; 
    Theta_tilde = 0.5*(Theta_tilde + Theta_tilde.');
    obj_f(1) = norm(A_i.'*Theta_tilde(:))^2; 
    B = A_i.';
    A_pre = B'/( B* B')*B; 
    for i_iter = 1:num_iter
        Theta_tilde = Theta_tilde(:)-A_pre*Theta_tilde(:); 
        Theta_tilde(abs(Theta_tilde)<1e-10) = 1;
        Theta_tilde = reshape(Theta_tilde,[N N]);
        [U,S,V] = svd(Theta_tilde);
        Theta_tilde = U*V'; 
        Theta_tilde = 0.5*(Theta_tilde + Theta_tilde.');
        obj_f(i_iter) = norm(A_i.'*Theta_tilde(:))^2; 
        if i_iter>1 && ((obj_f(i_iter) < epsilon) || (abs((obj_f(i_iter)-obj_f(i_iter-1))/obj_f(i_iter-1))<varepsilon))
            Theta = Theta_tilde;
            break;
        end
    end
    Theta = Theta_tilde;
else % Group-connected
    A_tilde = [];
    G = N/NG;
    idx = [];
    for i = 1:G
        W_g = W_norm(1+NG*(i-1):i*NG,1:K);
        H_g = H_norm(1:K,1+NG*(i-1):i*NG);
        A_tilde_g = kron(W_g.',H_g);
        A_tilde = [A_tilde A_tilde_g];
        idx = blkdiag(idx,ones(NG,NG)); 
    end
    for i = 1:size(A_tilde,1)
        [row,col] = ind2sub([K K],i);
        A_cell{row,col} = A_tilde(i,:).';
    end
    A_d = cell(1,K);
    [A_d{:}] = A_cell{~~eye(size(A_cell))};
    A_i = cell(1,K*(K-1));
    [A_i{:}] = A_cell{~eye(size(A_cell))};
    A_d = cat(2,A_d{:});
    A_i = cat(2,A_i{:});

    % AO

    Theta_tilde = Theta_init(~~idx);
    Theta_tilde2 = [];
    for i = 1:G
        Theta_tilde_g = Theta_tilde(1+NG^2*(i-1):i*NG^2);
        Theta_tilde_g = reshape(Theta_tilde_g,[NG NG]);
        [U,S,V] = svd(Theta_tilde_g);
        Theta_tilde_g = U*V';
        Theta_tilde_g = 0.5*(Theta_tilde_g + Theta_tilde_g.');
        Theta_tilde2 = blkdiag(Theta_tilde2,Theta_tilde_g);
    end
    Theta_tilde = Theta_tilde2(~~idx);
    obj_f(1) = norm(A_i.'*Theta_tilde(:))^2;
    B = A_i.';
    A_pre = B'/( B* B')*B; 
    for i_iter = 1:num_iter
        Theta_tilde = Theta_tilde(:)-A_pre*Theta_tilde(:);
        Theta_tilde(abs(Theta_tilde)<1e-10) = 1;
        Theta_tilde2 = [];
        for i = 1:G
            Theta_tilde_g = Theta_tilde(1+NG^2*(i-1):i*NG^2);
            Theta_tilde_g = reshape(Theta_tilde_g,[NG NG]);
            [U,S,V] = svd(Theta_tilde_g);
            Theta_tilde_g = U*V';
            Theta_tilde_g = 0.5*(Theta_tilde_g + Theta_tilde_g.');
            Theta_tilde2 = blkdiag(Theta_tilde2,Theta_tilde_g);
        end
        Theta_tilde = Theta_tilde2(~~idx);
        obj_f(i_iter) = norm(A_i.'*Theta_tilde(:))^2;
        if i_iter>1 && ((obj_f(i_iter) < epsilon) || (abs((obj_f(i_iter)-obj_f(i_iter-1))/obj_f(i_iter-1))<varepsilon))
            Theta = Theta_tilde2;
            break;
        end
    end
    Theta = Theta_tilde2;
end
