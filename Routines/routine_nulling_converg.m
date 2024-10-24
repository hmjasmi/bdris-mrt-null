% This routine is used to generate the results for 
% the Convergence Figure (Fig. 2)
% 100 realizations are considered
% 1E4 number of iterations for the algorithm
% epsilon = 1e-10
% varepsilon = 1e-6
% N = 144
% K = 8
% No pathloss is considered
clc
clear
close all
N = 144;
K = 8;
M = K;
PL = 10.^(-[zeros(1,K)].'./10);
c = 3e8;        % Speed of light
f = 2.4e9;      % Carrier frequency
lambda = c/f;   % Wavelength
d_BS_BDRIS = 50;
d_BDRIS_Uk = 10*ones(K,1);
d0 = 1;
rho = 2.7;
C_0 = (lambda/4/pi/d0)^2;
PL_BS_BDRIS = C_0*(d_BS_BDRIS/d0).^-rho;
PL_BDRIS_Uk = C_0*(d_BDRIS_Uk/d0).^-rho;
PL_dB = 10.*log10([PL_BS_BDRIS;PL_BDRIS_Uk]);

Iter = 100;
num_iter = 1e4;
epsilon = 1e-10;
varepsilon = 1e-6;
iter = 1;
while(iter<=Iter)
    disp([iter])
    W = sqrt(PL_BS_BDRIS)./sqrt(2).*(randn(N,M) + 1j*randn(N,M));
    H = sqrt(PL_BDRIS_Uk)./sqrt(2).*(randn(K,N) + 1j*randn(K,N));
    W_norm = W./sqrt(PL_BS_BDRIS);
    H_norm = H./sqrt(PL_BDRIS_Uk);
    [~,f_SC{iter}] = func_Nulling_GC(W_norm,H_norm,1,1,num_iter,epsilon,varepsilon);
    [~,f_NG_2{iter}] = func_Nulling_GC(W_norm,H_norm,2,1,num_iter,epsilon,varepsilon);
    [~,f_NG_4{iter}] = func_Nulling_GC(W_norm,H_norm,4,1,num_iter,epsilon,varepsilon);
    [~,f_FC{iter}] = func_Nulling_GC(W_norm,H_norm,N,1,num_iter,epsilon,varepsilon);
    iter = iter + 1;
end
