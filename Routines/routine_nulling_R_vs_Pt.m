clc
clear
close all
num_iter = 1e2;
epsilon = 1e-6;
varepsilon = 1e-6;
M_range = 5;
K_range = M_range;
N = 32; % {32, 64}
warning('off','all')
Iter = 1e2;
P_tx_dBm = 0:5:20;
P_tx = 10.^(P_tx_dBm./10)*1e-3;
N0_dBm = -80;
N0 = 10.^(N0_dBm./10)*1e-3;
data_length = 1;
SNR_dB = P_tx_dBm - N0_dBm;
SNR = 10.^(SNR_dB./10);
c = 3e8;        % Speed of light
f = 2.4e9;      % Carrier frequency
lambda = c/f;   % Wavelength
d_BS_BDRIS = 50;
d_BDRIS_Uk = 2.5;
d0 = 1;
rho = 2.2;
C_0 = 10.^(-30/10);
PL_BS_BDRIS = C_0*(d_BS_BDRIS/d0).^-rho;
PL_BDRIS_Uk = C_0*(d_BDRIS_Uk/d0).^-rho;
W_all = cell(1,Iter);
H_all = cell(1,Iter);
for i_iter = 1:Iter
    W_all{i_iter} = sqrt(PL_BS_BDRIS)./sqrt(2).*(randn(N,M_range(end)) + 1j*randn(N,M_range(end)));
    H_all{i_iter} = sqrt(PL_BDRIS_Uk)./sqrt(2).*(randn(K_range(end),N) + 1j*randn(K_range(end),N));
end

parfor i_loop = 1:length(P_tx_dBm)
    i_loop
    tic
    warning('off','all')
    iter = 0;
    r = 0;
    iter = 0;
    K = K_range;
    M = M_range;
    P_max = P_tx(i_loop);
    while(iter<Iter)
        W = W_all{iter+1}(1:N,1:M);
        H = H_all{iter+1}(1:K,1:N);
        W_norm = W./sqrt(PL_BS_BDRIS);
        H_norm = H./sqrt(PL_BDRIS_Uk);
        Theta_SC = func_Nulling_GC(W_norm,H_norm,1,0,num_iter,epsilon,varepsilon);
        Theta_NG_2 = func_Nulling_GC(W_norm,H_norm,2,0,num_iter,epsilon,varepsilon);
        Theta_NG_4 = func_Nulling_GC(W_norm,H_norm,4,0,num_iter,epsilon,varepsilon);
        Theta_FC = func_Nulling_GC(W_norm,H_norm,N,0,num_iter,epsilon,varepsilon);
        P = func_Prec_UP(P_max,K);
        r = r + [...
                 func_compute_sR(P,H*Theta_SC*W,N0);...
                 func_compute_sR(P,H*Theta_NG_2*W,N0);...
                 func_compute_sR(P,H*Theta_NG_4*W,N0);...
                 func_compute_sR(P,H*Theta_FC*W,N0);...
                ];
        iter = iter + 1;
    end
    R(:,i_loop) = r./Iter;
    toc
end
