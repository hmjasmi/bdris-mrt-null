clc
clear
close all
num_iter = 1e2;
epsilon = 1e-6;
varepsilon = 1e-6;
M_range = 5;
K_range = M_range;
N = 64; 
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
        Theta_MRT_SC = func_MRT_GC(W_norm,H_norm,1);
        Theta_MRT_NG_2 = func_MRT_GC(W_norm,H_norm,2);
        Theta_MRT_NG_4 = func_MRT_GC(W_norm,H_norm,4);
        Theta_MRT_FC = func_MRT_GC(W_norm,H_norm,N);
        Theta_Null_SC = func_Nulling_GC(W_norm,H_norm,1,0,num_iter,epsilon,varepsilon);
        Theta_Null_NG_2 = func_Nulling_GC(W_norm,H_norm,2,0,num_iter,epsilon,varepsilon);
        Theta_Null_NG_4 = func_Nulling_GC(W_norm,H_norm,4,0,num_iter,epsilon,varepsilon);
        Theta_Null_FC = func_Nulling_GC(W_norm,H_norm,N,0,num_iter,epsilon,varepsilon);
        
        P_UP = func_Prec_UP(P_max,K);

        P_RM_SC = func_Prec_RM(W,H,Theta_MRT_SC,P_max,N0,3);
        P_RM_NG_2 = func_Prec_RM(W,H,Theta_MRT_NG_2,P_max,N0,3);
        P_RM_NG_4 = func_Prec_RM(W,H,Theta_MRT_NG_4,P_max,N0,3);
        P_RM_FC = func_Prec_RM(W,H,Theta_MRT_FC,P_max,N0,3);

        P_ZF_SC = func_Prec_ZF(W,H,Theta_MRT_SC,P_max);
        P_ZF_NG_2 = func_Prec_ZF(W,H,Theta_MRT_NG_2,P_max);
        P_ZF_NG_4 = func_Prec_ZF(W,H,Theta_MRT_NG_4,P_max);
        P_ZF_FC = func_Prec_ZF(W,H,Theta_MRT_FC,P_max);

        P_WF_SC = func_Prec_WF(W,H,Theta_Null_SC,P_max,N0,4);
        P_WF_NG_2 = func_Prec_WF(W,H,Theta_Null_NG_2,P_max,N0,4);
        P_WF_NG_4 = func_Prec_WF(W,H,Theta_Null_NG_4,P_max,N0,4);
        P_WF_FC = func_Prec_WF(W,H,Theta_Null_FC,P_max,N0,4);

        r = r + [...
                 func_compute_sR(P_UP,H*Theta_MRT_SC*W,N0);...
                 func_compute_sR(P_UP,H*Theta_MRT_NG_2*W,N0);...
                 func_compute_sR(P_UP,H*Theta_MRT_NG_4*W,N0);...
                 func_compute_sR(P_UP,H*Theta_MRT_FC*W,N0);...

                 func_compute_sR(P_RM_SC,H*Theta_MRT_SC*W,N0);...
                 func_compute_sR(P_RM_NG_2,H*Theta_MRT_NG_2*W,N0);...
                 func_compute_sR(P_RM_NG_4,H*Theta_MRT_NG_4*W,N0);...
                 func_compute_sR(P_RM_FC,H*Theta_MRT_FC*W,N0);...

                 func_compute_sR(P_ZF_SC,H*Theta_MRT_SC*W,N0);...
                 func_compute_sR(P_ZF_NG_2,H*Theta_MRT_NG_2*W,N0);...
                 func_compute_sR(P_ZF_NG_4,H*Theta_MRT_NG_4*W,N0);...
                 func_compute_sR(P_ZF_FC,H*Theta_MRT_FC*W,N0);...

                 func_compute_sR(P_UP,H*Theta_Null_SC*W,N0);...
                 func_compute_sR(P_UP,H*Theta_Null_NG_2*W,N0);...
                 func_compute_sR(P_UP,H*Theta_Null_NG_4*W,N0);...
                 func_compute_sR(P_UP,H*Theta_Null_FC*W,N0);...

                 func_compute_sR(P_WF_SC,H*Theta_Null_SC*W,N0);...
                 func_compute_sR(P_WF_NG_2,H*Theta_Null_NG_2*W,N0);...
                 func_compute_sR(P_WF_NG_4,H*Theta_Null_NG_4*W,N0);...
                 func_compute_sR(P_WF_FC,H*Theta_Null_FC*W,N0);...
                ];
        iter = iter + 1;
    end
    R(:,i_loop) = r./Iter;
    toc
end
