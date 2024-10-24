function sR = func_sR(P_n,E,N0)
sigma = sqrt(N0);
[P_R, P_I] = func_signalPower(abs(E)*diag(sqrt(P_n)));
SINR = P_R./(sum(P_I,2) + sigma.^2);
R = log2(1 + SINR);
sR = sum(R,1);
