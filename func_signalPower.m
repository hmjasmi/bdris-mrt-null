function [P_R, P_I] = func_signalPower(x)
y = abs(x).^2;
P_R = diag(y);
K = length(P_R);
z = y.';
P_I = z(~eye(size(z)));
P_I = reshape(P_I,[K-1,K]).';


