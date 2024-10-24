% This is a code taken from the following paper:
% https://github.com/YijieLinaMao/BD-RIS-low-complexity/tree/main

% T. Fang and Y. Mao, "A Low-Complexity Beamforming Design for Beyond-Diagonal RIS Aided Multi-%User Networks," in IEEE Communications Letters, vol. 28, no. 1, pp. 203-207, Jan. 2024, doi: %10.1109/LCOMM.2023.3333411.



function out=alg_symuni(A)
% symmetric unitary projection
% input a square matrix
% outout a symmetric unitary matrix
[U,S,V]=svd(A+A.');
R=rank(S);
N=size(A,1);
U(:,R+1:N)=conj(V(:,R+1:N));
out=U*V';

end
