% The operators were originally constructed in the following paper
% @article{Mattsson2006,
% Author={K.~Mattsson and J. Nordstr\"{o}m},
% Title={High order finite difference methods for wave propagation in discontinuous media},
% Journal={J. Comput. Phys.},
% Volume={220},
% Year={2006},
% Pages={249--269},
% }
function [H,D1,D2,S] = SBP2(N,h)
H = speye(N,N);
H(1,1) = 1/2;
H(end,end) = 1/2;
H = H*h;
e = ones(N,1);
D2 = spdiags([e -2*e e], -1:1, N, N);
D2(1,1:3) = [1 -2 1]; 
D2(end,end-2:end) = [1 -2 1]; 


D1 = spdiags([-0.5*e -0*e 0.5*e], -1:1, N, N);
D1(N,N-1:N) = [-1,1];
D1(1,1:2) = [-1,1];

D2 = D2/h/h; D1 = D1/h;
S = speye(N,N);
S(1,1:3) = [3/2 -2 1/2]*(-1);
S(end,end-2:end) = [1/2 -2 3/2];
S = S/h;
