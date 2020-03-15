function [H,D2,S] = SBP_variable_2(m,h,c)

% Norm
Hv = ones(m,1);
Hv(1) = 1/2;
Hv(m:m) = 1/2;
Hv = h*Hv;
H = spdiags(Hv, 0, m, m);
HI = spdiags(1./Hv, 0, m, m);

% Boundary operators
e_l = sparse(m,1);
e_l(1) = 1;
e_r = rot90(e_l, 2);

d1_l = sparse(m,1);
d1_l(1:3) = 1/h*[-3/2 2 -1/2];
d1_r = -rot90(d1_l, 2);

S = sparse(m,m);
S(1,1:3) = [3/2 -2 1/2]*(-1);
S(end,end-2:end) = [1/2 -2 3/2];
S = S/h;

M=sparse(m,m);

scheme_width = 3;
scheme_radius = (scheme_width-1)/2;
r = (1+scheme_radius):(m-scheme_radius);

Mm1 = -c(r-1)/2 - c(r)/2;
M0  =  c(r-1)/2 + c(r)   + c(r+1)/2;
Mp1 =            -c(r)/2 - c(r+1)/2;

M(r,:) = spdiags([Mm1 M0 Mp1],0:2*scheme_radius,length(r),m);


M(1:2,1:2)=[c(1)/2 + c(2)/2 -c(1)/2 - c(2)/2; -c(1)/2 - c(2)/2 c(1)/2 + c(2) + c(3)/2;];
M(m-1:m,m-1:m)=[c(m-2)/2 + c(m-1) + c(m)/2 -c(m-1)/2 - c(m)/2; -c(m-1)/2 - c(m)/2 c(m-1)/2 + c(m)/2;];
M=M/h;

D2=HI*(-M-c(1)*e_l*d1_l'+c(m)*e_r*d1_r');

end