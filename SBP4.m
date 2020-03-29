% The operators were originally constructed in the following paper
% @article{Mattsson2006,
% Author={K.~Mattsson and J. Nordstr\"{o}m},
% Title={High order finite difference methods for wave propagation in discontinuous media},
% Journal={J. Comput. Phys.},
% Volume={220},
% Year={2006},
% Pages={249--269},
% }
function [H,D1,D2,S] = SBP4(m,h)
H=diag(ones(m,1),0);
H(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
H(m-3:m,m-3:m)=fliplr(flipud(diag([17/48 59/48 43/48 49/48])));   
  
D1=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)- ...
    8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));

D1(1:4,1:6)=[-24/17,59/34,-4/17,-3/34,0,0; -1/2,0,1/2,0,0,0; 4/43,-59/86,0,59/86,-4/43,0; 3/98,0,-59/98,0,32/49,-4/49];
D1(m-3:m,m-5:m)=flipud( fliplr(-D1(1:4,1:6)));
D1=D1/h;


D2=(-1/12*diag(ones(m-2,1),2)+16/12*diag(ones(m-1,1),1)+16/12* ...
    diag(ones(m-1,1),-1)-1/12*diag(ones(m-2,1),-2)-30/12*diag(ones(m,1),0));

m55=-5/2;
s15=0;

D2(1:5,1:7)=[154/17+48/17*m55-48/17*s15, -565/17-192/17*m55+192/17*s15, ...
	     788/17+288/17*m55-288/17*s15, -497/17-192/17*m55+192/17*s15, ...
	     120/17+48/17*m55-48/17*s15, 0, 0; -421/59-192/59*m55, ...
	     1802/59+768/59*m55, -2821/59-1152/59*m55, 1920/59+768/59*m55, ...
	     -480/59-192/59*m55, 0, 0; 716/43+288/43*m55, -2821/43-1152/43*m55, 4210/43+1728/43*m55, -2821/43-1152/43*m55, 716/43+288/43*m55, 0, 0; -481/49-192/49*m55, 1920/49+768/49*m55, -403/7-1152/49*m55, 1802/49+768/49*m55, -416/49-192/49*m55, -4/49, 0; 5/2+m55, -10-4*m55, 179/12+6*m55, -26/3-4*m55, m55, 4/3, -1/12];


D2(m-4:m,m-6:m)=flipud( fliplr(D2(1:5,1:7) ) );


D2=D2/(h^2);

BS=zeros(m,m);
BS(1,1:5)=-[s15-11/6,-4*s15+3,6*s15-3/2,-4*s15+1/3,s15];
BS(m,m-4:m)=fliplr(-[s15-11/6,-4*s15+3,6*s15-3/2,-4*s15+1/3,s15]);
BS=BS/h;


H=h*H;
HI=inv(H);

M=BS-H*D2;
Q=H*D1;

S = BS;
S(1,:)=S(1,:)*(-1);
for i = 2:size(S,1)-1
    S(i,i) = 1/h;
end
S = sparse(S);
H = sparse(H);
D2 = sparse(D2);
