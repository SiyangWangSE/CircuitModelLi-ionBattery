r_distance = (0:h_s_neg:r_s_neg)';
ny = length(r_distance);
% SBP operator D2 approximates second derivative
% H: discrete norm 
% S: approximates first derivative on the boundary
if order == 4
    [H,D2,S] = SBP_variable_4(ny,h_s_neg,r_distance.^2);
elseif order == 2
    [H,D2,S] = SBP_variable_2(ny,h_s_neg,r_distance.^2);
else
    error('order can only be 2 or 4')
end
% build the SBP operator with the SAT boundary treatment
B = sparse(ny,ny); B(1,1) = 0; B(ny,ny) = r_s_neg.^2;
invHB_cs_neg = inv(H)*B;
Laplacian_cs_neg = D2-invHB_cs_neg*S;
g_neg = zeros(ny,1);
M = diag(r_distance.^2);
options_c_s_neg = odeset('Mass',M); % used in ode15s

%%% construct the finite difference operator for function dcsdt = calculate_cs_pos(~,cs)
r_distance = (0:h_s_pos:r_s_pos)';
ny = length(r_distance);
% SBP operator D2 approximates second derivative
% H: discrete norm 
% S: approximates first derivative on the boundary
if order == 4
    [H,D2,S] = SBP_variable_4(ny,h_s_pos,r_distance.^2);
elseif order == 2
    [H,D2,S] = SBP_variable_2(ny,h_s_pos,r_distance.^2);
else
    error('order can only be 2 or 4')
end
% build the SBP operator with the SAT boundary treatment
B = sparse(ny,ny); B(1,1) = 0; B(ny,ny) = r_s_pos.^2;
Laplacian_cs_pos = D2-inv(H)*B*S;
invHB_cs_pos = inv(H)*B;
g_pos = zeros(ny,1);
M = diag(r_distance.^2);
options_c_s_pos = odeset('Mass',M); % used in ode15s

%%% constrct the finite difference operator for function dcedt = calculate_ce(~,ce)
% variable coefficient
e1 = eps_l(1).^Bruggeman_neg;
e2 = eps_l(n_mesh_neg+3).^Bruggeman_sep;
e3 = eps_l(end).^Bruggeman_pos;

% SBP operator
if order == 4
    [H1,~,D21,S1] = SBP4(n_mesh_neg+1,h_neg);
    [H2,~,D22,S2] = SBP4(n_mesh_sep+1,h_sep);
    [H3,~,D23,S3] = SBP4(n_mesh_pos+1,h_pos);
else
    [H1,~,D21,S1] = SBP2(n_mesh_neg+1,h_neg);
    [H2,~,D22,S2] = SBP2(n_mesh_sep+1,h_sep);
    [H3,~,D23,S3] = SBP2(n_mesh_pos+1,h_pos);
end
% difference operator
D2 = [e1*D21, sparse(n_mesh_neg+1, n_mesh_sep+n_mesh_pos+2)
    sparse(n_mesh_sep+1,n_mesh_neg+1), e2*D22, sparse(n_mesh_sep+1,n_mesh_pos+1)
    sparse(n_mesh_pos+1,n_mesh_neg+n_mesh_sep+2), e3*D23];

% pick-up matrices
E11 = sparse(n_mesh_neg+1,n_mesh_neg+1); E11(1,1) = 1;
E1n = sparse(n_mesh_neg+1,n_mesh_neg+1); E1n(end,end) = 1;
E21 = sparse(n_mesh_sep+1,n_mesh_sep+1); E21(1,1) = 1;
E2n = sparse(n_mesh_sep+1,n_mesh_sep+1); E2n(end,end) = 1;
E12m = sparse(n_mesh_neg+1,n_mesh_sep+1); E12m(n_mesh_neg+1,1) = 1;
E21m = E12m';
E31 = sparse(n_mesh_pos+1,n_mesh_pos+1); E31(1,1) = 1;
E3n = sparse(n_mesh_pos+1,n_mesh_pos+1); E3n(end,end) = 1;
E23m = sparse(n_mesh_sep+1,n_mesh_pos+1); E23m(n_mesh_sep+1,1) = 1;
E32m = E23m';

% boundaries
SATB1 = +e1*H1\E11*S1;
SATB3 = -e3*H3\E3n*S3;

% interfaces penalty parameters
a1 = -0.5; a2 = -0.5; b1 = 0.5; b2 = 0.5;

% interface condition 1 & 2
SAT11 = a1*e1*inv(H1)*S1'*E1n + 1*a2*e1*inv(H1)*E1n*S1;
SAT12 = -a1*e1*inv(H1)*S1'*E12m - 1*a2*e2*inv(H1)*E12m*S2;
SAT22 = b1*e2*inv(H2)*S2'*E21 + b2*e2*inv(H2)*E21*S2;
SAT21 = -b1*e2*inv(H2)*S2'*E21m - b2*e1*inv(H2)*E21m*S1;

% interface condition 2 & 3
SAT22 = SAT22 + a1*e2*inv(H2)*S2'*E2n + a2*e2*inv(H2)*E2n*S2;
SAT23 = -a1*e2*inv(H2)*S2'*E23m - a2*e3*inv(H2)*E23m*S3;
SAT33 = b1*e3*inv(H3)*S3'*E31 + b2*e3*inv(H3)*E31*S3;
SAT32 = -b1*e3*inv(H3)*S3'*E32m - b2*e2*inv(H3)*E32m*S2;

SAT = [SAT11+SATB1, SAT12, sparse(n_mesh_neg+1,n_mesh_pos+1)
    SAT21, SAT22, SAT23
    sparse(n_mesh_pos+1,n_mesh_neg+1),SAT32,SAT33+SATB3];
% SBP operator with the SAT boundary treatment
Laplacian_cl = D2+SAT;