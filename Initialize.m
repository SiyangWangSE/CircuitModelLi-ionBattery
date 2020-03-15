% Initialize the battery
%==========================================================================
% electrolyte
c_l = ones(1,n_mesh_neg+n_mesh_sep+n_mesh_pos+3)*1000; % one vector (row)
phi_l = zeros(1,n_mesh_neg+n_mesh_sep+n_mesh_pos+3);   % one vector (row)

%==========================================================================
% negative electrode
c_s_neg = ones(n_mesh_s+1,n_mesh_neg+1)*c_s_max_neg*soc_init_neg; % particle domain meshing -by- electrolyte domain meshing
c_s_surface_neg = c_s_neg(end,:); % the first row is the center concentration and the last row is the surface concentration

phis_eq_neg= ones(1,n_mesh_neg+1)*pchip(theta_lookup_neg, OCV_lookup_neg, soc_init_neg); % based on the average concentration
phis_U_neg = phis_eq_neg; % based on the suurface concentration

R_s_neg = h_neg/(eps_s_neg^Bruggeman_neg*sigma_s_neg)*ones(1,n_mesh_neg); % per electrode area
R_l_neg = h_neg/(eps_l_neg^Bruggeman_neg*sigma_l)*ones(1,n_mesh_neg); % per electrode area

%linear B-V
i_0_neg = F*k_ct_neg*(c_s_max_neg-c_s_surface_neg).^0.5.*c_s_surface_neg.^0.5.*c_l(1:n_mesh_neg+1).^0.5;
R_ct_neg = R*T./(F*i_0_neg)/(Sa_neg*h_neg); % per electrode area

% initialize current distribution
for i = 1:10
    R_x_neg = spdiags([-R_ct_neg(2:end)' (R_s_neg + R_ct_neg(1:end-1)+R_ct_neg(2:end)+R_l_neg)' -R_ct_neg(1:end-1)'],-1:1,n_mesh_neg,n_mesh_neg);
    V_x_neg = phis_U_neg(1:end-1)-phis_U_neg(2:end)+I*R_l_neg;
    V_x_neg(1) = V_x_neg(1)+I*R_ct_neg(1);
    I_x_neg = R_x_neg\V_x_neg';
    
    I_ct_neg = [-diff([I;I_x_neg]);I_x_neg(end)]'; % charge transfer current per electrode area
    I_loc_neg = I_ct_neg./(Sa_neg*volume_length_neg); % charge transfer current per active area
    eta_neg = 2*R*T/F*asinh(I_loc_neg./2./i_0_neg);
    R_ct_neg = eta_neg./I_ct_neg;
end
%==========================================================================
% positive electrode
c_s_pos = ones(n_mesh_s+1,n_mesh_pos+1)*c_s_max_pos*soc_init_pos; % particle domain meshing -by- electrolyte domain meshing
c_s_surface_pos = c_s_pos(end,:); % the first row is the center concentration and the last row is the surface concentration

phis_eq_pos= ones(1,n_mesh_pos+1)*pchip(theta_lookup_pos, OCV_lookup_pos, soc_init_pos); % based on the average concentration
phis_U_pos = phis_eq_pos; % based on the suurface concentration

R_s_pos = h_pos/(eps_s_pos^Bruggeman_pos*sigma_s_pos)*ones(1,n_mesh_pos); % per electrode area
R_l_pos = h_pos/(eps_l_pos^Bruggeman_pos*sigma_l)*ones(1,n_mesh_pos); % per electrode area

%linear B-V
i_0_pos = F*k_ct_pos*(c_s_max_pos-c_s_surface_pos).^0.5.*c_s_surface_pos.^0.5.*c_l(n_mesh_neg+n_mesh_sep+3:end).^0.5;
R_ct_pos = R*T./(F*i_0_pos)/(Sa_pos*h_pos); % per electrode area

% initialize current distribution
for i = 1:10
    R_x_pos = spdiags([-R_ct_pos(2:end)' (R_s_pos + R_ct_pos(1:end-1)+R_ct_pos(2:end)+R_l_pos)' -R_ct_pos(1:end-1)'],-1:1,n_mesh_pos,n_mesh_pos);
    V_x_pos = phis_U_pos(1:end-1)-phis_U_pos(2:end)+I*R_l_pos;
    V_x_pos(end) = V_x_pos(end)+I*R_ct_pos(end);
    I_x_pos = R_x_pos\V_x_pos';
    
    I_ct_pos = [I_x_pos(1);diff([I_x_pos;I])]'; % charge transfer current per electrode area
    I_loc_pos = I_ct_pos./(Sa_pos*volume_length_pos); % charge transfer current per active area
    eta_pos = 2*R*T/F*asinh(I_loc_pos./2./i_0_pos);
    R_ct_pos = eta_pos./I_ct_pos;
end


