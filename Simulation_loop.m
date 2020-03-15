% main simulation loop

time_now = 0;

% initialize data storage
time = 0:time_step:time_max; % only valid for a constant time step
Ah_time = zeros(time_max/time_step+1,1);
V_batt_time = zeros(time_max/time_step+1,1);V_batt_time(1) = phis_eq_pos(end)-phis_eq_neg(1);
c_s_surface_to_compare = zeros(time_max/time_step+1,n_mesh_neg+n_mesh_pos+2);c_s_surface_to_compare(1,:) = [c_s_surface_neg c_s_surface_pos];
c_l_to_compare=zeros(time_max/time_step+1,n_mesh_neg+n_mesh_sep+n_mesh_pos+3);c_l_to_compare(1,:) = c_l;

k = 1;
while time_now < time_max
    
    %======================================================================
    % calculate electrolyte concentration
    R_li = [I_ct_neg./volume_length_neg/F  zeros(1,n_mesh_sep+1)  -I_ct_pos./volume_length_pos/F]; % lithium ion generation per volume (electrode area)
    op = odeset('reltol',1e-3);
    [t,c_l_temp] = ode15s(@(t,c_l_temp) ((D_l.*Laplacian_cl*c_l_temp+R_li'*(1-tplus))./eps_l'), [0 time_step], c_l', op);
    c_l = c_l_temp(end,:); % update
    
    %======================================================================
    % calculate negative electrode concentratoin
    for n_particle = 1:n_mesh_neg+1
        c_s_particle = c_s_neg(:,n_particle); % distribution inside one particle
        g_neg(end) = -I_loc_neg(n_particle)/D_s_neg/F; % boundry condition
        [t,c_s_temp] = ode15s(@(t,c_s_temp)(D_s_neg*(Laplacian_cs_neg*c_s_temp+invHB_cs_neg*g_neg)), [0 time_step],c_s_particle,options_c_s_neg);
        c_s_neg(:,n_particle) = c_s_temp(end,:)';
    end
    substance = sum((c_s_neg(1:end-1,:)+c_s_neg(2:end,:))/2*4/3*pi.*(r_distance_neg(2:end).^3-r_distance_neg(1:end-1).^3)');
    c_s_average_neg = substance/(4/3*pi*r_s_neg^3); % only for data analysis, not needed for simulation
    c_s_surface_neg = c_s_neg(end,:);
    phis_eq_neg= pchip(theta_lookup_neg, OCV_lookup_neg, c_s_average_neg/c_s_max_neg); % based on the average concentration, only for data analysis, not needed for simulation
    phis_U_neg = pchip(theta_lookup_neg, OCV_lookup_neg, c_s_surface_neg/c_s_max_neg); % based on the suurface concentration
    
    %======================================================================
    % calculate positive electrode concentratoin
    for n_particle = 1:n_mesh_pos+1
        c_s_particle = c_s_pos(:,n_particle); % distribution inside the particle
        g_pos(end) = I_loc_pos(n_particle)/D_s_pos/F; % boundry condition
        [t,c_s_temp] = ode15s(@(t,c_s_temp)(D_s_pos*(Laplacian_cs_pos*c_s_temp+invHB_cs_pos*g_pos)), [0 time_step],c_s_particle,options_c_s_pos);
        c_s_pos(:,n_particle) = c_s_temp(end,:)';
    end
    substance = sum((c_s_pos(1:end-1,:)+c_s_pos(2:end,:))/2*4/3*pi.*(r_distance_pos(2:end).^3-r_distance_pos(1:end-1).^3)');
    c_s_average_pos = substance/(4/3*pi*r_s_pos^3); % only for data analysis, not needed for simulation
    c_s_surface_pos = c_s_pos(end,:);
    phis_eq_pos= pchip(theta_lookup_pos, OCV_lookup_pos, c_s_average_pos/c_s_max_pos); % based on the average concentration, only for data analysis, not needed for simulation
    phis_U_pos = pchip(theta_lookup_pos, OCV_lookup_pos, c_s_surface_pos/c_s_max_pos); % based on the suurface concentration
    %======================================================================
    % calculate electrolyte potential
    I_e = [I-I_x_neg;ones(n_mesh_sep+3,1)*I;I-I_x_pos]';
    % take the electrolyte at the negative current collector as the
    % reference
    phi_lx = (I_e-2*sigma_l_eff.*R*T/F*(1+f_a)*(1-tplus).*[gradient(c_l(1:n_mesh_neg+1)) gradient(c_l(n_mesh_neg+2:n_mesh_neg+n_mesh_sep+2)) gradient(c_l(n_mesh_neg+n_mesh_sep+3:end))]./h_electrolyte./c_l)./(-sigma_l_eff);      
    phi_l(1:n_mesh_neg+1) = cumtrapz(phi_lx(1:n_mesh_neg+1)).*h_neg;
    phi_l(n_mesh_neg+2:n_mesh_neg+n_mesh_sep+2) = cumtrapz(phi_lx(n_mesh_neg+2:n_mesh_neg+n_mesh_sep+2)).*h_sep+ phi_l(n_mesh_neg+1);
    phi_l(n_mesh_neg+n_mesh_sep+3:end) = cumtrapz(phi_lx(n_mesh_neg+n_mesh_sep+3:end)).*h_pos+ phi_l(n_mesh_neg+n_mesh_sep+2);
    
    % calculate electrode potential
    phis_neg = phi_l(1:n_mesh_neg+1)+phis_U_neg+I_ct_neg.*R_ct_neg;
    phis_pos = phi_l(end-n_mesh_pos:end)+phis_U_pos-I_ct_pos.*R_ct_pos;

    
    %======================================================================
    % update resistance matrix and current distribution
    % R_s_neg and R_s_pos remina unchange
    
    % negative electrode
    R_l_neg = -phi_lx(1:n_mesh_neg).*h_neg./I_e(1:n_mesh_neg);
    i_0_update = F*k_ct_neg*(c_s_max_neg-c_s_surface_neg).^0.5.*c_s_surface_neg.^0.5.*c_l(1:n_mesh_neg+1).^0.5;
    eta_neg = 2*R*T/F*asinh(I_loc_neg./2./i_0_update);
    R_ct_neg = eta_neg./I_ct_neg;
    
    R_x_neg = spdiags([-R_ct_neg(2:end)' (R_s_neg + R_ct_neg(1:end-1)+R_ct_neg(2:end)+R_l_neg)' -R_ct_neg(1:end-1)'],-1:1,n_mesh_neg,n_mesh_neg);
    V_x_neg = phis_U_neg(1:end-1)-phis_U_neg(2:end)+I*R_l_neg;
    V_x_neg(1) = V_x_neg(1)+I*R_ct_neg(1);
    I_x_neg = R_x_neg\V_x_neg';
    
    I_ct_neg = [-diff([I;I_x_neg]);I_x_neg(end)]'; % charge transfer current per electrode area
    I_loc_neg = I_ct_neg./(Sa_neg*volume_length_neg); % charge transfer current per active area
    
    
    % positive electrode
    R_l_pos = -phi_lx(end-n_mesh_pos+1:end).*h_pos./I_e(end-n_mesh_pos+1:end);
    i_0_update = F*k_ct_pos*(c_s_max_pos-c_s_surface_pos).^0.5.*c_s_surface_pos.^0.5.*c_l(n_mesh_neg+n_mesh_sep+3:end).^0.5;
    eta_pos = 2*R*T/F*asinh(I_loc_pos./2./i_0_update);
    R_ct_pos = eta_pos./I_ct_pos;
    
    R_x_pos = spdiags([-R_ct_pos(2:end)' (R_s_pos + R_ct_pos(1:end-1)+R_ct_pos(2:end)+R_l_pos)' -R_ct_pos(1:end-1)'],-1:1,n_mesh_pos,n_mesh_pos);
    V_x_pos = phis_U_pos(1:end-1)-phis_U_pos(2:end)+I*R_l_pos;
    V_x_pos(end) = V_x_pos(end)+I*R_ct_pos(end);
    I_x_pos = R_x_pos\V_x_pos';
    
    I_ct_pos = [I_x_pos(1);diff([I_x_pos;I])]'; % charge transfer current per electrode area
    I_loc_pos = I_ct_pos./(Sa_pos*volume_length_pos); % charge transfer current per active area
    %======================================================================
    % record the data
    k = k+1;
    Ah_time(k) = Ah_time(k-1)+I*time_step/3600*0.1; % mAh/cm2
    V_batt_time(k) = phis_pos(end)-phis_neg(1);
    c_l_to_compare(k,:) = c_l;
    c_s_surface_to_compare(k,:) = [c_s_surface_neg c_s_surface_pos];
    %======================================================================
    % stop criteria
    if V_batt_time(k) > V_max || V_batt_time(k) < V_min
        disp('Maximum or minimum voltage limite is reached')
        Ah_time(k+1:end) = [];
        V_batt_time(k+1:end) = [];
        time(k+1:end) = [];
        break
    end
    
    if any(c_l<1) || any(c_s_surface_neg)<1 || any(c_s_surface_neg)>c_s_max_neg ...
            || any(c_s_surface_pos)<1 || any(c_s_surface_pos)>c_s_max_pos
        disp('Concentration exceed limit')
        Ah_time(k+1:end) = [];
        V_batt_time(k+1:end) = [];
        time(k+1:end) = [];
        break
    end
    
    time_now =time_now + time_step;
end



