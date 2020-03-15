% setup parameters
%==========================================================================
% constant
F = 96485;
R = 8.3145;
T = 293.15;

%==========================================================================
% geometry ----- graphite | sep | NMC
Area = 1; % m^2
L_neg = 46.6e-6; % um, thickness
L_sep = 18.7e-6; % um, thickness
L_pos = 43e-6; % um, thickness
r_s_neg = 6.3e-6; % um, particle radius
r_s_pos = 2.13e-6; % um, particle radius

%==========================================================================
% meshing
n_mesh_neg = 8;
n_mesh_sep = 8;
n_mesh_pos = 8;
n_mesh_s = 30; % mesh inside the particle, same for both electrode for now

h_neg = L_neg/n_mesh_neg;
h_sep = L_sep/n_mesh_sep;
h_pos = L_pos/n_mesh_pos;
h_electrolyte = [ones(1,n_mesh_neg+1)*h_neg ones(1,n_mesh_sep+1)*h_sep ones(1,n_mesh_pos+1)*h_pos];
volume_length_neg = ones(1,n_mesh_neg+1)*h_neg;
volume_length_neg(1) = h_neg/2;
volume_length_neg(end) = h_neg/2;
volume_length_pos = ones(1,n_mesh_pos+1)*h_pos;
volume_length_pos(1) = h_pos/2;
volume_length_pos(end) = h_pos/2;

h_s_neg = r_s_neg/n_mesh_s;
h_s_pos = r_s_pos/n_mesh_s;

r_distance_neg = 0:h_s_neg:r_s_neg;
r_distance_pos = 0:h_s_pos:r_s_pos;
%==========================================================================
% porosity related parameter
eps_l_neg = 0.29; % porosity
eps_l_sep = 0.4; % porosity
eps_l_pos = 0.21; % porosity
eps_l = [ones(1,n_mesh_neg+1)*eps_l_neg  ones(1,n_mesh_sep+1)*eps_l_sep ones(1,n_mesh_pos+1)*eps_l_pos];

eps_s_neg = 0.49; % volume fraction of the active material
eps_s_pos = 0.57; % volume fraction of the active material

Sa_neg = 3*eps_s_neg/r_s_neg; % spacific surface area m^2/m^3^, assume sphere shape particle
Sa_pos = 3*eps_s_pos/r_s_pos; % spacific surface area m^2/m^3^, assume sphere shape particle

Bruggeman_neg = 1.52; %Bruggeman constant
Bruggeman_sep = 1.62;
Bruggeman_pos = 1.44;
Bruggeman = [ones(1,n_mesh_neg+1)*Bruggeman_neg  ones(1,n_mesh_sep+1)*Bruggeman_sep ones(1,n_mesh_pos+1)*Bruggeman_pos];
%==========================================================================
% negative solid material parameter
sigma_s_neg = 100;  % conductivity S/m
D_s_neg = 3e-14;     % diffusion coefficient m^2/s
c_s_max_neg = 31390; % maximum concentration mol/m3
k_ct_neg = 2e-11;   % rate constant

% OCV curve
filename = 'COMSOL_Graphite_OCV.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, '%f%f%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
theta_lookup_neg = dataArray{1};
OCV_lookup_neg = dataArray{2};

%==========================================================================
% positive solid material parameter
sigma_s_pos = 10;  % conductivity S/m
D_s_pos = 5e-15;     % diffusion coefficient m^2/s
c_s_max_pos = 48390; % maximum concentration mol/m3
k_ct_pos = 2e-11;   % rate constant

% OCV curve
filename = 'COMSOL_NMC_OCV.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, '%f%f%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
theta_lookup_pos = dataArray{1};
OCV_lookup_pos = dataArray{2};

%==========================================================================
% electrolyte parameter
sigma_l = 1;               % conductivity S/m
D_l =2.5e-10;                  % diffusion coefficient m^2/s
tplus = 0.26;               % transport number of lithium ion
f_a = 0.2;                   % activity coefficient

sigma_l_eff = eps_l.^Bruggeman.*sigma_l; % effective conductivity S/m
D_l_eff = eps_l.^Bruggeman*D_l; % effective diffusion coefficient of the electrolyte m^2/s
