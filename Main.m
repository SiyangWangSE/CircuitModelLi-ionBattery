% Example script to simulate a discharge process
close all
clear
clc

% initial condition
soc_init_neg = 0.9;
soc_init_pos = 0.1;

% time and voltage limations
time_max = 3600*4; % seconds
V_max = 4.2; % cut-over voltage
V_min = 3; % cut-off voltage

time_step = 10;
order = 4;
I = 20; %A/m^2

Setup_parameters; % The cell geometry, material properties, meshing are defined here
SpatialDiscretization; % create the spatial discretization
Initialize;

tic
Simulation_loop;
toc

%% plot
figure
hold on
plot(Ah_time,V_batt_time,'LineWidth',2);
grid on
xlabel('Charge (mAh/cm^2)')
ylabel('Voltage (V)')

