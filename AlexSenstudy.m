%Ben Kramer
% LU Aero, DBF '212-'23 Sensitivity Study
%9/14/22

%Succ the Zucc

%% Parameters Defined

clc; clear; close all;

%General Parameters/constants
rho_air = .0765 * (1/12)^3; %lb/in^3
g = 32.2; %ft/s^2
e = .8;

c_do = 0.00679; %zero-lift drag coeff for NACA 63-210


%Extreme Parameters
bounding_box = 62; %in
max_gross_weight = 50; %lbs
electronic_pachage_min_dims = [3,3,6]; %in [length, width, height]

b_max = 62; %in
c_max = 7.75; %in
MTOW = 9; %Max takeoff weight - lbs

L_D_max = 3; %Assumed best L/D
C_L_cruise = .5; %Lift Coefficient at cruise

eta_prop = .6; %Propeller efficiency

Energy_max = 100; %Max energy [W-hr]
Energy_max_W_s = 100/60; %Max energy [W-s]


%Actual Parameters of the plane
D_prop = 12; % Propeller diameter [in]
%v_prop = 880; %Velocity at propeller [in/s]
A_prop = (D_prop/2)^2*pi; %Effective propellor area [in^2]

b = b_max; %in
c = c_max; %in
S = b*c; %in^2
AR = b/c;

k = 1/(pi*e*AR);




%% Mission 1 - Takeoff, 3 laps, landing (no payload)

%time_max = 5*60; %Max time [s]
%P_prop = Energy_max_W_s * time_max %Propeller Power [W]
%P_prop_hp = 
%thrust_static = k*(2*rho_air*A_prop*P_prop^2)^(1/3); %Propeller Thrust []



f0 = 3100; %Drag Area (no payload) [in^2]
v = 880; % Cruising Velocity [in/s]
q = (1/2)*rho_air*v^2; %Dynamic pressure [lb/(in-s^2)], 

L = q*S*C_L_cruise; %Lift [(lb-in)/s^2]
D = f0*q+(L^2)/(pi*q*b^2); %Drag [lb-in/s^2]

T = D; %Thrust Needed for constant flight [lb-in/s^2]
P = (T^3/(k*2*rho_air*A_prop)); 
