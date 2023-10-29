% LU Aero DBF 2023-24 Sensitivity Study

clc;
clear; 
close all;

%MTOW = ; %lbs
rho = 0.002377
b = 2.5 * 12; %in from wingspan feet
c = 6 %in from chord length
S = b * c;
AR = 6; %Typically from 5-6
e = 0.8; %Oswald Eff. Factor
k = 1/(pi * e * AR); %K-factor
c_do = 0.4 %zero-lift drag coeff from zero AOA airfoil
g = 32.2; %ft/s^2

n_prop = 0.83 %prop eff. 
P_a = 2 * 550 * n_prop; %Power av.

%n = L/W;
n_g = 4; %Place holder for now

m2_maxT = 5 * 60; %Five minute in seconds
m3_maxT = 5 * 60; 

N_packages_max = 7;
N_packages = [1:7];
W_package = 1;
W_tot = 0;
V_max = 0;
P = 0;
P_round = 0;
Lap_times = zeros(1,7);
total_time = zeros(1,7);
time_sum = 0;
actual_max_laps = 0;

coeff1 = .5 * rho * S * c_do;

battery_capacity = 1; % TODO:
total_flight_time = 5 * 60; % 5 minutes
weight_plane_empty = 5; % pounds

m3 = zeros(1, N_packages_max);
for count = 1:1:N_packages_max
    % compute weight for `count` passangers
    W_tot = weight_plane_empty + count * W_package;
    coeff2 = 2 * k * W_tot^2 / (rho * S);

    % Trust function P w.r.t. speed v
    % h(v)  = P(v) - P_a
    h = @(v) coeff1 * v^3 + coeff2 / v - P_a;
    V_max = fzero(h, 4);
    fprintf('V_max = %f\n', V_max);

    % compute minimum turning radius
    R_min = (V_max^2) / (g * sqrt(n_g^2-1));
    fprintf('min_R = %f\n', R_min);

    % 
    Lap_distance = 2000 + 4 * R_min * pi;
    Lap_time = Lap_distance / V_max;
    fprintf('Lap time = %f\n', Lap_time);

    number_of_laps = floor(total_flight_time / Lap_time)
    m3(count) = 2 + number_of_laps * count / battery_capacity
end

min_obj_weight = 0.5; % minimum weight change, i.e., by adding an object
number_of_laps_m2 = 3;

weights = weight_plane_empty:min_obj_weight:max_weight;
m2 = zeros(size(weights));
i = 1;
for weight = weights
    coeff2 = 2 * k * weight^2 / (rho * S);

    % Trust function P w.r.t. speed v
    % h(v)  = P(v) - P_a
    h = @(v) coeff1 * v^3 + coeff2 / v - P_a;
    V_max = fzero(h, 4);
    fprintf('V_max = %f\n', V_max);

    % compute minimum turning radius
    R_min = (V_max^2) / (g * sqrt(n_g^2-1));
    fprintf('min_R = %f\n', R_min);

    % 
    Lap_distance = 2000 + 4 * R_min * pi;
    Lap_time = Lap_distance / V_max;
    %fprintf('Lap time = %f\n', Lap_time);

    time = number_of_laps_m2 * Lap_time;
    payload_weight = weight - weight_plane_empty;

    m2(i) = 1 + payload_weight / time;
    i = i + 1; 
end