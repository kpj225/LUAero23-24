% LU Aero DBF 2023-24 Sensitivity Study

clc;
clear; 
close all;

%MTOW = ; %lbs
rho = 0.002377
%b = %in from wingspan feet
%c = %in from chord length
S = b * c;
AR = 6; %Typically from 5-6
e = 0.8; %Oswald Eff. Factor
k = 1/(pi * e * AR); %K-factor
%c_do = %zero-lift drag coeff from zero AOA airfoil
g = 32.2; %ft/s^2
