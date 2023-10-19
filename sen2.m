% LU Aero, DBF '21-'22 Sensitivity Study
clc;
clear;
close all;

max_gross_weight = 15; %lbs
rho = .002377;
b = 8*12; %in
c = 16; %in
S = b*c; %in^2
AR = 6;
e = .8;
k = 1/(pi*e*AR);
c_do = 0.4; %zero-lift drag coeff
n_prop = 0.83; %prop efficiency
P_a = 3*4400*n_prop; %power available (W), Avian 8085-160Kv Outrunner Brushless Motor
g = 32.2;
n_g = 4;
m3_maxT = 10*60;
m2_maxT = 5*60;
land_tO_time = 10;
deploy_time = 10;
N_packages_max = 8;
N_packages = [1:8];
W_package = 8/16;
W_tot = 0;
V_max = 0;
P = 0;
P_round = 0;
Lap_times = zeros(1,8);
total_time = zeros(1,8);
time_sum = 0;
count = 1;
actual_max_laps = 0;

coeff1 = .5*rho*S*c_do;
coeff2 = (2*k*((max_gross_weight)^2)/(rho*S));

    
for i = 1:1:N_packages_max

         for j = 1:0.00001:40
             P = coeff1*(j^3) + coeff2*(1/j);
             P_round = round(P);
             if(P_round == P_a)
                 V_max = j;
             end
         end

        fprintf('V_max = %f\n', V_max);
        R_min = (V_max^2)/(g*sqrt((n_g^2)-1));
        fprintf('min_R = %f\n', R_min);
        Lap_distance = 2000+4*(R_min*pi);
        Lap_time = Lap_distance/V_max;
        %fprintf('Lap time = %f\n', Lap_time);

        m3_nominal_lap = land_tO_time + Lap_time + deploy_time;
        Lap_times(count) = m3_nominal_lap;
        time_sum = time_sum + m3_nominal_lap;
        total_time(count) = time_sum;
        W_tot = max_gross_weight - ((count-1)*W_package);
        coeff2 = (2*k*((W_tot)^2)/(rho*S));
        count = count +1;
    if time_sum <= m3_maxT
        actual_max_laps = actual_max_laps +1;
    end
end

fprintf('Max laps = %f\n', actual_max_laps);
figure
hold on
plot(total_time, N_packages)
xline(m3_maxT)
yline(actual_max_laps)
xlabel('Time (s)')
ylabel('Packages deployed')
ylim([0 8])
title('Packages deployed vs Time')
hold off

N_syringes = 10*N_packages;
min_syringes = 10*actual_max_laps;
W_payload = actual_max_laps*W_package;
W_empty = max_gross_weight - W_payload;
W_syringe = 20/454;
max_syringes = floor(W_payload/W_syringe);
fprintf('Max # of syringes = %f\n', max_syringes);

W_tot2 = W_empty;
time_sum2 = 0;
count2 = 1;
coeff22 = (2*k*((max_gross_weight)^2)/(rho*S));
P2 = 0;
V_max2 = 0;
P_round2 = 0;
R_min2 = 0;
Lap_distance2 = 0;
Lap_time2 = 0;
cols = length(10:max_syringes);
m2_times = zeros(1, cols);

m2_points = [];
m3_points = [];
actual_max_syringes = 9;
for i = 10:1:max_syringes
        W_tot2 = W_empty + ((i)*W_syringe);
        coeff22 = (2*k*((W_tot2)^2)/(rho*S));

         for j = 1:0.00001:25
             P2 = coeff1*(j^3) + coeff22*(1/j);
             P_round2 = round(P2);
             if(P_round2 == P_a)
                 V_max2 = j;
             end
         end
        %fprintf('V_max2 = %f\n', V_max2);
      
        R_min2 = (V_max2^2)/(g*sqrt((n_g^2)-1));
        
        Lap_distance2 = 2000+4*(R_min2*pi);
        %fprintf('Lap distance = %f\n', Lap_distance2);
        Lap_time2 = Lap_distance2/V_max2;
       % fprintf('Lap time = %f\n', Lap_time2);
        time_sum2 = land_tO_time - 5 + 3*Lap_time2;
        m2_times(count2) = time_sum2;
        %fprintf('Mission time = %f\n', time_sum2);
       
    if time_sum2 <= m2_maxT
        actual_max_syringes = actual_max_syringes +1;
    end
    m2_points(count2) = actual_max_syringes/time_sum2;
    m3_points(count2) = floor(actual_max_syringes/10);
    count2 = count2 +1;
end
%fprintf('Actual max syringes = %f\n', actual_max_syringes);

m2_syringes = 10:actual_max_syringes;

figure
hold on
plot(m2_syringes, m2_points, '-r')
plot(m2_syringes, m3_points, '-b')
xlabel('Payload Size (# of syringes)')
ylabel('Mission Points')
legend('Mission 2', 'Mission 3') 
title('Mission Points vs Payload')
hold off

percent_del_syringe = [];
percent_del_mission2 = [];
percent_del_mission3 = [];
m3_points_adjust = linspace(1,5,47);
for i = 1:(length(m2_syringes)-1)
%     if(i==1)
%         percent_del_syringe(i) = m2_syringes(1)/10;
%         percent_del_mission2(i) = m2_points(1)/1;
%         percent_del_mission3(i) = m3_points(1)/2;
%     else
        percent_del_syringe(i) = 100*(m2_syringes(i+1)- m2_syringes(i))/m2_syringes(i);
        percent_del_mission2(i) = 100*(m2_points(i+1)- m2_points(i))/m2_points(i);
        percent_del_mission3(i) = 100*(m3_points_adjust(i+1)- m3_points_adjust(i))/m3_points_adjust(i);
%     end
end

figure
hold on
plot(percent_del_syringe, percent_del_mission2, '--r')
plot(percent_del_syringe, percent_del_mission3, '-b')
xlabel('% Change in Score Variable')
ylabel('% Change in Final Score')
legend('M2', 'M3') 
%title('Mission Points vs Payload')
hold off