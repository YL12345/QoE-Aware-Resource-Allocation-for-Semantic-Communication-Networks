clc 
clear

% comparison between the proposed mathod and single-cell optimization method
% N_cell=3, N_channels=6
% x-axis: (N_S,N_Bi):(6,12),(6,18),(9,18)
% Monte carlo times: 500

figure;
y=[14.0847,13.2116; 15.7155,14.3883; 16.4691,15.0225];
b=bar(y);
