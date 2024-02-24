% Zane Grothe
% AERO 6526
% HW 4
% 3/11/21

% Problem #4

clear all
close all
clc

% Define Givens
GAM=1.6;
J=.25;
L=1;

% Plot p/p2 as a function of distance from head end
x=linspace(0,L);
p_p2=1+(GAM^2)*(J^2)*((L-x)/L);

plot(x,p_p2)
xlabel('Ratio of Distance from Head End to Full Grain Length (x/L)')
ylabel('Ratio of Pressure to Aft End Pressure p(x)/p2')
title('p/p2 as a Funtion of Distance from Head End x')


