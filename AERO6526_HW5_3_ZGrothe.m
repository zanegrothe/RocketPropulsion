% Zane Grothe
% AERO 6526
% HW 5
% 3/31/21

% Problem #3

clear all
close all
clc

% Givens----------

cf=0.0035;
cp=20.8;
Pr=0.7;
gam=1.24;
hg=4600;
Me=3;
k=5;

% Equations----------

M=linspace(0,3.5);
rho_r=(1+(gam-1)/2*M.^2).^(-1/(gam-1));
v_r=M/Me.*(1+(gam-1)/2*M.^2).^(-1/2);

plot(M,rho_r,M,v_r,M,rho_r.*v_r.*k)
xlabel('Mach Number')
title('Nozzle Heat Transfer Estimation')
legend('rho/rho0','V/Vmax','K(rho)(V) Heat Transfer Rate')


