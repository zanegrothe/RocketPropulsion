% Zane Grothe
% AERO 6526
% Final Paper
% 4/28/21

clear all
close all
clc

R=linspace(0,5,50);
u=15;

v=50;
while v<=200;
    f=v*sqrt(u*R/4/pi)*60;
    plot(R,f)
    hold on
    v=v+50;
end

xlabel('r/m, s/ft^3')
ylabel('Pump Frequency, Rpm')
text(.25,28000,'v is rotor tip speed')
text(.25,26500,'r is propellant density, lb/ft^3')
text(.25,25000,'m is propellant mass flow rate, lb/s')
text(2,3500,'v=50ft/s')
text(2.5,9000,'v=100ft/s')
text(3,16000,'v=150ft/s')
text(3.5,24000,'v=200ft/s')

