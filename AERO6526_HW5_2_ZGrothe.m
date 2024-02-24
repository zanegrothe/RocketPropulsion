% Zane Grothe
% AERO 6526
% HW 5
% 3/31/21

% Problem #2

clear all
close all
clc

% Givens----------

At=10; %in^2
eps=20; %Ae/At
xe=20; %inches (length of nozzle)
theta_n=28*pi/180; %radians (pulled from graph for Lf=86% and eps=20)

% Calculate Parameters----------

Ae=eps*At; %in^2
rt=sqrt(At/pi); %inches
re=sqrt(Ae/pi); %inches
ra=.382*rt; %inches
sig=ra/rt;

% Equations----------

xp=sig*rt*sin(theta_n);
yp=rt*(1+sig)-sig*rt*cos(theta_n);
k1=xe-xp;
k2=2*tan(theta_n);
ro=(yp^2+yp*k1*k2-re^2)/(2*yp+k1*k2-2*re);
b=xp-(yp-ro)/k2;
a=4*(xp-b)*tan(theta_n)^2;

x=linspace(xp,xe,50);
r=ro+sqrt(a*(x-b)); % parabolic nozzle
z=linspace(0,xp,25);
rr=-sqrt(0.4624-z.^2)+2.46; % arc throat

% Plotting----------

plot(x,r)
hold on
plot(z,rr)
xlim([0,xe+xe/10])
ylim([0,max(r)+max(r)/10])
xlabel('x (inches)')
ylabel('r (inches, from slip line)')
title('Nozzle Wall Boundary')

