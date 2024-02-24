% Zane Grothe
% AERO 6526
% HW 3
% 3/1/21

% Problem #2

clear all
close all
clc


% Givens
a=.0285; %(in/s)/psi^n
n=.35;
rho=.064; %lbm/in^3
gam=1.24;
To=6156; %Rankine
MW=24.7;
Ru=1545; %ft*lbf/(lbm*mole*Rankine)
g=32.174; %ft/s^2
pa=14.7; %psi

% Assumptions
zeta=.75;
phi=.8;


% Original motor dimentions
Lm=100; %in
Dm=16.7; %in
Dt=6.82; %in


% Part 1 ----------------------

% 1/2 motor dimentions (h=half)
Lm_h=Lm/2; %in
Dm_h=Dm/2; %in
Dt_h=Dt/2; %in

Vm_h=Lm_h*pi*(Dm_h^2)/4; %in^3
mp_h=Vm_h*rho*phi; %lbm (mass of propellant)

Wp_h=mp_h; %lbf (weight of propellant)
Wo_h=Wp_h/zeta; %lbf (initial weight)
mo_h=Wo_h; %lbm (initial mass)

De_h=Dm_h; %in
At_h=pi*Dt_h^2/4; %in^2
Ae_h=pi*De_h^2/4; %in^2
AeAt_h=Ae_h/At_h; % (area ratio)

Ab_h=Ae_h; %in^2
cstar_h=sqrt([((gam+1)/2)^((gam+1)/(gam-1))]*Ru*g*To/MW/gam); %ft/s (characteristic velocity)

syms M1
AA1=(1/M1)*[2/(gam+1)*((1+(gam-1)/2*M1^2))]^((gam+1)/2/(gam-1))==AeAt_h;
Me_h=double(solve(AA1,M1)); % (exit Mach number)
po_h=(cstar_h*Ab_h*a*rho/At_h/32.174)^(1/(1-n));
pe_h=po_h/[(1+(gam-1)/2*Me_h^2)^(gam/(gam-1))]; %psi (chamber pressure)

U1=(2*(gam^2)/(gam-1));
V1=((2/(gam+1))^((gam+1)/(gam-1)));
W1=(1-((pe_h/po_h)^((gam-1)/gam)));
Cf_h=sqrt(U1*V1*W1)+AeAt_h*(pe_h/po_h-pa/po_h); % (thrust coefficient)

Ft_h=po_h*At_h*Cf_h; %lbf (thrust)
Isp_h=Cf_h*cstar_h/g; %seconds (specific impulse)
I_h=Wp_h*Isp_h; %lbf-seconds (total impulse)
ta_h=I_h/Ft_h; %seconds (action time)


% Part 2 ----------------------

% Original motor dimentions
Vm=Lm*pi*(Dm^2)/4; %in^3
mp=Vm*rho*phi; %lbm (mass of propellant)

Wp=mp; %lbf (weight of propellant)
Wo=Wp/zeta; %lbf (initial weight)
mo=Wo; %lbm (initial mass)

De=Dm; %in
At=pi*Dt^2/4; %in^2
Ae=pi*De^2/4; %in^2
AeAt=Ae/At; % (area ratio)

Ab=Ae; %in^2
cstar=sqrt([((gam+1)/2)^((gam+1)/(gam-1))]*Ru*g*To/MW/gam); %ft/s (characteristic velocity)

syms M2
AA2=(1/M2)*[2/(gam+1)*((1+(gam-1)/2*M2^2))]^((gam+1)/2/(gam-1))==AeAt;
Me=double(solve(AA2,M2)); % (exit Mach number)
po=(cstar*Ab*a*rho/At/32.174)^(1/(1-n));
pe=po/[(1+(gam-1)/2*Me^2)^(gam/(gam-1))]; %psi (chamber pressure)

U2=(2*(gam^2)/(gam-1));
V2=((2/(gam+1))^((gam+1)/(gam-1)));
W2=(1-((pe/po)^((gam-1)/gam)));
Cf=sqrt(U2*V2*W2)+AeAt*(pe/po-pa/po); % (thrust coefficient)

Ft=po*At*Cf; %lbf (thrust)
Isp=Cf*cstar/g; %seconds (specific impulse)
I=Wp*Isp; %lbf-seconds (total impulse)
ta=I/Ft; %seconds (action time)


% Part 3 ----------------------

% Double motor dimentions (d=double)
Lm_d=Lm*2; %in
Dm_d=Dm*2; %in
Dt_d=Dt*2; %in

Vm_d=Lm_d*pi*(Dm_d^2)/4; %in^3
mp_d=Vm_d*rho*phi; %lbm (mass of propellant)

Wp_d=mp_d; %lbf (weight of propellant)
Wo_d=Wp_d/zeta; %lbf (initial weight)
mo_d=Wo_d; %lbm (initial mass)

De_d=Dm_d; %in
At_d=pi*Dt_d^2/4; %in^2
Ae_d=pi*De_d^2/4; %in^2
AeAt_d=Ae_d/At_d; % (area ratio)

Ab_d=Ae_d; %in^2
cstar_d=sqrt([((gam+1)/2)^((gam+1)/(gam-1))]*Ru*g*To/MW/gam); %ft/s (characteristic velocity)

syms M3
AA3=(1/M3)*[2/(gam+1)*((1+(gam-1)/2*M3^2))]^((gam+1)/2/(gam-1))==AeAt_d;
Me_d=double(solve(AA3,M3)); % (exit Mach number)
po_d=(cstar_d*Ab_d*a*rho/At_d/32.174)^(1/(1-n));
pe_d=po_d/[(1+(gam-1)/2*Me_d^2)^(gam/(gam-1))]; %psi (chamber pressure)

U3=(2*(gam^2)/(gam-1));
V3=((2/(gam+1))^((gam+1)/(gam-1)));
W3=(1-((pe_d/po_d)^((gam-1)/gam)));
Cf_d=sqrt(U3*V3*W3)+AeAt_d*(pe_d/po_d-pa/po_d); % (thrust coefficient)

Ft_d=po_d*At_d*Cf_d; %lbf (thrust)
Isp_d=Cf_d*cstar_d/g; %seconds (specific impulse)
I_d=Wp_d*Isp_d; %lbf-seconds (total impulse)
ta_d=I_d/Ft_d; %seconds (action time)


% Part 4 ----------------------

% Collect data and plot
Dmm=[Dm_h,Dm,Dm_d];

mpm=[mp_h,mp,mp_d];
mom=[mo_h,mo,mo_d];
pom=[po_h,po,po_d];
Ftm=[Ft_h,Ft,Ft_d];
Cfm=[Cf_h,Cf,Cf_d];
tam=[ta_h,ta,ta_d];
Ispm=[Isp_h,Isp,Isp_d];
Im=[I_h,I,I_d];
cstarm=[cstar_h,cstar,cstar_d];

figure(1)
plot(Dmm,mpm)
xlabel('Body Diameter (in)')
ylabel('Propellant Mass (lbm)')
title('Propellant Mass vs Body Diameter')

figure(2)
plot(Dmm,mom)
xlabel('Body Diameter (in)')
ylabel('Initial Mass (lbm)')
title('Initial Mass vs Body Diameter')

figure(3)
plot(Dmm,pom)
xlabel('Body Diameter (in)')
ylabel('Chamber Pressure (psi)')
title('Chamber Pressure vs Body Diameter')

figure(4)
plot(Dmm,Ftm)
xlabel('Body Diameter (in)')
ylabel('Thrust (lbf')
title('Thrust vs Body Diameter')

figure(5)
plot(Dmm,Cfm)
xlabel('Body Diameter (in)')
ylabel('Thrust Coefficient')
title('Thrust Coefficient vs Body Diameter')

figure(6)
plot(Dmm,tam)
xlabel('Body Diameter (in)')
ylabel('Action Time (seconds)')
title('Action Time vs Body Diameter')

figure(7)
plot(Dmm,Ispm)
xlabel('Body Diameter (in)')
ylabel('Specific Impulse (seconds)')
title('Specific Impulse vs Body Diameter')

figure(8)
plot(Dmm,Im)
xlabel('Body Diameter (in)')
ylabel('Total Impulse (lbf-seconds)')
title('Total Impulse vs Body Diameter')

figure(9)
plot(Dmm,cstarm)
xlabel('Body Diameter (in)')
ylabel('Characteristic Velocity (ft/s)')
title('Characteristic Velocity vs Body Diameter')






% I didn't know how to find chamber pressure. I incorrectly assumed an
% end burning grain to get a po, so I could move forward with the problem.
% Hopefully everything after that is calculated and plotted correctly (just
% with the wrong values).