% Zane Grothe
% AERO 6526
% HW 3
% 3/1/21

% Problem #5

clear all
close all
clc

% Givens------------------------------------

% Design Parameters
Ro=4; %in
Rp=2.5; %in
f=0.5; %in
epsilon=0.9;
N=9;
L=24; %in
po=800; %psia
Web=3; %in

% Propellant Parameters
a=0.0563; %in/sec/psin
n=0.33;
cstar=5343; %ft/sec
rho=0.064; %lbm/in^3
Tc=6260; %Rankine
MW=25.3;
gam=1.24;

% Equations------------------------------------

% Independent
pen=pi*epsilon/N;
H=Rp*sin(pen);

% Calculate Ri for neutral Phase I burn
syms Rii
TH_2i=atan(H*tan(pen)/(H-Rii*tan(pen)));
coef=pi/2-TH_2i+pi/N-cot(TH_2i)==0;
Ri=double(solve(coef,Rii));

TH_2=atan(H*tan(pen)/(H-Ri*tan(pen)));
beta=(pi/2-TH_2+pen);
yo=H/cos(TH_2);

% Calculate Throat Area
S1i=H/sin(TH_2)-(0+f)*cot(TH_2);
S2i=(0+f)*beta;
S3i=(Rp+0+f)*(pi/N-pen);
Si=2*N*(S1i+S2i+S3i);
Abi=Si*L;
At=Abi*a*rho*cstar/32.2/po^(1-n);


% Calculate pressure change throughout burn

% Phase I
y=Web/100;
while y+f<=yo
    S1=H/sin(TH_2)-(y+f)*cot(TH_2);
    S2=(y+f)*beta;
    S3=(Rp+y+f)*(pi/N-pen);
    S=2*N*(S1+S2+S3);
    Ab=S*L;
    pc1=((Ab/At)*a*rho*cstar/32.2)^(1/(1-n));
    if y==Web/100
        pc=pc1;
    else
        pc=[pc,pc1];
    end
    y=y+Web/100;
end

web1=Rp*sin(pen)/cos(TH_2)-f;
r1=a*po^n;
t1=web1/r1;
tb=t1;
web2=Web-web1;

% Phase II
while Web+f>y+f>yo
    psi=atan(sqrt((y+f)^2-H^2)/H)-TH_2;
    S2=(y+f)*(beta-psi');
    S3=(Rp+y+f)*(pi/N-pen);
    S=2*N*(S2+S3);
    Ab=S*L;
    pc2=((Ab/At)*a*rho*cstar/32.2)^(1/(1-n));
    pc=[pc,pc2];
    r2=a*pc2^n;
    t2i=.03/r2;
    tb=tb+t2i;
    y=y+Web/100;
end

% Check loop outputs
size(pc);
pc(:,1:10);
range=[min(pc),max(pc)];


% Plot Chamber Pressure vs. Time
t=linspace(0,tb);
plot(t,pc)
hold on
xlim([0,7])
ylim([0,2000])
xlabel('Burn Time')
ylabel('Chamber Pressure (psi)')
title('Neutral Star (Phase I) Chamber Pressure vs. Time')

plot(t1,pc)
text(t1-.5,1700,'End of Phase I')
plot(tb,pc)
text(tb-.5,700,'End of Burn')

