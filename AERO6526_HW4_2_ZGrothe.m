% Zane Grothe
% AERO 6526
% HW 4
% 3/11/21

% Problem #2 (star grain)

clear all
close all
clc

% Givens------------------------------------

% Design Parameters
Ro=5; %in
Rp=2.6; %in
Ri=1.2; %in
f=0.2; %in
epsilon=0.8;
N=5;
L=24; %in
TH_2=60*pi/180; %radians
Dt=2.0; %in
pa=14.7; %psi
pe=pa; %psi
AeAt=8.0; 
Web=5; %in

% Propellant Parameters (HTPB/AP/AL)
a=0.025; %in/sec/psin
n=0.35;
cstar=5364; %ft/sec
rho=0.064; %lbm/in^3
Tc=5700; %Rankine
MW=24.7;
gam=1.24;

% Equations------------------------------------

% Independent
At=pi*(Dt^2)/4;
pen=pi*epsilon/N;
H=Rp*sin(pen);
beta=(pi/2-TH_2+pen);
yo=H/cos(TH_2);

syms M
AA=(1/M)*[2/(gam+1)*((1+(gam-1)/2*M^2))]^((gam+1)/2/(gam-1))==AeAt;
Me=double(solve(AA,M)); % (exit Mach number)
po=pe*[(1+(gam-1)/2*Me^2)^(gam/(gam-1))]; %psi (chamber pressure)

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

% Phase II
while Web>y
    psi=atan(sqrt((y+f)^2-H^2)/H)-TH_2;
    S2=(y+f)*(beta-psi');
    S3=(Rp+y+f)*(pi/N-pen);
    S=2*N*(S2+S3);
    Ab=S*L;
    pc2=((Ab/At)*a*rho*cstar/32.2)^(1/(1-n));
    pc=[pc,pc2];
    r2=a*pc2^n;
    t2i=Web/100/r2;
    tb=tb+t2i;
    y=y+Web/100;
end

% Check loop outputs
[b,c]=size(pc);
pc(:,1:10);
range=[min(pc),max(pc)];

% Maximum chamber pressure
MaximumChamberPressure=max(pc)

% Phase I Port Area ---------
Sp=H/sin(TH_2)-(f)*cot(TH_2);
Ap1=1/2*H*(Rp*cos(pen)+H*tan(TH_2))-1/2*Sp^2*tan(TH_2);
Ap2=1/2*(f)^2*beta;
Ap3=1/2*(Rp+f)^2*(pi/N-pen);
Ap=2*N*(Ap1+Ap2+Ap3);
PortArea=Ap
