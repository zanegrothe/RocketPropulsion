% Zane Grothe
% AERO 6526
% HW 4
% 3/11/21

% Problem #3 (star grain)

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

% Calculations------------------------------------

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
tb1=web1/r1;
tb=tb1;

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
    tb2=Web/100/r2;
    tb=tb+tb2;
    y=y+Web/100;
end

% Check loop outputs
pc(:,1:10);
range=[min(pc),max(pc)];

% Calculate Ignition Transient
GAM=sqrt(gam)*((gam+1)/2)^-((gam+1)/2/(gam-1));
V=2443.32; %in^3
S11=H/sin(TH_2)-(f)*cot(TH_2);
S22=(f)*beta;
S33=(Rp+f)*(pi/N-pen);
SS=2*N*(S11+S22+S33);
Abb=SS*L; %in^2
K=Abb/At;

ti=(1/(1-n))*(V/GAM^2/cstar/At/12)*log((rho*a*cstar*K*min(pc)^(1-n))/(rho*a*cstar*K-min(pc)^(1-n)));

tki=ti/25;
while tki<ti
    syms pcii
    tit=(1/(1-n))*(V/GAM^2/cstar/At/12)*log((rho*a*cstar*K*pcii^(1-n))/(rho*a*cstar*K-pcii^(1-n)))==tki;
    pcig=double(solve(tit,pcii));
    if tki==ti/25
        pci=pcig;
    else
        pci=[pci,pcig];
    end
    tki=tki+ti/25;
end

% Calculate Ignition Tailoff
tt=0;
pcto=max(pc);
while pcto>(min(pc)*.1)
    pcto=max(pc)*exp(-(GAM^2*At*cstar*(tt)*12/V));
    if pcto==max(pc)
        pct=pcto;
    else
        pct=[pct,pcto];
    end
    tt=tt+.01;
end

%pc=[pci,pc,pct];
[pcir,pcic]=size(pci);
[pcr,pcc]=size(pc);
[pctr,pctc]=size(pct);
tot=ti+tb+tt;

% Plot Chamber Pressure vs. Time
t1=linspace(0,ti,pcic);
plot(t1,pci)
hold on
t2=linspace(ti,tb+ti);
plot(t2,pc)
t3=linspace(tb+ti,tot,pctc);
plot(t3,pct)

xlim([0,tot+tot/5])
ylim([0,max(pc)+max(pc)/5])
xlabel('Burn Time')
ylabel('Chamber Pressure (psi)')
title('Chamber Pressure vs. Time')

plot(ti+tb1,pc)
text(ti+(tb1-tb1/10),max(pc)+max(pc)/10,'End of Phase I')


