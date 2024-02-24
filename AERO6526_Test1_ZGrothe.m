% Zane Grothe
% AERO 6526
% Test 1
% 3/12/21

% Problem 31

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
po=600; %psia
Web=3; %in

% Propellant Parameters
a=0.00475; %in/sec/psin
n=0.6;
cstar=3803; %ft/sec
rho=0.063; %lbm/in^3
Tc=2759; %Rankine
MW=21.8;
gam=1.26;

% Equations------------------------------------

% Independent
pen=pi*epsilon/N;
H=Rp*sin(pen);

% Calculate Ri for neutral Phase I burn
syms Rii
TH_2i=atan(H*tan(pen)/(H-Rii*tan(pen)));
coef=pi/2-TH_2i+pi/N-cot(TH_2i)==0;
Ri=double(solve(coef,Rii))

TH_2=atan(H*tan(pen)/(H-Ri*tan(pen)));
beta=(pi/2-TH_2+pen);
yo=H/cos(TH_2);

% Calculate Throat Area
S1i=H/sin(TH_2)-(0+f)*cot(TH_2);
S2i=(0+f)*beta;
S3i=(Rp+0+f)*(pi/N-pen);
Si=2*N*(S1i+S2i+S3i);
Abi=Si*L;
At=Abi*a*rho*cstar/32.2/po^(1-n)

