% Zane Grothe
% AERO 6526
% HW 4
% 3/11/21

% Problem #1

clear all
close all
clc

% Calculating Burn Area and Port Area for Star Grain

% Givens
Ri=4;
Rp=7.5;
epsilon=.9;
N=5;
f=1;
L=60;

% Preliminary Calculations ----------

pen=pi*epsilon/N;
H=Rp*sin(pen);
TH_2=atan(H*tan(pen)/(H-Ri*tan(pen)));
beta=(pi/2-TH_2+pen);

% Phase I ----------

S1=H/sin(TH_2)-(f)*cot(TH_2);
S2=(f)*beta;
S3=(Rp+f)*(pi/N-pen);
S=2*N*(S1+S2+S3);
Ab=S*L

Ap1=1/2*H*(Rp*cos(pen)+H*tan(TH_2))-1/2*S1^2*tan(TH_2);
Ap2=1/2*(f)^2*beta;
Ap3=1/2*(Rp+f)^2*(pi/N-pen);
A=2*N*(Ap1+Ap2+Ap3)

