% Zane Grothe
% AERO 6526
% HW 3
% 3/1/21

% Problem #1

clear all
close all
clc


% Pressure data in atm
ANp=[68,69,71,102,103,115,140,141,180,181,205,206,235,236,280,282,320,325];
RDXp=[1,6.5,13,20,21,22,30,31,72,74,107,120,130,160,200,250,251,360,460,700,760];
ADNp=[14,20,27,34,40,41,55,60,70,80,100,120,150,160];
APp=[27,28,35,40,41,58,59,70,82,85,101,102,110,120,130,140];
CL20p=[7,14,27,55,70,100];
GAPp=[10,20,30,40,50,60];

% Burn rate data in cm/sec
ANr=[.18,.21,.20,.25,.29,.32,.38,.40,.45,.52,.49,.54,.63,.68,.80,.90,.95,1.1];
RDXr=[.039,.15,.32,.39,.40,.42,.64,.69,1.2,1.5,1.6,1.9,2.4,2.9,3.0,3.5,3.9,4.9,7.0,7.9,10.0];
ADNr=[2.6,2.5,3.4,2.4,3.0,2.8,3.3,3.1,3.5,3.2,3.6,4.0,4.5,4.8];
APr=[.39,.39,.5,.6,.58,.71,.69,.80,.98,.91,1.0,1.1,1.1,1.2,1.3,1.4];
CL20r=[.46,.63,1.1,2.0,2.4,3.0];
GAPr=[.51,.71,.85,.95,1.1,1.3];

% Find coefficients for a linear polynomial and convert to n&a
lc_AN=polyfit(log(ANp),log(ANr),1);
ANn=lc_AN(1,1);
ANa=exp(lc_AN(1,2));

lc_RDX=polyfit(log(RDXp),log(RDXr),1);
RDXn=lc_RDX(1,1);
RDXa=exp(lc_RDX(1,2));

lc_ADN=polyfit(log(ADNp),log(ADNr),1);
ADNn=lc_ADN(1,1);
ADNa=exp(lc_ADN(1,2));

lc_AP=polyfit(log(APp),log(APr),1);
APn=lc_AP(1,1);
APa=exp(lc_AP(1,2));

lc_CL20=polyfit(log(CL20p),log(CL20r),1);
CL20n=lc_CL20(1,1);
CL20a=exp(lc_CL20(1,2));

lc_GAP=polyfit(log(GAPp),log(GAPr),1);
GAPn=lc_GAP(1,1);
GAPa=exp(lc_GAP(1,2));

n_a=[ANn,ANa;RDXn,RDXa;ADNn,ADNa;APn,APa;CL20n,CL20a;GAPn,GAPa];
disp('      n         a')
disp(n_a)

% AN
% RDX
% ADN >10atm
% AP <150atm
% CL20
% GAP
