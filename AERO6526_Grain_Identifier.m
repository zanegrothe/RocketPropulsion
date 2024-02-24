% Zane Grothe
% 3/2021

% SRM Grain Geometry Identification

% Based on givens, this script will determine whether your grain is a star,
% a short spoke wagon wheel, or a long spoke wagon wheel.

clear all
close all
clc

% Givens----------

Rp=;
Ri=;
ep=; %epsilon
N=;
th_2=; %theta/2 (radians)

% Evaluate----------

pen=pi*ep/N;
h=(Rp*cos(pen)-Rp*sin(pen)/tan(th_2)-Ri)*sin(th_2); %star or wagon

if h<0
    disp('h is negative indicating you have a star grain')
else
    l_s=Rp*cos(pen)-Ri-Rp*sin(pen)/sin(th_2); %long or short spoke
    if l_s<0
        disp('h is positive but l_s is negative indicating you have a short spoke wagon wheel')
    else
        disp('h is positive and l_s is positive indicating you have a long spoke wagon wheel')
    end
end

