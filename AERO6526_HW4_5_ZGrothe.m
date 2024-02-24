% Zane Grothe
% AERO 6526
% HW 4
% 3/11/21

% Problem #5

clear all
close all
clc

% Givens----------

Ru=1545.43; %ft lbf/lbm mole R

% Propellant (PBAA/AP/AN)
a=.0285; %in/sec/psin
n=.35;
rho=.064; %lbm/in^3
T=5700; %Fahrenheit
To=T+460; %Rankine
gam=1.24;
cstar=5364; %ft/s
MW=24.7;

% Denison/Baum model
A=12;
B=.7;

% Calculations----------

omega=0;
k=1;
while omega<pi
    Omega=omega/10;
    lam=1+sqrt(1+4*1i*Omega);
    R=(A*B*n)/(lam+(A/lam)-(1+A)+(A*B));
    theta=angle(R);
    Lstar=(R*sin(theta))/(cstar*MW/Ru/32.2/To)/omega;
    alpha=(R*cos(theta)-1)/(cstar*MW/Ru/32.2/To)/Lstar;
    if alpha>0
        if k==1
            disp('Motor is unstable below Lstar=')
            disp(real(Lstar/12))
            k=k+1;
        end
    end
    if omega==0
        Rnm=R/n;
        thetam=theta;
        Lstarm=Lstar;
        alpham=alpha;
    else
        Rnm=[Rnm,R/n];
        thetam=[thetam,theta];
        Lstarm=[Lstarm,Lstar];
        alpham=[alpham,alpha];
    end
    omega=omega+.01;
end

rangeR=[min(Rnm),max(Rnm)];
rangeT=[min(thetam),max(thetam)];

[Rr,Rc]=size(Rnm);
[Thr,Thc]=size(thetam);

O=linspace(0,Omega,Rc);
Rnm=real(Rnm);
thetam=real(thetam);
Lstarm=real(Lstarm);
alpham=real(alpham);
Lstar=real(Lstar);

% Plots----------

figure(1)
plot(O,Rnm)
xlabel('Omega (radians)')
ylabel('Theta                            R/n              ')
hold on
plot(O,thetam)
hold off

figure(2)
plot(O,Rnm)
xlabel('Omega (radians)')
ylabel('R/n')

figure(3)
plot(O,thetam)
xlabel('Omega (radians)')
ylabel('Theta')

Lstarm=Lstarm/12;
Lstar=Lstar/12;
alpham=alpham*180/pi;

figure(4)
plot(Lstarm,alpham)
xlabel('L* (ft)')
ylabel('Alpha (degrees)')
hold on

y_zero=linspace(Lstar,Lstarm(1,2));
plot(y_zero,0)


