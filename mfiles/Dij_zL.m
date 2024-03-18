

% read and analyze Dij for different z/L runs for alfalfa.
% Boulding Alfalfa

clear all
close all

Dij_00=csvread('D:\CanAlfalfa\Dij_Alfalfa.csv');
Dij_25=csvread('D:\CanAlfalfa\Dij_zlm025.csv');
Dij_50=csvread('D:\CanAlfalfa\Dij_zlm050.csv');
Dij_100=csvread('D:\CanAlfalfa\Dij_zlm100.csv');
Dij_300=csvread('D:\CanAlfalfa\Dij_zlm300.csv');
Dij_p50=csvread('D:\CanAlfalfa\Dij_zlp050.csv');
Dij_p100=csvread('D:\CanAlfalfa\Dij_zlp100.csv');

zL=[0,-0.25,-0.5,-1,-3,0.5,1.];

figure(2)
clf
plot(Dij_00,Dij_25,'.')
hold on

plot(Dij_00,Dij_50,'.')
hold on
plot(Dij_00,Dij_100,'.')

hold on
plot(Dij_00,Dij_300,'.')
hold on
plot(Dij_00,Dij_p50,'.')
hold on
plot(Dij_00,Dij_p100,'.')



p1=polyfit(Dij_00,Dij_25,1);

p2=polyfit(Dij_00,Dij_50,1);

p3=polyfit(Dij_00,Dij_100,1);

p4=polyfit(Dij_00,Dij_300,1);

p5=polyfit(Dij_00,Dij_p50,1);

p6=polyfit(Dij_00,Dij_p100,1);

y=[1,p1(1),p2(1),p3(1),p4(1),p5(1),p6(1)]

figure(1)
clf;

plot(zL(1:5),y(1:5), '.-')
hold on
plot(zL(6:7),y(6:7), '+-')
xlabel('z/L')
ylabel('Dij(0)/Dij(z/L)')
title('Bouldin Alfalfa, h=1 m')


% non linear fit

Y=y(1:6);
X=zL(1:6);

%modelfun = @(b,x)(b(1)+b(2)*exp(-b(3)*x));
modelfun= @(b,x) (b(1) .*b(2)) ./(b(2)+X);
beta0=[0.1 0.1];

beta = nlinfit(X,Y,modelfun,beta0);
