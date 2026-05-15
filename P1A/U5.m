clc    
clearvars  
close all
%% U4 
% Utifrån givna initialvärden från U3C, relativ tolerans 1e-8:
% - approximera körbanan med inbyggda ode45 lösare,
%   jämför med analytisk lösning.

%% Initial värden
b = 1;  aR = 1.2; aL = 0.4; wR = 3; wL = 1;

% Initial villkor
theta0 = 0; x0 = 0; y0 = 0;

% Hämta startkordinater, radien och tiden för ett varv
[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0);
Tlap_ny = 10;

s = [xbar, ybar, theta0]';

%% ODE45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Refine', 1);
[tODE,s] = ode45( @(tODE,s) fvel(tODE,s,b,aL,aR,wL,wR)', [0 Tlap_ny], s, options);

%% Plot
hold on
plot(s(1,1), s(1,2), 'x')
plot(s(:,1), s(:,2))
plot(s(end,1), s(end,2), 'x')
hold off

legend('startpkt', 'ode45', 'slutpkt')
axis equal
axis image

