clc    
clearvars  
close all
%% U4 
% Utifrån givna initialvärden från U3C, relativ tolerans 1e-8:
% - approximera körbanan med inbyggda ode45 lösare,
%   jämför med analytisk lösning.

%% Initial värden
b = 1;  aL = 0; aR = 0; wL = 4; wR = 2;

% Initial villkor
theta0 = 0; x0 = 0; y0 = 0;

% Hämta startkordinater, radien och tiden för ett varv
[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0);

s = [xbar, ybar, theta0]';

%% Fördela tiden med 100 tidsteg
N = 100;
h = Tlap/N;
tid = (0:h:Tlap)';

%% ODE45
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Refine', 1);
[tODE,s] = ode45( @(tODE,s) fvel(tODE,s,b,aL,aR,wL,wR)', [0 Tlap], s, options);

%% Analytisk lösning
D = (wR - wL) / b;
thetaA = D*tODE + theta0;
xA = x0 + R * (sin(thetaA));
yA = y0 - R * (cos(thetaA));

%% Fel i x och y som funktion av tiden
felX = abs(xA - s(:,1));
felY = abs(yA - s(:,2));

%% Plotta
fprintf('ode45 tar %i tidsteg \n', length(tODE))

subplot(2,1,1)
    plot(tODE, felX, 'r')
    legend('fel i X')

subplot(2,1,2)
    plot(tODE, felY, 'g')
    legend('fel i Y')

xlabel('t = t_{ode45} d.v.s 37 tidsteg')
ylabel('Absolutbeloppet av felet')





