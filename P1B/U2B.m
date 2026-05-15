clc
clearvars
clf
%% Initialdata
% Konstanter
m1 = 465; 
m2 = 55;
k1ref = 5350; 
k2ref = 136100; 
c1 = 310; 
c2 = 1250; 
v = 63/3.6; 
H = 0.27; 
L = 1.1; 

% Initialvärden
z1 = 0; 
z2 = 0; 
dz1 = 0; 
dz2 = 0; 
V0 = [z1, z2, dz1, dz2]';

% Tid
Tlap = L/v + 1;

%% ODE45, Relativ tolerans 10^-6
V = V0;
options = odeset('RelTol', 1e-6, 'Refine', 1);
[t,Vode45] = ode45( @(t,V) quartercar(t, V, m1, m2, k1ref, k2ref, c1, c2, H, L, v), [0 Tlap], V, options);

%% Plottar
subplot(2,1,1)
hold on
plot(t, Vode45(:,1)) % z1
plot(t, Vode45(:,2)) % z2
legend('z1', 'z2')
title('ode45, RelTol = 1e-6')
hold off

subplot(2,1,2)
plot(t, zeros(size(t)),'-o')
xlim([0 Tlap])
title('ode45 steglängd')

%% Kommentar
% Ode45 är en adaptiv lösare som tar korta tidssteg när dervitatan är stor, 
%   och stora steg då derivatan är liten.
% Det syns i graferna där i början är avvikelsen stor och många korta tidssteg 
%   tas och i slutet är avvikelsen liten och få långa tidssteg tas. 
% Ytterligare observation är att skulle enbart z1 mätas med ode skulle
%   väldigt få tidssteg tas då den är väldigt stabil över mätnings perioden.