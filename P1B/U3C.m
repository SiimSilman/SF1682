clc
clearvars
clf
%% Initialdata
% Konstanter
m1 = 465; 
m2 = 55;
k1 = 641; 
k2 = 185118; 
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
[t,Vode45] = ode45( @(t,V) quartercar(t, V, m1, m2, k1, k2, c1, c2, H, L, v), [0 Tlap], V, options);

%% Plott
hold on
plot(t, Vode45(:,1)) % z1
plot(t, Vode45(:,2)) % z2
legend('z1', 'z2')
title('ode45, RelTol = 1e-6, optimerad k_1 och k_2')
hold off

%% Kommentar
z1Max = max(Vode45(:,1));
z2Max = max(Vode45(:,2));
fprintf('z1 max utslag %.4fm \n', z1Max)
fprintf('z2 max utslag %.4fm \n', z2Max)

% För bättrad resultat, vi ser jämförelsevis att utslaget vid hjulen blir
%   något större, 0.3m jämfört med 0.25m, dock blir utslaget vid chassit
%   något mindre, 0.008m jämfört med 0.025m.
% Det vill säga att körsäkerheten blir faktor 1/5 sämre (om linjärt samtband) 
%   men komforten blir faktor 3 bättre. 