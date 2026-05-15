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
N = 1000;
h = Tlap/N;

dt1 = 5e-3; 
dt2 = 5e-4;

tid = (0:h:Tlap)';
tid1 = (0:dt1:Tlap)';
tid2 = (0:dt2:Tlap)';

N1 = length(tid1);
N2 = length(tid2);


%% ODE45, Relativ tolerans 10^-6
V = V0;
options = odeset('RelTol', 1e-6, 'Refine', 1);
[t,Vode45] = ode45( @(t,V) quartercar(t, V, m1, m2, k1ref, k2ref, c1, c2, H, L, v), [0 Tlap], V, options);

%% Eulers metod, dt = 5e-3
z2EF1 = zeros(1,N1)';
V = [z1, z2, dz1, dz2]';
for i = 1:N1
    dV = quartercar(tid1(i), V, m1, m2, k1ref, k2ref, c1, c2, H, L, v);
    V = V + dt1 * dV;
    z2EF1(i) = V(2);
end

%% Eulers metod, dt = 5e-4
z2EF2 = zeros(1,N2)';
V = [z1, z2, dz1, dz2]';
for i = 1:N2
    dV = quartercar(tid2(i), V, m1, m2, k1ref, k2ref, c1, c2, H, L, v);
    V = V + dt2 * dV;
    z2EF2(i) = V(2);
end

%% Plottar
hold on
plot(t, Vode45(:,2)) % z2
plot(tid1, z2EF1)
plot(tid2, z2EF2)
legend('ode45', 'EF: 5e-3', 'EF: 5e-4')
title('z2 för olika metoder')
hold off

%% Kommentar
% Ode45 och EF med dt = 5e-4 ger väldigt liknande lösningar vilket är
%   förväntad då ode45 ska vara väldigt stark metod samt att EF blir 
%   noggrannare med fler tidssteg.
% EF med dt = 5e-3 är något sämre approximation om man jämför den med de
%   andra, då får man överväga ifall man vill ha sämre approximation men
%   spara datorkraft eller vice versa.

