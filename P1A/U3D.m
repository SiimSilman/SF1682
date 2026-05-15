clc    
clearvars  
close all
%% U3 
% Utifrån givna initialvärden:
% - approximera körbanan med Euler-Fram med 100 tidsteg
% - approximera körbanan med Euler-Fram med 1000 tidsteg
% - approximera körbanan med Runge-Kutta 4 med 100 tidsteg

% Slutsats:
% Vet att RK4 är mycket noggrann metod, med ökade tidsteg för Euler-Fram
% närmar sig approximationen RK4.

%% Initial värden
b = 1;  aL = 0; aR = 0; wL = 4; wR = 2;

% Initial villkor
theta0 = 0; x0 = 0; y0 = 0;

% Hämta startkordinater, radien och tiden för ett varv
[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0);

s0 = [xbar, ybar, theta0];

%% Fördela tiden med X tidsteg
N1 = 100;
h1 = Tlap/N1;
tid1 = 0:h1:Tlap;

N2 = 1000;
h2 = Tlap/N2;
tid2 = 0:h2:Tlap;

%% Allokera minne för sparade värden
% Euler Fram, 100 steg
dx1 = zeros(N1+1, 1);
dy1 = zeros(N1+1, 1);
dtheta1 = zeros(N1+1, 1);

% Euler Fram, 1000 steg
dx2 = zeros(N2+1, 1);
dy2 = zeros(N2+1, 1);
dtheta2 = zeros(N2+1, 1);

% Runge-Kutta 4, 100 steg
dx3 = zeros(N1+1, 1);
dy3 = zeros(N1+1, 1);
dtheta3 = zeros(N1+1, 1);

%% Euler Fram med 100 steg
s = s0;
for i = 1:N1+1
    dx1(i) = s(1);
    dy1(i) = s(2);
    dtheta1(i) = s(3);
    
    ds = fvel(tid1(i), s, b, aL, aR, wL, wR);
    s = s + ds*h1;
end

%% Euler Fram med 1000 steg
s = s0;
for i = 1:N2+1
    dx2(i) = s(1);
    dy2(i) = s(2);
    dtheta2(i) = s(3);
    
    ds = fvel(tid2(i), s, b, aL, aR, wL, wR);
    s = s + ds*h2;
end

%% Runge Kutta 4 med 100 steg
s = s0;
for i = 1:N1+1
    dx3(i) = s(1);
    dy3(i) = s(2);
    dtheta3(i) = s(3);

    k1 = fvel(tid1(i), s, b, aL, aR, wL, wR);
    k2 = fvel(tid1(i) + h1/2, s + (h1/2) * k1, b, aL, aR, wL, wR);
    k3 = fvel(tid1(i) + h1/2, s + (h1/2) * k2, b, aL, aR, wL, wR);
    k4 = fvel(tid1(i) + h1, s + h1 * k3, b, aL, aR, wL, wR);
    
    ds = (k1 + 2*k2 + 2*k3 + k4) / 6;
    s = s + ds*h1;
end
%% Plottar
hold on
plot(0, 0, 'o')         % Mittpunkt
plot(dx1, dy1, 'b')     % EF 100 steg
plot(dx2, dy2, 'g')     % EF 1000 steg
plot(dx3, dy3, 'r')     % RK4 100 steg

legend('Mittp', 'EF100', 'EF1000', 'RK4')

axis equal
axis image
hold off