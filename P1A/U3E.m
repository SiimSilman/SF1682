clc    
clearvars  
close all
%% U3 E
% Utifrån givna initialvärden och Tlap_ny = 7pi/8:
% - Plotta analytisk lösning
% - approximera körbanan med Euler-Fram med 100 tidsteg
% - approximera körbanan med Runge-Kutta 4 med 100 tidsteg
% - Felet i x-kord för tlap jämfört med analytisk lösning?

% Slutsats:
% Bekräftar tidigare påståendet om noggrannhet av RK4

%% Initial värden
b = 1;  aL = 0; aR = 0; wL = 4; wR = 2;

% Initial villkor
theta0 = 0; x0 = 0; y0 = 0;

% Hämta startkordinater, radien och tiden för ett varv
[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0);

s0 = [xbar, ybar, theta0];

%% Fördela tiden med X tidsteg
Tlap_ny = 7*pi/8;

N = 100;
h = Tlap_ny/N;
tid = 0:h:Tlap_ny;

%% Allokera minne för sparade värden
% Euler Fram, 100 steg
xEF = zeros(N+1, 1);
yEF = zeros(N+1, 1);
thetaEF = zeros(N+1, 1);

% Runge-Kutta 4, 100 steg
xRK = zeros(N+1, 1);
yRK = zeros(N+1, 1);
thetaRK = zeros(N+1, 1);

%% Euler Fram med 100 steg
s = s0;
for i = 1:N+1
    xEF(i) = s(1);
    yEF(i) = s(2);
    thetaEF(i) = s(3);
    
    ds = fvel(tid(i), s, b, aL, aR, wL, wR);
    s = s + ds*h;
end

%% Runge Kutta 4 med 100 steg
s = s0;
for i = 1:N+1
    xRK(i) = s(1);
    yRK(i) = s(2);
    thetaRK(i) = s(3);

    k1 = fvel(tid(i), s, b, aL, aR, wL, wR);
    k2 = fvel(tid(i) + h/2, s + (h/2) * k1, b, aL, aR, wL, wR);
    k3 = fvel(tid(i) + h/2, s + (h/2) * k2, b, aL, aR, wL, wR);
    k4 = fvel(tid(i) + h, s + h * k3, b, aL, aR, wL, wR);
    
    ds = (k1 + 2*k2 + 2*k3 + k4) / 6;
    s = s + ds*h;
end

%% Analytisk lösning
D = (wR - wL) / b;
for i = 1:N+1
    thetaA = D*tid + theta0;
    xA = x0 + R * ( sin(thetaA) );
    yA = y0 - R * ( cos(thetaA) );
end

%% Felet i x-led vid t = 7*pi/8 jämfört med analytisk lösning
felEF = abs(xA(end) - xEF(end));
felRK = abs(xA(end) - xRK(end));

fprintf('Felet av respektive metod jämfört med analytisk lösning: \n')
fprintf('Felet för EF: %.4f \n', felEF)
fprintf('Felet för RK: %.9f \n', felRK)

%% Plottar
hold on
plot(0, 0, 'o')         % Mittpunkt
plot(xEF, yEF, 'g')     % EF 100 steg
plot(xRK, yRK, 'r')     % RK4 100 steg
plot(xA, yA, 'c')       % Analytisk 100 steg

legend('Mittp', 'EF100', 'RK4', 'Analytisk')

axis equal
axis image
hold off