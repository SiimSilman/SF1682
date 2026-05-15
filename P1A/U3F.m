clc    
clearvars  
close all
%% U3 F
% Utifrån givna initialvärden och Tlap_ny = 7pi/8:
% - approximera körbanan med Euler-Fram och Runge-Kutta 4 med
%   100 tidsteg dt, dt/2 och dt/4. 
% - Vad blir felet i x-kord vid t = Tlap för de två metoderna?

% Slutsats:
% EF fel minskar med faktor 2 vid halverad tidsteg.
% RK fel minskas med ca faktor 16 vid halverad tidsteg.

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

h1 = Tlap_ny/N;
h2 = h1/2;
h4 = h1/4;

tid1 = 0:h1:Tlap_ny;
tid2 = 0:h2:Tlap_ny;
tid4 = 0:h4:Tlap_ny;

%% Allokera minne för sparade värden
% Euler Fram, 100 steg
xEF1 = zeros(N+1, 1);
yEF1 = zeros(N+1, 1);
thetaEF1 = zeros(N+1, 1);

% Euler Fram, 200 steg
xEF2 = zeros(2*N+1, 1);
yEF2 = zeros(2*N+1, 1);
thetaEF2 = zeros(2*N+1, 1);

% Euler Fram, 400 steg
xEF4 = zeros(4*N+1, 1);
yEF4 = zeros(4*N+1, 1);
thetaEF4 = zeros(4*N+1, 1);

% Runge-Kutta 4, 100 steg
xRK1 = zeros(N+1, 1);
yRK1 = zeros(N+1, 1);
thetaRK1 = zeros(N+1, 1);

% Runge-Kutta 4, 200 steg
xRK2 = zeros(2*N+1, 1);
yRK2 = zeros(2*N+1, 1);
thetaRK2 = zeros(2*N+1, 1);

% Runge-Kutta 4, 400 steg
xRK4 = zeros(4*N+1, 1);
yRK4 = zeros(4*N+1, 1);
thetaRK4 = zeros(4*N+1, 1);

%% 100 steg
sEF = s0;
sRK = s0;
for i = 1:N+1
    % Euler Fram
    xEF1(i) = sEF(1);
    yEF1(i) = sEF(2);
    thetaEF1(i) = sEF(3);
    
    dsEF = fvel(tid1(i), sEF, b, aL, aR, wL, wR);
    sEF = sEF + dsEF * h1;

    % Runge Kutta
    xRK1(i) = sRK(1);
    yRK1(i) = sRK(2);
    thetaRK1(i) = sRK(3);

    k1 = fvel(tid1(i), sRK, b, aL, aR, wL, wR);
    k2 = fvel(tid1(i) + h1/2, sRK + (h1/2) * k1, b, aL, aR, wL, wR);
    k3 = fvel(tid1(i) + h1/2, sRK + (h1/2) * k2, b, aL, aR, wL, wR);
    k4 = fvel(tid1(i) + h1, sRK + h1 * k3, b, aL, aR, wL, wR);
    
    dsRK = (k1 + 2*k2 + 2*k3 + k4) / 6;
    sRK = sRK + dsRK * h1;
end

%% 200 steg
sEF = s0;
sRK = s0;
for i = 1:2*N+1
    % Euler Fram
    xEF2(i) = sEF(1);
    yEF2(i) = sEF(2);
    thetaEF2(i) = sEF(3);
    
    dsEF = fvel(tid2(i), sEF, b, aL, aR, wL, wR);
    sEF = sEF + dsEF * h2;

    % Runge Kutta
    xRK2(i) = sRK(1);
    yRK2(i) = sRK(2);
    thetaRK2(i) = sRK(3);

    k1 = fvel(tid2(i), sRK, b, aL, aR, wL, wR);
    k2 = fvel(tid2(i) + h2/2, sRK + (h2/2) * k1, b, aL, aR, wL, wR);
    k3 = fvel(tid2(i) + h2/2, sRK + (h2/2) * k2, b, aL, aR, wL, wR);
    k4 = fvel(tid2(i) + h2, sRK + h2 * k3, b, aL, aR, wL, wR);
    
    dsRK = (k1 + 2*k2 + 2*k3 + k4) / 6;
    sRK = sRK + dsRK * h2;
end

%% 400 steg
sEF = s0;
sRK = s0;
for i = 1:4*N+1
    % Euler Fram
    xEF4(i) = sEF(1);
    yEF4(i) = sEF(2);
    thetaEF4(i) = sEF(3);
    
    dsEF = fvel(tid4(i), sEF, b, aL, aR, wL, wR);
    sEF = sEF + dsEF * h4;

    % Runge Kutta
    xRK4(i) = sRK(1);
    yRK4(i) = sRK(2);
    thetaRK4(i) = sRK(3);

    k1 = fvel(tid4(i), sRK, b, aL, aR, wL, wR);
    k2 = fvel(tid4(i) + h4/2, sRK + (h4/2) * k1, b, aL, aR, wL, wR);
    k3 = fvel(tid4(i) + h4/2, sRK + (h4/2) * k2, b, aL, aR, wL, wR);
    k4 = fvel(tid4(i) + h4, sRK + h4 * k3, b, aL, aR, wL, wR);
    
    dsRK = (k1 + 2*k2 + 2*k3 + k4) / 6;
    sRK = sRK + dsRK * h4;
end

%% Analytisk lösning
D = (wR - wL) / b;
for i = 1:4*N+1
    thetaA = D*tid1 + theta0;
    xA = x0 + R * ( sin(thetaA) );
    yA = y0 - R * ( cos(thetaA) );
end

%% Noggrannhetsordning
% Felet, 100 steg
felEF1 = abs(xA(end) - xEF1(end));
felRK1 = abs(xA(end) - xRK1(end));

% Felet, 200 steg
felEF2 = abs(xA(end) - xEF2(end));
felRK2 = abs(xA(end) - xRK2(end));

% Felet, 400 steg
felEF4 = abs(xA(end) - xEF4(end));
felRK4 = abs(xA(end) - xRK4(end));

% Beräkna kvoten av felen
kvotEF12 = felEF1 / felEF2;
kvotEF24 = felEF2 / felEF4;

kvotRK12 = felRK1 / felRK2;
kvotRK24 = felRK2 / felRK4;

% Beräkna noggrannhetsordningen
pEF12 = log(felEF1 / felEF2) / log(2);
pEF24 = log(felEF2 / felEF4) / log(2);

pRK12 = log(felRK1 / felRK2) / log(2);
pRK24 = log(felRK2 / felRK4) / log(2);


%% Plottar
fprintf('Felet och noggrannhetsordningen:\n\n')

fprintf('Felet för EF med 100 steg: %.9f \n', felEF1)
fprintf('Felet för EF med 200 steg: %.9f \n', felEF2)
fprintf('Felet för EF med 400 steg: %.9f \n\n', felEF4)

fprintf('Felet för RK med 100 steg: %.12f \n', felRK1)
fprintf('Felet för RK med 200 steg: %.12f \n', felRK2)
fprintf('Felet för RK med 400 steg: %.12f \n\n', felRK4)

fprintf('Konvergens för EF med 100 resp 200 steg: %.9f \n', pEF12)
fprintf('Konvergens för EF med 200 resp 400 steg: %.9f \n\n', pEF24)

fprintf('Konvergens för RK med 100 resp 200 steg: %.9f \n', pRK12)
fprintf('Konvergens för RK med 200 resp 400 steg: %.9f \n\n', pRK24)

hold on
plot(xEF1, yEF1)        % EF 100 steg
plot(xRK1, yRK1)        % RK4 100 steg

plot(xEF2, yEF2)        % EF 200 steg
plot(xRK2, yRK2)        % RK4 200 steg

plot(xEF4, yEF4)        % EF 400 steg
plot(xRK4, yRK4)        % RK4 400 steg

plot(xA, yA)            % Analytisk 400 steg
hold off

legend('EF100', 'RK100', 'EF200', 'RK200', 'EF400', 'RK400', 'Analytisk')
axis equal
axis image