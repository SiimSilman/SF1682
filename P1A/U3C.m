clc    
clearvars  
close all
%% U3 C
% Utifrån given initialvärden:
% - approximera körbanan med Euler-Fram med 100 tidsteg

%% Initial värden
b = 1;  
aL = 0; aR = 0; wL = 4; wR = 2;

% Initial villkor
theta0 = 0; x0 = 0; y0 = 0;

% Hämta startkordinater, radien och tiden för ett varv
[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0);

s = [xbar, ybar, theta0];

%% Fördela tiden med 100 tidsteg
N = 100;
h = Tlap/N;
tid = 0:h:Tlap;

%% Allokera minne för sparade värden
dx = ones(N+1, 1);
dy = zeros(N+1, 1);
dtheta = zeros(N+1, 1);

%% Euler Fram
for i = 1:N+1
    dx(i) = s(1);
    dy(i) = s(2);
    dtheta(i) = s(3);
    
    ds = fvel(tid(i), s, b, aL, aR, wL, wR);
    s = s + ds*h;
end

%% Plotta
hold on
plot(dx, dy, 'r'); 
plot(0, 0, 'o'); 
hold off

legend('EF', 'Mittp')
axis equal
axis image





