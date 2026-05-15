%% P2A U2d
% Au = q -> u = A\q
% u RV kända, dvs lös okända u(2:end-1)
clearvars, clc, close all

% Parametrar
gamma = 1;
delta = 1;
d     = 5;

A     = (delta - d/12);
B     = gamma / 4;
C     = d/12;
U     = @(x) A*cos(4*x) + B*sin(4*x) + C*cos(2*x);
F     = @(x) d*cos(2*x);

% Diskretisering
x0    = 0;
xN    = 2*pi;
N     = 400;
h     = (xN - x0) / N;
x     = linspace(x0, xN, N+1);
xi    = x(2:end-1);

% System Matris (formel 3)
diagOffC = (1/h^2)      * ones(1, N-2);   % Diagonal, Off-Center
diagC    = (16 - 2/h^2) * ones(1, N-1);   % Diagonal, Center
A1       = diag(diagOffC, -1);
A2       = diag(diagC, 0);
A3       = diag(diagOffC, +1);
A        = A1 + A2 + A3;

% Neumann RV (formel 6)
A(1,1) = diagC(1) + (4/3)*diagOffC(1);
A(1,2) = diagOffC(1) - (1/3)*diagOffC(1);

% Dirichlet RV
beta   = U(2*pi);

% Definera högerled
q      = F(xi)';
q(1)   = q(1) + (2*h*gamma/3)*diagOffC(1);
q(end) = q(end) - diagOffC(end)*beta;

% Numerisk och Analytisk lösning
u      = A\q;
w      = U(x)';
error  = abs(u - w(2:end-1));

%% Plotta Numerisk och Analytisk lösning
figure;
hold on
plot(xi, u, '--x');
plot(x, w);
hold off
legend("Numerisk", "Analytisk")
xlabel("[0 2pi]")
ylabel("amplitud")

% Plotta felet mellan Numerisk och Analytisk lösning
figure;
plot(xi, error)
xlabel("[0 2pi]")
ylabel("amplitud")
title("Error = |Analytisk - Numerisk|")