%% P2A U2f
% Ny f(x) = exp(-(x-pi)^2)
% RV Fall 1: gamma = 0, delta = 0
% RV Fall 2: gamma = 1, delta = 0
clearvars, clc, close all

% Parametrar
gamma_list = [0, 1];
delta_list = [0, 0];
F     = @(x) exp(-(x-pi).^2);

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

hold on
for i = 1:length(gamma_list)
    % Hämta aktuell parameter
    gamma = gamma_list(i);
    delta = delta_list(i);
    
    % Definera högerled
    q      = F(xi)';
    q(1)   = q(1) + (2*h*gamma/3)*diagOffC(1);
    q(end) = q(end) - diagOffC(end)*delta;
    
    % Beräkning
    u      = A\q;
    plot(xi, u)
end
hold off
legend("Fall 1", "Fall 2")
xlabel("[0 2pi]")
ylabel("Amplitud")
title("Jämför två fall av RV")