%% P2A U4f
% Analysera strukturen vid senare tidssteg
clearvars, clc, close all

% Parametrar
D    = 1;
T_L  = @(t)   0; 
T_R  = @(t)   sin(exp(t));
F    = @(x,t) exp(t).*x.*cos(exp(t).*x) + exp(2*t).*sin(exp(t).*x);
U0   = @(x)   sin(x);
W    = @(x,t) sin(exp(t).*x);

% Diskretisering av tid
dt     = 10^-3;
T_slut = 4;
T      = 0:dt:T_slut;

% Diskretisering av rum
N      = 100;
x0     = 0;
xN     = 1;
dx     = (xN - x0) / N;
x      = linspace(x0, xN, N+1);
xi     = x(2:end-1)';

% Matris
sigma = (dt * D) / (2 * dx^2);
A     = diag( 1+2*sigma * ones(N-1, 1)) ...  % Center
      + diag( -sigma * ones(N-2, 1), 1) ...  % Off-Center
      + diag( -sigma * ones(N-2, 1), -1);    % Off-Center

B     = diag( 1-2*sigma * ones(N-1, 1)) ...  % Center
      + diag( +sigma * ones(N-2, 1), 1) ...  % Off-Center
      + diag( +sigma * ones(N-2, 1), -1);    % Off-Center

% Nummerisk lösning
U = zeros(length(xi), length(T));   % Allokera minne
U(:, 1) = U0(xi);                   % Begynnelse villkor

for i = 1:length(T)-1
    f1 = F(xi, T(i));
    f2 = F(xi, T(i+1));

    HL = B * U(:,i) + (dt/2) * (f1 + f2);

    HL(1) = HL(1) + sigma*( T_L(T(i)) + T_L(T(i+1)) );
    HL(end) = HL(end) + sigma*( T_R(T(i)) + T_R(T(i+1)) );
    
    U(:, i+1) = A \ HL;
end

% Analytisk lösning
w = W(xi, T);

% Felet
r = abs(w(:, end) - U(:, end));
max(r)

% Resultat
surfl(U)
colormap(pink)
shading interp
title("Nummerisk lösning")
xlabel("tid")
ylabel("rum")
zlabel("temperatur")
