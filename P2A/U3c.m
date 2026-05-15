%% P2A U3c
% Jämför Analytisk och Nummerisk lösning, om korrekt genomförd
% förväntas högst på fel ca 10^-4 med dt = 10^-5;
clearvars, clc, close all

% Parametrar
D    = 1;
T_L  = @(t) 0; 
T_R  = @(t) 0;
F    = @(x,t) 0;
U0   = @(x) sin(5*pi*x);
W    = @(x,t) sin(5*pi*x) .* exp(-D*(5*pi)^2 * t);

% Diskretisering av tid
dt     = 10^-5;
T_slut = 0.01;
T      = 0:dt:T_slut;

% Diskretisering av rum
N      = 100;
x0     = 0;
xN     = 1;
dx     = (xN - x0) / N;
x      = linspace(x0, xN, N+1);
xi     = x(2:end-1);

% Matris
sigma = D / dx^2;
A     = diag( -2*sigma * ones(N-1, 1)) ...  % Off-Center
      + diag( sigma * ones(N-2, 1), 1) ...  % Center
      + diag( sigma * ones(N-2, 1), -1);    % Off-Center

% Nummerisk lösning: U genom rum (y-led) och tid (x-led)
U = zeros(length(xi), length(T));

U(:, 1) = U0(xi);
for i = 1:length(T)-1
   U(:, i+1) = U(:, i) + dt*( A*U(:,i) );
end

% Felet
xi_col = xi(:);
err = W(xi_col, T_slut) - U(:, end);
max(err)

% 3D plott nummerisk lösning
surfl(U)
colormap(pink)
shading interp
title("Nummerisk lösning")
xlabel("tid")
ylabel("rum")
zlabel("temperatur")

% Plott av felet
figure;
plot(xi, err)
ylabel("error amplitud")
xlabel("rum [0 1]")