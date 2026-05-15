%% P2A U2c
% Ny N = [50, 100, 200, 400, 800]
% Diskreta 2-felet (rms) och noggrannhetsordning
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

N_list = [50, 100, 200, 400, 800];  % N lista med olika antal diskretiseringar
r2     = zeros(1,length(N_list));   % Allokera minne åt 2-felet
 
for i = 1:length(N_list)
    % Diskretisering
    N     = N_list(i); 
    x0    = 0;
    xN    = 2*pi;
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
    
    % Diskreta 2-felet (rms av felet)
    r      = w(2:end-1) - u;
    r2(i)  = ( 1/(N-1) * sum(r.^2) )^(1/2);

    % Noggrannhetsordning + Print
    if i >= 2
        p = log(r2(i-1) / r2(i)) / log(2);
        fprintf("N = %i -> r2 = %.5f -> p = %.4f \n", N, r2(i), p)
    else
        fprintf("N = %i  -> r2 = %.5f \n", N, r2(i))
    end
end