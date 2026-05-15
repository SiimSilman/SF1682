%% P2A U4e
%  Konvergensordning för Crank–Nicolson med dt = C*dx
clearvars, clc, close all

% Parametrar
D    = 1;
T_L  = @(t)   0; 
T_R  = @(t)   sin(exp(t));
F    = @(x,t) exp(t).*x.*cos(exp(t).*x) + exp(2*t).*sin(exp(t).*x);
U0   = @(x)   sin(x);
W    = @(x,t) sin(exp(t).*x);

% Tidsintervall
t0 = 0;
tM = 1;
C = 0.5; % Välj C (dt = C*dx)

% Lista N för konvergensstudie
N_list = [25, 50, 100, 200];

% Allokera minne
r2 = zeros(size(N_list));
dx_list = zeros(size(N_list));

for j = 1:length(N_list)
    % Diskretisering av rum
    N  = N_list(j);
    x0 = 0; xN = 1;
    dx = (xN - x0)/N;
    x  = linspace(x0, xN, N+1);
    xi = x(2:end-1)';
    
    dx_list(j) = dx;
    
    % Diskretisering av tid
    dt = C * dx;
    T  = t0:dt:tM;
    
    % Matris
    sigma = (dt * D) / (2 * dx^2);

    A = diag( 1 + 2*sigma * ones(N-1,1) ) ...
      + diag( -sigma * ones(N-2,1),  1 ) ...
      + diag( -sigma * ones(N-2,1), -1 );

    B = diag( 1 - 2*sigma * ones(N-1,1) ) ...
      + diag(  sigma * ones(N-2,1),  1 ) ...
      + diag(  sigma * ones(N-2,1), -1 );

    % Nummerisk lösning
    U = zeros(N-1, length(T));
    U(:,1) = U0(xi);

    for i = 1:length(T)-1
        f1 = F(xi, T(i));
        f2 = F(xi, T(i+1));

        HL = B * U(:,i) + (dt/2) * (f1 + f2);

        HL(1)   = HL(1)   + sigma*(T_L(T(i)) + T_L(T(i+1)));
        HL(end) = HL(end) + sigma*(T_R(T(i)) + T_R(T(i+1)));

        U(:,i+1) = A \ HL;
    end

    % Analytisk lösning vid slutet
    w_end = W(xi, tM);

    % Diskret 2-norm
    r = U(:,end) - w_end;
    r2(j) = sqrt( (1/(N-1)) * sum( r.^2 ) );
end

% Konvergensordning
p = zeros(length(N_list)-1,1);
for k = 1:length(p)
    p(k) = log(r2(k)/r2(k+1)) / log(dx_list(k)/dx_list(k+1));
end

% Utskrift
disp('   N          dx        Error     p');
disp([N_list' dx_list' r2' [NaN; p]])
