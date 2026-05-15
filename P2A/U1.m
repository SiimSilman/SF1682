%% P2A U1b
clearvars, clc, close all

% Definera parametrar
x = 0.75;               % Utgångspunkt
k = 1:1:8;              % Faktor för steglängds multiplikation

% Definera funktioner
H     = @(k) 2.^(-k);                   % Steglängd 
U     = @(x) sin(exp(x));               % Funktion
Uprim = @(x) cos(exp(x)) * exp(x);      % Funktionens derivata

U2 = @(x, h) ( + 1*U(x+h)   - 1*U(x-h)              ) / (2*h); % Derivata approx (2)
U6 = @(x, h) ( - 1*U(x+2*h) + 4*U(x+h) - 3*U(x)     ) / (2*h); % Derivata approx (6)
U7 = @(x, h) ( + 3*U(x)     - 4*U(x-h) + 1*U(x-2*h) ) / (2*h); % Derivata approx (7)

% Allokera minne
err2 = zeros(1, length(k));
err6 = zeros(1, length(k));
err7 = zeros(1, length(k));

%% Genomförande
ref = Uprim(x);

for i = k
    % Hämta aktuell steglängd
    h = H(i);

    % Bestäm felet mellan referens och approx
    err2(i) = abs(ref - U2(x, h));
    err6(i) = abs(ref - U6(x, h));
    err7(i) = abs(ref - U7(x, h)); 
end

hold on
loglog(H(k), err2)
loglog(H(k), err6)
loglog(H(k), err7)
hold off
xlabel("stepsize, log(h)")
ylabel("error, log( |ref - approx| )")
legend("err2", "err6", "err7")