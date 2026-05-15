%% P2B U3a
clearvars, clc, close all

% Parametrar
b  = 6;         % beta
N  = 16;        % # diskretiseringar
D = 1;

% Diskretisering
t0 = 0;
dt = (1/2)*10^-3;
tN = 0.1;
t  = 0:dt:tN;

k  = (-N/2: N/2-1)';
x  = (0:N-1)/N;

% Funktioner
F     = @(x) cos(b*pi*x);
G     = @(x, t) ((b*pi)^2 -1) * F(x) .* exp(-t); 
U     = @(x, t) F(x) .* exp(-t);                  

FT    = @(fj) fftshift(fft(fj));      % Fourier transform
IFT   = @(fk) ifft(ifftshift(fk));    % Invers Fourier transform
norm2 = @(f_ana, f_num) sqrt(mean(sum(abs( f_ana - f_num ).^2 ))); % d2-normen

% Initialvärden
u0 = U(x,0);
uk = (1/N) * FT(u0)';

% Allokera minne & Ansätt initialvärden
uj_matris = zeros(N, length(t));
uj_matris(:,1) = u0;

uk_matris = zeros(N, length(t));
uk_matris(:,1) = uk;

%% Lös Tidsrummet numeriskt (RK4)
for i = 1:length(t)-1
    gj = G(x, t(i));    % Hämta g(x,t)
    gk = (1/N)*FT(gj)'; % Hämta gk(t)

    uk_matris(:, i+1) = RK4(uk_matris(:,i), gk, t(i), dt, D, k); % uk(t)
    uj_matris(:, i+1) = N*IFT(uk_matris(:,i+1));                 % u(x,t) 
end

% Felet vid tN = 0.1s
fel = U(x', t(end)) - uj_matris(:, end);
fel_norm = norm2(U(x', t(end)), uj_matris(:, end));
plot(x, fel)
title(['Fel, d2norm = ', num2str(fel_norm)])
xlabel("x diskretiserad m.a.p N")
ylabel("fel vid t = 0.1")

% 3D plott av rum och tid
figure
surf(uj_matris)
colormap default
shading interp
title("Tid och rum")
xlabel("tid 0:1/2:100 [ms]")
ylabel("(0:N-1)/N")

%% RK4 funktion
function uk_next = RK4(uk, gk, t, dt, D, k)
    % du/dt = u_prim
    Uk_prim  = @(t, uk) -D*(2*pi.*k).^2 .* uk + gk;

    % Kurva
    k1 = Uk_prim(t, uk);                  % Början av intervall
    k2 = Uk_prim(t + dt/2, uk + k1*dt/2); % Mitten av intervall
    k3 = Uk_prim(t + dt/2, uk + k2*dt/2); % Mitten av intervall
    k4 = Uk_prim(t + dt, uk + dt*k3);     % Slutet av intervall
    
    % u_k i nästa tidssteg
    uk_next = uk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end