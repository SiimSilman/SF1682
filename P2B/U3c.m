%% P2B U3c
clearvars, clc, close all

% Parametrar
N  = 16;        
D = 1;
delta = 0.05;

% Diskretisering
t0 = 0;
dt = 1.14*10^-3;
tN = 0.1;
t  = 0:dt:tN;

k  = (-N/2: N/2-1)';
x  = (0:N-1)/N;

% Funktioner
F     = @(x) 0;
G     = @(x, t) 100 * eta(2*abs(x-1/2), delta) .* exp(-100*t); 
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

% 3D plott av rum och tid
figure
waterfall(uj_matris)
colormap default
shading interp
title("dt = 1.14ms, N = 16")
xlabel("rum")
ylabel("tid")


%% Funktioner

%  RK4 funktion
function uk_next = RK4(uk, gk, t, dt, D, k)
    Uk_prim  = @(t, uk) -D*(2*pi.*k).^2 .* uk + gk; % du/dt = u_prim

    % Kurva
    k1 = Uk_prim(t, uk);                  % Början av intervall
    k2 = Uk_prim(t + dt/2, uk + k1*dt/2); % Mitten av intervall
    k3 = Uk_prim(t + dt/2, uk + k2*dt/2); % Mitten av intervall
    k4 = Uk_prim(t + dt, uk + dt*k3);     % Slutet av intervall
    
    uk_next = uk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4); % u_k i nästa tidssteg
end

% Eta funktion
function out = eta(y, delta)
    h = @(s) exp((2*exp(-1./(s+eps))) ./ (s - 1 - eps));
    a = 1 - 2*delta;
    out = zeros(size(y));

    mask1 = (y <= delta);
    mask3 = (y >= 1-delta);
    mask2 = (~mask1 & ~mask3);

    out(mask1) = 1;
    out(mask3) = 0;

    if any(mask2)
        s = (y(mask2) - delta)/a;
        out(mask2) = h(s);
    end
end