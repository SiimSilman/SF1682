% P2B U2c
clearvars, clc, close all

% Parametrar
D = 1;
g = 0;
N = 16;
t = 0.1;

% Diskretisering          
x = (0:N-1) / N;    
k = -N/2:N/2-1;  

% Funktioner
F       = @(x) sin(2 * pi .*x);                % Begynnelse villkor
U       = @(x) F(x) .* exp(-D*(2*pi)^2*t);     % Analytisk lösning
FT      = @(fj, N) (1/N) * fftshift(fft(fj));  % Fourier transform
iFT     = @(fk, N) N * ifft(ifftshift(fk));    % Invers Fourier transform

% Genomförande
fj = F(x);
uj = U(x);

fk   = FT(fj, N);
uk   = fk .* exp(-D*(2*pi.*k).^2*t);
uk_i = iFT(uk, N);

fel = abs(uj - uk_i);


hold on
plot(x, uj, '-o')
plot(x, uk_i, '-x')
hold off
legend("Analytisk", "Fourier")
xlabel("x")
ylabel("amplitud")
title("Analytisk lösning och FFT")

figure;
plot(x, fel)
title("Fel")
xlabel("x")
ylabel("Fel")
