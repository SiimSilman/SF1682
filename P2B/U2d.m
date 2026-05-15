% P2B U2c
clearvars, clc, close all

% Parametrar
D = 1;
t = 0.1;
delta = 0.05;

% Funktioner
F     = @(x) eta(2*abs(x-1/2), delta);       % Begynnelse villkor
FT    = @(fj, N) (1/N) * fftshift(fft(fj));  % Fourier transform
IFT   = @(fk, N) N * ifft(ifftshift(fk));    % Invers Fourier transform
norm2 = @(f_ana, f_num) sqrt(mean(sum(abs( f_ana - f_num ).^2 ))); % d2-normen

% Referens (hög-upplösning av numerisk lösning)
N_ref = 2^11;
x_ref = (0:N_ref-1) / N_ref;    
k_ref = -N_ref/2:N_ref/2-1;  

fj_ref = F(x_ref);
fk_ref = FT(fj_ref, N_ref);
uk_ref = fk_ref .* exp(-D * (2*pi*k_ref).^2 * t);
u_ref  = real(IFT(uk_ref, N_ref));

% Jämförelse (låg-upplösning av numerisk lösning)
m     = 4:9;
fel   = zeros(1, length(m));

for i = 1:length(m)
    % Diskretisera
    N_num = 2^m(i);
    x_num = (0:N_num-1) / N_num;
    k_num = -N_num/2:N_num/2-1;
    
    % Hämta funktioner
    fj_num = F(x_num);
    fk_num = FT(fj_num, N_num);
    uk_num = fk_num .* exp(-D * (2*pi*k_num).^2 * t);
    u_num  = real(IFT(uk_num, N_num));

    % Diskreta 2-normen
    steg = N_ref / N_num;
    u_ref_skalad = u_ref(1:steg:end);
    fel(i) = norm2(u_ref_skalad, u_num);
end

hold on
plot(x_ref, fj_ref)
plot(x_ref, u_ref)
legend("f(x)", "u_{ref}(x)")
title("f(x) = eta(2|x-1/2|, 0.05)")
xlabel("x")
ylabel("f(x)")

figure;
plot(log2(2.^m), log10(fel), '-x')
title("loglog av diskreta 2-normen")
xlabel("log_2(N)")
ylabel("log_{10}(fel)")

% Eta funktionen
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