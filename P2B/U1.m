% P2B U1
clearvars, clc, close all

% Parametrar
m = 4:8;
N = 2.^m;
a = 400;

% Funktioner, analytiska lösningar
F      = @(x) exp(-a * (x - pi/5).^2);
F_prim = @(x) F(x) .* (-2*a*(x - pi/5)).^1;
F_bis  = @(x) F(x) .* (-2*a*(x - pi/5)).^2 + F(x) .* (-2*a); 

% Funktioner, fourier transform
Fk      = @(fj, N) (1/N) * fftshift(fft(fj)); % Fourier transform
Fk_prim = @(fk, k) (1i*2*pi.*k).^1 .*fk;
Fk_bis  = @(fk, k) (1i*2*pi.*k).^2 .*fk;
Fk_i    = @(fk, N) N * ifft(ifftshift(fk)); % Invers Fourier transform

% Funktion, diskreta 2-fel
d2fel   = @(f_ana, f_num, N) sqrt( (1/N) .* sum( abs( f_ana - f_num ).^2 ) ); 

% Allokera minne
fel_prim = zeros(1, length(N));
fel_bis = zeros(1, length(N));

for i = 1:length(N)
    % Diskretisera
    Ni = N(i);             
    xj = (0:Ni-1) / Ni;    
    k  = -Ni/2 : Ni/2 - 1;  
    
    % Analytisk lösning
    fj      = F(xj);  
    fj_prim = F_prim(xj);
    fj_bis  = F_bis(xj); 
    
    % Fourier transform
    fk      = Fk(fj, Ni);          
    fk_prim = Fk_prim(fk, k); 
    fk_bis  = Fk_bis(fk, k);  

    % Invers Fourier transform
    fk_prim_i = Fk_i(fk_prim, Ni);
    fk_bis_i  = Fk_i(fk_bis, Ni);

    % Diskreta 2-fel
    fel_prim(i) = d2fel(fj_prim, fk_prim_i, Ni);
    fel_bis(i)  = d2fel(fj_bis, fk_bis_i, Ni);
end

% Plot av faktiska felet
hold on
plot(N, fel_prim, '-o')
plot(N, fel_bis, '-x')
hold off
legend("prim", "bis")
title("Diskreta 2-fel")
xlabel("N")
ylabel("Fel")
xlim([0 256])

% Plot av loggade felet
figure;
hold on
plot(log2(N), log10(fel_prim), '-o')
plot(log2(N), log10(fel_bis), '-x')
hold off
legend("prim", "bis")
title("Diskreta 2-fel")
xlabel("log_2(N)")
ylabel("log_{10}(Fel)")