clc
clearvars
clf
%% Initialdata
% Konstanter
k1ref = [5e2,       5350,       5e4,        5e5];  
k2ref = [13.6e2,    13.6e3,     13.6e4,     136100]; 
ck = 0.7;
cs = 1;
tolerans = 1e-6;

for j = 1:length(k1ref)                     % Testa olika startvärden
    %% Newton-Raphson med tolerans 1e-6
    k = [k1ref(j), k2ref(j)]';
    konv = [];
    error = [];
    error(1) = 1;
    i = 1;
    while error(i) > tolerans
        i = i + 1;                                      % Iteration
    
        F = transfer_functions(k(1), k(2), ck, cs);     % Hämta övergångsfnk
        J = Jacobian_transfer_functions(k(1), k(2));    % Hämta jacobianen
    
        error(i) = norm(k - (k - J \ F));               % Räkna felet
        k = k -J \ F;                                   % Räkna förbättrad approx
    
        if i >= 3                                       % Räkna konvergensen
            konv(i-2) = log(error(i) / error(i-1)) / log(error(i-1) / error(i-2));
        end
    end
    
    %% Kommentar
    fprintf('Initalvärden k1_ref = %i och k2_ref = %i \n', k1ref(j), k2ref(j))
    fprintf('Nya k1 värde blir %i kg/s^2\n', round(k(1)))
    fprintf('Nya k2 värde blir %i kg/s^2\n', round(k(2)))
    disp('Konvergensen blir')
    fprintf('%.4f \n', konv)
    fprintf('\n\n')
end