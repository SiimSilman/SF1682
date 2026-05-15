clc
clearvars
clf
%% Initialdata
% Konstanter
k1ref = 5350;  
k2ref = 136100; 
ck = 0.7;
cs = 1;
tolerans = 1e-6;

%% Newton-Raphson med tolerans 1e-6
k = [k1ref, k2ref]';
error = 1;
while error > tolerans
    F = transfer_functions(k(1), k(2), ck, cs);
    J = Jacobian_transfer_functions(k(1), k(2));

    error = norm(k - (k - J \ F));
    k = k -J \ F;
end

%% Kommentar
% Optimerad k1 och k2 värde sådan att komforten maximeras samtidigt som
%   körsäkerheten inte kompromissas.
fprintf('Initalvärdet k1_ref = %i kg/s^2\n', k1ref)
fprintf('Initalvärdet k2_ref = %i kg/s^2\n', k2ref)
fprintf('Optimerade värde k1 värde blir %i kg/s^2 \n', round(k(1)))
fprintf('Optimerade värde k2 värde blir %i kg/s^2 \n', round(k(2)))
