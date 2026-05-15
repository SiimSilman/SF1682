clc
clearvars
clf
%% Initialdata
% Konstanter
m1 = 465; 
m2 = 55;
k1ref = 5350; 
k2ref = 136100; 
c1 = 310; 
c2 = 1250; 
v = 63/3.6; 
H = 0.27; 
L = 1.1; 

% Initialvärden
z1 = 0; 
z2 = 0; 
dz1 = 0; 
dz2 = 0; 
V0 = [z1, z2, dz1, dz2]';

% Tid
Tlap = L/v + 0.5;
a = [0.9, 1, 1.1, 1.5];
dt_max = 0.0111;
dt = dt_max * a;

%% Euler Framåt, dt_max = 0.0111*a
hold on
for j = 1:4                     % Iterera dt_max*a
    V = V0;                     % Nollställ V
    tid = (0:dt(j):Tlap)';      % Anpassa längden tidslistan
    N = length(tid);            % Anpassa antal tidssteg
    z2_EF = zeros(1,N)';        % Allokera lista
    for i = 1:N                 % Euler Fram
        dV = quartercar(tid(i), V, m1, m2, k1ref, k2ref, c1, c2, H, L, v);
        V = V + dt(j) * dV;
        z2_EF(i) = V(2);
    end
    plot(tid, z2_EF)
end
legend('dt_{max}*0.9', 'dt_{max}', 'dt_{max}*1.1', 'dt_{max}*1.5')
hold off