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

k1 = k1ref;
k2 = 100*k2ref;

% Initialvärden
Tlap = L/v + 0.1;
z1 = 0; 
z2 = 0; 
dz1 = 0; 
dz2 = 0; 
V0 = [z1, z2, dz1, dz2]';

%% Tag fram dt_max
% Systematris, H
rad1        = [0,          0,                  1,          0           ];
rad2        = [0,          0,                  0,          1           ];
rad3        = [-k1/m1,     k1/m1,              -c1/m1,     c1/m1       ];
rad4        = [k1/m2,      -(k1+k2)/m2,        c1/m2,      -(c1+c2)/m2 ];
H_matris    = [rad1; rad2; rad3; rad4];

% Egenvärden
eigvalue = eig(H_matris);

% Printa dt_max
for i = 1:4
    t_max = 2 * abs(real(eigvalue(i))) / abs(eigvalue(i)^2);
    fprintf('λ%i = %.4f %+.4fi, Re(λ) = %.4f\n', i, real(eigvalue(i)), imag(eigvalue(i)), real(eigvalue(i)));
    fprintf('t_max%i = %.4f\n', i, t_max)
end

%% Räkna med dt_max
a = [1, 10, 100];
dt_max = 0.0001;
dt = dt_max*a;

hold on
for j = 1:3
    V = V0;                     % Nollställ V
    tid = (0:dt(j):Tlap)';      % Anpassa längden tidslistan
    N = length(tid);            % Anpassa antal tidssteg
    z2_EF = zeros(1,N);          % Allokera lista

    for i = 1:N               
        V = trapets(tid(i), dt(j), V, m1, m2, k1, k2, c1, c2, H, L, v);
        z2_EF(i) = V(2);
    end
    plot(tid, z2_EF)
end
title('dt = dt_{max} = 0.0001')
legend('dt', '10*dt', '100*dt')
hold off


%% Funktion för trapets metoden
function dV = trapets(t, dt, V, m1, m2, k1, k2, c1, c2, H, L, v)
    
    if t <= L/v
        h = (H/2) * (1-cos(2*pi*v*t/L));
        dh = (H/2) * (2*pi*v/L) * sin(2*pi*v*t/L);
    elseif t > L/v
        h = 0;
        dh = 0;
    end
    g0 = [0, 0, 0, (k2*h + c2*dh)/m2]';

    t = t + dt;
    if t <= L/v
        h = (H/2) * (1-cos(2*pi*v*t/L));
        dh = (H/2) * (2*pi*v/L) * sin(2*pi*v*t/L);
    elseif t > L/v
        h = 0;
        dh = 0;
    end
    g1 = [0, 0, 0, (k2*h + c2*dh)/m2]';

    rad1 = [0, 0, 1, 0];
    rad2 = [0, 0, 0, 1];
    rad3 = [-k1/m1, k1/m1, -c1/m1, c1/m1];
    rad4 = [k1/m2, -(k1+k2)/m2, c1/m2, -(c1+c2)/m2];
    H = [rad1; rad2; rad3; rad4];

    I = eye(4);

    dV = (I - dt/2*H) \ (V + dt/2*(H*V + g0 + g1));
end