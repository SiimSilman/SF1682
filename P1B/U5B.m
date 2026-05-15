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
z1 = 0; 
z2 = 0; 
dz1 = 0; 
dz2 = 0; 
V0 = [z1, z2, dz1, dz2]';

% Tiden
Tlap = 0.05;
dtmax = 0.001;
a = [1, 1/2, 1/4, 1/8]';
dt = a*dtmax;

hold on
error = zeros(1,length(a));
for j = 1:length(a)
    %% Approx: Trapets, dt = dt_max*a
    tid = (0:dt(j):Tlap)';
    N = length(tid);
    z2_Tr = zeros(N,1);

    %% Approx: Trapets 
    V = V0;
    for i = 1:N    
        z2_Tr(i) = V(2);
        V = trapets(tid(i), dt(j), V, m1, m2, k1, k2, c1, c2, H, L, v);
    end
    
    %% Exakt: ode45
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Refine', 1);
    [t_ode, V_ode] = ode45(@(t,V) quartercar(t, V, m1, m2, k1, k2, c1, c2, H, L, v), tid, V0, options);
    z2_ode_interp = interp1(t_ode, V_ode(:,2), tid);
    
    %% Plot och Error
    error(j) = max(abs(z2_Tr - z2_ode_interp));
    
    plot(tid, z2_Tr)
    if j == length(a)
        plot(t_ode, V_ode(:,2))
        legend('dt', 'dt/2', 'dt/4', 'dt/8', 'ode45')
    end
end

%% Konvergensstudie
for i = 1:length(error)-1
    konvergens = log(error(i)/error(i+1)) / log(2);
    fprintf('Konvergens: %.4f \n', konvergens);
end

%% Funktion för trapets metoden (samma som förra)
function dV = trapets(t, dt, V, m1, m2, k1, k2, c1, c2, H, L, v)
    
    if t <= L/v
        h = (H/2)*(1-cos(2*pi*v*t/L));
        dh = (H*v*pi/L) * sin(2*pi*v*t/L);
    elseif t > L/v
        h = 0;
        dh = 0;
    end
    g0 = [0, 0, 0, (k2*h + c2*dh)/m2]';

    t = t + dt;
    if t <= L/v
        h = (H/2)*(1-cos(2*pi*v*t/L));
        dh = (H*v*pi/L) * sin(2*pi*v*t/L);
    elseif t > L/v
        h = 0;
        dh = 0;
    end
    g1 = [0, 0, 0, (k2*h + c2*dh)/m2]';

    Ar1 = [0, 0, 1, 0];
    Ar2 = [0, 0, 0, 1];
    Ar3 = [-k1/m1, k1/m1, -c1/m1, c1/m1];
    Ar4 = [k1/m2, -(k1+k2)/m2, c1/m2, -(c1+c2)/m2];
    A = [Ar1; Ar2; Ar3; Ar4];

    I = eye(4);

    dV = (I - dt/2*A) \ (V + dt/2*(A*V + g0 + g1));
end