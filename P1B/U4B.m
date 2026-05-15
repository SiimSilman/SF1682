clc
clearvars
clf
%% Initialdata
m1 = 465; 
m2 = 55;
k1ref = 5350; 
k2ref = 136100; 
c1 = 310; 
c2 = 1250; 
v = 63/3.6; 
L = 1.1; 

%% Systematris, H
rad1 = [0,          0,                  1,          0           ];
rad2 = [0,          0,                  0,          1           ];
rad3 = [-k1ref/m1,  k1ref/m1,           -c1/m1,     c1/m1       ];
rad4 = [k1ref/m2,   -(k1ref+k2ref)/m2,  c1/m2,      -(c1+c2)/m2 ];
H = [rad1; rad2; rad3; rad4];

%% Egenvärden, dt_max tas fram
eigvalue = eig(H);
for i = 1:4
    t_max = 2 * abs(real(eigvalue(i))) / abs(eigvalue(i)^2);
    fprintf('λ%i = %.4f %+.4fi, Re(λ) = %.4f\n', i, real(eigvalue(i)), imag(eigvalue(i)), real(eigvalue(i)));
    fprintf('t_max%i = %.4f\n', i, t_max)
end