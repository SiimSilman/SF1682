 function dV = quartercar(t, V, m1, m2, k1, k2, c1, c2, H, L, v)
    
    % Definitionsmängden av h(t), enligt ekv (2)
    if t <= L/v
        h   = (H/2) * (1 - cos(2*pi*v*t/L));           
        dh  = (H/2) * (2*pi*v/L) * sin(2*pi*v*t/L);     
    elseif t > L/v
        h   = 0;
        dh  = 0;
    end

    % Beräknad tidsberoende term g(t), enligt ekv (3)
    G = [0, 0, 0, (k2*h + c2*dh)/m2]';                 

    % Beräknad matris H, enligt ekv (3)
    rad1 = [0, 0, 1, 0];
    rad2 = [0, 0, 0, 1];
    rad3 = [-k1/m1, k1/m1, -c1/m1, c1/m1];
    rad4 = [k1/m2, -(k1+k2)/m2, c1/m2, -(c1+c2)/m2];
    H = [rad1; rad2; rad3; rad4];

    % Första ordningens system i tiden, enligt ekv (3)
    dV = H*V + G;
end