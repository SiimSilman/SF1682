% Transponerad fvel för att anpassa till ode45
function ds = fvel(t, s, b, aL, aR, wL, wR)
    A = (aR + aL) / 2;
    B = (wR + wL) / 2;
    C = (aR - aL) / (2*b);
    D = (wR - wL) / b;

    theta = s(3);
    dx = (A*t + B)*cos( theta );
    dy = (A*t + B)*sin( theta );
    d0 = 2*C*t + D;
    
    ds = [dx, dy, d0];
end