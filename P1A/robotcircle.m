function [xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0)
    %A = (aR + aL)/2;
    B = (wR + wL)/2;
    %C = (aR - aL)/(2*b);
    D = (wR - wL)/b;

    R = B/D;
    Tlap = 2*pi/D;
    xbar = x0 + B/D * sin(theta0);
    ybar = y0 - B/D * cos(theta0);
end