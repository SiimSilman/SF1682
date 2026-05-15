clc 
clearvars
close all
%% U6B: Startgissningsvärde och förbättring med Newton Raphson
% Samma som U6A fast T = 2

%% Initialdata
h = 1;
aR = (0:h:10);
aL = 0.4;
wR = 3; 
wL = 1; 
b = 1; 

%% Funktioner
RT = 1;
R   = @(aR) b.*( 2*(aR+aL)+(wR+wL) )    ./ (2.*( 2*(aR-aL)+(wR-wL) ) ) - RT;
dR  = @(aR) -2*b.*(2*aL+wL)             ./ ( 2*(aR-aL)+(wR-wL) ).^2;

%% Hitta start-gissnings-värde
Rlist = R(aR);                          % alla y-led värden
[Rmin, Rmin_index] = min(abs(Rlist));   % index där y är närmast 0
aR_start = aR(Rmin_index);              % index i aR där y är närmast 0

%% Newton-Raphson
aRlist = zeros();   aRlist(1) = aR_start;
error = zeros();
konvergens = zeros();

for i = 1:5
    aRlist(i+1) = aRlist(i) - R(aRlist(i)) / dR(aRlist(i)); % NR
    error(i) = abs(aRlist(i+1)-aRlist(i));                  % Fel
    
    if i >= 3
        konvergens(i-2) = log(error(i)/error(i-1)) / log(error(i-1)/error(i-2));
    end
end

%% Plottar
fprintf('R(aR) - R(T) så nära 0 som möjligt när aR = %.2f\n', aR_start);
fprintf('Förbättrad värde på aR = %.2f\n', aRlist(end));
fprintf('Konvergensen för NR är: \n')
fprintf('%.2f\n', konvergens(1,:))

hold on
    plot(aR, Rlist)
    plot(aR_start, Rmin, 'x')
    legend('R(a_R) - R(T)', 'R(a_R) nära 0')
    xlabel('a_R')
    ylabel('R(a_R) - R(T)')
    title('aR fördelad som 1:1:10')
hold off