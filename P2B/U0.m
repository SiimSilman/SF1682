% Fråga: Vad menas med ett exponentiellt avtagande? 
% I en log-log plot: 
% 1) plotta e^(−qN) mot N för något q 
% 2) plotta N^(−q) mot N. 
% Hur ser kurvorna ut?

clearvars, clc, close all

q = 3;
N = 1:10;

f1 = @(N) exp(-q.*N);
f2 = @(N) N.^(-q);

subplot(2,1,1)
loglog(N, f1(N))
title("loglog e^{-qN}")

subplot(2,1,2)
loglog(N, f2(N))
title("loglog N^{-q}")
xlabel("N = 1:10")

% Med epononetiellt avtagande menas:
% Att potentialvärdet ökas exponentiellt, vilket syns i översta figur.
% I nedersta figur ser vi potentialvärdet är konstant.