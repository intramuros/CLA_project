% Script regarding the concept of coherence. Average coherence of random
% matrices
close all;
clear all;
n=1000; % number of matrices

avg_coherence=0;

for i=1:1000
    A = rand(1000,50);
    [Q, R] = qr(A,0);
    cohr(i) = max(sum (Q .^2,2));
    avg_coherence = avg_coherence +cohr(i);
end

avg_coherence= avg_coherence/n;


%% PLOTS
set(0,'DefaultAxesFontSize',12)
close all;
hist(cohr)
xlabel('Coherence \mu(A)')