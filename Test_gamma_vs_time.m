close all;
clear all;

% set matrix dimensioins
m = 20000;
n = 400;

% set coherence
is_coherent = true;

% Make an ill-conditioned matrix
A = make_matrix(m,n, is_coherent);


b = rand(m,1); 

% Create an array of different gammas

gamma_start = 2;
gamma_max = 10;
gamma_step = 1;

gammavec = gamma_start:gamma_step:gamma_max;
gammavec = [1.09 gammavec];

%
average_over = 1;

[timevec_minres, timevec_lsqr] = calculate_gamma_vs_time(A, b, gammavec, average_over);


%% Plotting 

if is_coherent == false 
    coh_str = ' incoherent ';
elseif is_coherent == true
    coh_str = ' coherent ';
end

size_str = [int2str(m) '-by-' int2str(n)];
title_str = [size_str coh_str ', ill-conditioned matrix'];

plot(gammavec, timevec_minres, '-o', gammavec, timevec_lsqr, '-*');
grid on
title(title_str);
xlabel('\gamma'); ylabel('Time (sec)'); 
legend('MINRES', 'LQSR');
