%% Problem set up
close all;
clear all

m=40000;
n= 800;

%Parameters to change
incoherent =1; %switch between cases
gamma = 4;     % Tune for gamma incoherent 1.5 91 iterations

% Construct matrices
b= rand(m,1);
    %% Incoherent, ill-conditioned matrix
    rng(11);
    U = orth(rand(m, n));
    S = diag(linspace(1, 1e5, n));
    V = orth(rand(n));
    A = U*S*V';
    
    
    %% Solve incoherent


[x, conv_rateMR_in] = blendenpik_iter(A,b, gamma, 'MINRES');


[x, conv_rateLSQR_in] = blendenpik_iter(A,b, gamma, 'LSQR');
    
    %% Coherent, ill-conditioned matrix
    rng(11);
    A = [diag(linspace(1,1e5,n)); zeros(m-n,n)];
    A = A+ 1e-8*ones(m,n);

    %% Solve coherent


[x, conv_rateMR_cor] = blendenpik_iter(A,b, gamma, 'MINRES');


[x, conv_rateLSQR_cor] = blendenpik_iter(A,b, gamma, 'LSQR');




%% Solve with MINRES


[x, conv_rateMR] = blendenpik_iter(A,b, gamma, 'MINRES');


[x, conv_rateLSQR] = blendenpik_iter(A,b, gamma, 'LSQR');
% Convergence in how many steps?
%% PLOTS
set(0,'DefaultAxesFontSize',12)
close all;
semilogy(conv_rateMR_in, '-x', 'LineWidth', 1.5);
hold on;
semilogy(conv_rateLSQR_in, '-x', 'LineWidth', 1.5)
semilogy(conv_rateMR_cor, '-x', 'LineWidth', 1.5)
semilogy(conv_rateLSQR_cor, '-x', 'LineWidth', 1.5)

grid on;
%set(gca,'yminorgrid','off')
legend('Incoherent- MINRES', 'Incoherent- LSQR', 'Coherent- MINRES', 'Coherent- LSQR')
xlabel('Iteration');
ylabel('Norm');
axis([0 80 1e-17 1])

