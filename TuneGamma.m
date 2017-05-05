%% Problem set up
close all;
clear all

m=40000;
n= 2000;

%Parameters to change
incoherent =1; %switch between cases
gamma = 8;     % Tune for gamma incoherent 1.5 91 iterations
p= gamma*n/m
% Construct matrices
if incoherent
    %% Incoherent, ill-conditioned matrix
    rng(11);
    U = orth(rand(m, n));
    S = diag(linspace(1, 1e5, n));
    V = orth(rand(n));
    A = U*S*V';
    
else
    %% Coherent, ill-conditioned matrix
    rng(11);
    A = [diag(linspace(1,1e5,n)); zeros(m-n,n)];
    A = A+ 1e-8*ones(m,n);
end


b= rand(m,1);


%% Solve with MINRES
tic
[x, flag, iter, resvec] = blendenpik(A,b, gamma, 'MINRES');
toc

tic
x_backslash =A\b;
toc

rminres= norm(b-A*x)/norm(b);
r_backslash = norm(b-A*x_backslash)/norm(b);

% Convergence in how many steps?
semilogy(resvec)

