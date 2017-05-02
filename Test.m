%% Test runs
close all;
%clear all;
m=20000;
n= 400;

gamma = 4;
incoherent = 1;
%% Incoherent, ill-conditioned matrix
if incoherent
    rng(11);
    U = orth(rand(m, n));
    S = diag(linspace(1, 1e5, n));
    V = orth(rand(n));
    A = U*S*V';
    
    
    
    %% Coherent, ill-conditioned matrix
else
    rng(11);
    A = [diag(linspace(1,1e5,n)); zeros(m-n,n)];
    A = A+ 1e-8*ones(m,n);
end
b= rand(m,1);
x= blendenpik(A,b, gamma, 'MINRES');

rminres= norm(b-A*x)/norm(b);

x= blendenpik(A,b, gamma, 'LSQR');

rLSQR= norm(b-A*x)/norm(b);