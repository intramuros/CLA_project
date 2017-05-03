function A = make_matrix(m,n, coherence)

if coherence == 0
    rng(11);
    U = orth(rand(m, n));
    S = diag(linspace(1, 1e5, n));
    V = orth(rand(n));
    A = U*S*V';    
elseif coherence == 1
    rng(11);
    A = [diag(linspace(1,1e5,n)); zeros(m-n,n)];
    A = A + 1e-8*ones(m,n);    
else
    error('Coherence must equal 0 or 1');
end      
    

end