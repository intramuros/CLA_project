
%   b) Calculates average coherence of 1000 random matrices
%   
%
mu_vect = zeros(1000,1);

for i=1:1000
   A = rand(1000, 50); 
   mu_vect(i,1) = coherence(A);
end

average_coherence = mean(mu_vect);