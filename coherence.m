function mu = coherence(A)
%COHERENCE Coherence of a matrix A
%
%   mu = coherence(A) returns the coherence of matrix A.
%   The columns of A are first orthonormalized.
%
A = orth(A);

row_norms = sqrt(sum(abs(A).^2,2));
mu = max(row_norms);


end