function [ x, flag,  iter, resvec ] = blendenpik( A, b, gamma, type)


if  strcmp(type, 'LSQR')
    MINRES_flag=0;
elseif strcmp(type, 'MINRES')
    MINRES_flag=1;
else
    disp('Error: type must be either "LSQR" or "MINRES"! ');
end
maxit_bp = 3;
maxit = 100; % randomly chosen, should ideally never be reached
tol= 1e-14;
if (size(A,2) ~= size(b))
    sprintf( 'Warning: Dimensions A and b must agree!');
    
end
n = size(A,2);
m = size(A,1);
m_tilde = ceil(m/1000)*1000;
M = [A; zeros((m_tilde-m),n) ]; % Matlabs DCT can also do padding
% What does the tuning business mean? What do you tune on? What do you
% optimize for?? Does it transfer from fftw to matlab? Or is it just speed
% up fftw

for i=1:maxit_bp
    D = spdiags(sign(rand(m_tilde,1)- 0.5), 0, m_tilde, m_tilde);
    
    
    M = dct(D*M);
    M(1,:)= M(1,:)/sqrt(2);
    
    prob = gamma*n/m_tilde;
    s=rand(m_tilde,1) < prob; % row sampling
    SM= M(s,:);
    
    [~, R] = qr(SM, 0);
    size(R)
    kappa_rec =rcond(R);
    if kappa_rec > 5*eps % Compare with Matlabs machine accuracy
        if(MINRES_flag)
            [x,flag,~,iter, resvec] = minres((A'*A),(A'*b),tol*norm(A'*b), maxit,R', R);
            return;
        else
            [x,flag,~,iter, resvec] = lsqr(A,b,tol,maxit,R);
            return;
        end
    end
    
    
    
end
end
