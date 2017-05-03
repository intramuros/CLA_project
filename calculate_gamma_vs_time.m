function [timevec_minres, timevec_lsqr] = calculate_gamma_vs_time(A, b, gammavec, average_over)
   %% If iterations is > 1 calculates average time over iterations
   timevec_minres = zeros(length(gammavec), length(itervec));
   timevec_lsqr = zeros(length(gammavec), length(itervec));
   
if average_over < 1
    error('average_over must be larger or equal 1');
end
   
   
for iter = 1:average_over
    for i = 1:length(gammavec)
        tic;
        [~, ~, ~, ~] = blendenpik(A,b, gammavec(i), 'MINRES');
        timevec_minres(i, iter) = toc; 
        tic;
        [~, ~, ~, ~] = blendenpik(A,b, gammavec(i), 'LSQR');
        timevec_lsqr(i, iter) = toc; 
    end  
end


end