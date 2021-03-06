function [timevec_minres, timevec_lsqr] = calculate_gamma_vs_time(A, b, gammavec, average_over)
   %% If average_over is > 1 calculates average time over iterations
   
   
   timevec_minres = zeros(length(gammavec), average_over);
   timevec_lsqr = zeros(length(gammavec), average_over);
   
if average_over < 1
    error('average_over must be larger than or equal to 1');
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

if average_over > 1
   timevec_minres = mean(timevec_minres, 2);
   timevec_lsqr = mean(timevec_lsqr, 2);
end