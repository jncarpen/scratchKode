%% find percentage of cells that are significant
HD_sig_NP = zeros(100,1); RH_sig_NP = zeros(100,1); 
for i = 1:100
    % grab the modulation strengths for each unit
    MS_HD_now = A(i).out.measures.TS.HD;
    MS_RH_now = A(i).out.measures.TS.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = A(i).shift.mshd;
    RH_shuf_now = A(i).shift.msr;
    
    % find confidence interval for the shuffled distribution
    HD_ci = [mean(HD_shuf_now) - 2*std(HD_shuf_now), mean(HD_shuf_now) + 2*std(HD_shuf_now)];
    RH_ci = [mean(RH_shuf_now) - 2*std(RH_shuf_now), mean(RH_shuf_now) + 2*std(RH_shuf_now)];
    
    HD_sig(i) = MS_HD_now < HD_ci(1) | MS_HD_now > HD_ci(2);
    RH_sig(i) = MS_RH_now < RH_ci(1) | MS_RH_now > RH_ci(2);
    
    % sort the shuffled distribution and real values
    HD_sort = sort([HD_shuf_now, MS_HD_now], 'ascend');
    RH_sort = sort([RH_shuf_now, MS_RH_now], 'ascend');
    
    % use a non-parametric method to determine significance
    if find(HD_sort == MS_HD_now) > 950 | find(HD_sort == MS_HD_now) < 50; HD_sig_NP(i) = 1; end
    if find(RH_sort == MS_RH_now) > 950 | find(RH_sort == MS_RH_now) < 50; RH_sig_NP(i) = 1; end
end

disp('COMPUTING MODEL SENSITIVITY MEASURES...')
disp(['HD TS (parametric):', num2str(sum(HD_sig)), '%']);
disp(['HD TS (non-parametric):', num2str(sum(HD_sig_NP)), '%']);
disp(['RH TS (parametric):', num2str(sum(RH_sig)), '%']);
disp(['RH TS (non-parametric):', num2str(sum(RH_sig_NP)), '%']);