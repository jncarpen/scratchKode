%% find percentage of cells that are significant
pop = placehd1000; fp = placehd1000FP;

for un = 1:length(fp)
    unitNow(un) = fp(un).B.unit;
end

HD_sig_NP = zeros(length(pop),1).*nan; RH_sig_NP = zeros(length(pop),1).*nan; 
HD_sig = zeros(length(pop),1).*nan; RH_sig = zeros(length(pop),1).*nan; 
for i = 1:length(fp)
    pop_idx = unitNow(i);
    % grab the modulation strengths for each unit
    MS_HD_now = pop(pop_idx).out.measures.TS.HD;
    MS_RH_now = pop(pop_idx).out.measures.TS.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = fp(i).B.mshd;
    RH_shuf_now = fp(i).B.msrh;
    
    % find confidence interval for the shuffled distribution
    HD_ci = [mean(HD_shuf_now) - 2*std(HD_shuf_now), mean(HD_shuf_now) + 2*std(HD_shuf_now)];
    RH_ci = [mean(RH_shuf_now) - 2*std(RH_shuf_now), mean(RH_shuf_now) + 2*std(RH_shuf_now)];
    
    HD_sig(pop_idx) = MS_HD_now < HD_ci(1) | MS_HD_now > HD_ci(2);
    RH_sig(pop_idx) = MS_RH_now < RH_ci(1) | MS_RH_now > RH_ci(2);
    
    % sort the shuffled distribution and real values
    HD_sort = sort([HD_shuf_now, MS_HD_now], 'ascend');
    RH_sort = sort([RH_shuf_now, MS_RH_now], 'ascend');
    
    % use a non-parametric method to determine significance
    if find(HD_sort == MS_HD_now) > 950 | find(HD_sort == MS_HD_now) < 50; 
        HD_sig_NP(pop_idx) = 1; 
    else
        HD_sig_NP(pop_idx) = 0;
    end
    
    if find(RH_sort == MS_RH_now) > 950 | find(RH_sort == MS_RH_now) < 50; 
        RH_sig_NP(pop_idx) = 1;
    else
        RH_sig_NP(pop_idx) = 0;
    end
end

disp('COMPUTING MODEL SENSITIVITY MEASURES...')
disp(['HD TS (parametric):', num2str(nansum(HD_sig)), '%']);
disp(['HD TS (non-parametric):', num2str(nansum(HD_sig_NP)), '%']);
disp(['RH TS (parametric):', num2str(nansum(RH_sig)), '%']);
disp(['RH TS (non-parametric):', num2str(nansum(RH_sig_NP)), '%']);

sig_idx = find(RH_sig_NP==1);
nsig_idx = find(RH_sig_NP==0);


%% # POINTS IN ENVIRONMENT
percentEgo = (length(find(xref(sig_idx)<150&xref(sig_idx)>0&yref(sig_idx)<150&yref(sig_idx)>0))...
    ./length(RH_sig_NP)).*100;

sig_in = find(xref(sig_idx)<150 & xref(sig_idx)>0 & yref(sig_idx)<150 & yref(sig_idx)>0);
sig_in = sig_idx(sig_in);
sig_out = 1:1000;
ix = union(sig_in, nsig_idx);
sig_out(ix) = nan;
sig_out = sig_out(find(~isnan(sig_out)));