%% run the RH model on Jan Sigurd's dataset (via IDUN)

addpath(genpath('/cluster/home/jordannc'))
load('/cluster/home/jordannc/data/jansig/dset.mat')

for sess = 1:86
    P = dset(sess).P; t = P(:,1);
    Z = get_hd(P);
    for u = 1:length(dset(sess).unit)
        disp(['session ', num2str(sess), ' unit ', num2str(u), '...'])
        ST = dset(sess).unit(u).ST;
        
        % fit the model
        totalruns = 100;
        for run = 1:totalruns
            out = modelMe(P, ST, Z);
            monte(run).model = out; 
            err4compare(run) = out.model.error;
        end
        % find run that yieled smallest error value
        [~, errorValsMin] = nanmin(err4compare);
        % save the run that gives the global min
        modelData = monte(errorValsMin).model;
        
        % false positive
        total_shifts = 1000;
        for shift_iter = 1:total_shifts
            disp(['iter: ', num2str(iter), ', shift: ', num2str(shift_iter)]);
            % grab the current shuffle value from vector
            shift_val = floor(60/tpf + (600/tpf-60/tpf).*rand(1));
            % circularly shift the elements in array A by K positions
            Z_shift = circshift(Z, shift_val);
            % run the model on the shuffled data
            [out_shift] = modelMe(P, ST, Z_shift);
            % SAVE OUTPUT FROM EACH RUN        
            mshd_shift(shift_iter) = out_shift.measures.TS.HD;
            msrh_shift(shift_iter) = out_shift.measures.TS.RH;
        end
        
        % sort the shuffled distribution and real values
        RH_sort = sort([msrh_shift, modelData.measures.TS.RH], 'ascend');
        
        if find(RH_sort == modelData.measures.TS.RH) > 950
            RHsig = 1; 
        else
            RHsig = 0;
        end
        
        % package output
        dsetRH(sess).unit(u).out = modelData;
        dsetRH(sess).unit(u).sig = RHsig;
        
        if mod(sess,10) == 0
            disp('saving dsetRH.mat...');
            save('/cluster/home/jordannc/output/jansig/dsetRH.mat', 'dsetRH', '-v7.3');
        end
    end
end

disp('final save...')
save('/cluster/home/jordannc/output/jansig/dsetRH.mat', 'dsetRH', '-v7.3');

