%% FALSE POSITIVE RATE (WITH SIMULATED CELLS)
%   January 25, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = ego100;
warning('off','all')
for simcell = 1:100
    % tell user which iteration youre on
    dispMe = strcat('running simulation ', {' '}, sprintf('%.f', simcell), ' of 100...');
    disp(dispMe{1,1});
    % information about current cell
    P = A(simcell).param.P; tpf = mode(diff(P(:,1)));
    ST = A(simcell).ST; Z = A(simcell).param.Z;
    
    %% shuffle the head direction values
    total_shifts = 1000;
    for shift_iter = 1:total_shifts
        % grab the current shuffle value from vector
        shift_val = floor(60/tpf + (600/tpf-60/tpf).*rand(1));
        % circularly shift the elements in array A by K positions
        Z_shift = circshift(Z, shift_val);
        % run the model on the shuffled data
        [out_shift] = modelMe(P, ST, Z_shift);
        % SAVE OUTPUT FROM EACH RUN        
        % (1) best fit parameters
        g_shift(shift_iter) = out_shift.model.fitParams.g;
        theta_shift(shift_iter) = out_shift.model.fitParams.thetaP;
        x_shift(shift_iter) = out_shift.model.fitParams.xref;
        y_shift(shift_iter) = out_shift.model.fitParams.yref;
        % (2) variance explained (place & model)
        vep_shift(shift_iter) = out_shift.measures.VE.place;
        verh_shift(shift_iter) = out_shift.measures.VE.RH;
        % (3) modulation/tuning strength
        mshd_shift(shift_iter) = out_shift.measures.TS.HD;
        msrh_shift(shift_iter) = out_shift.measures.TS.RH;
        % (4) error
        err_shift(shift_iter) = out_shift.model.error;
    end
    
    % save stuff
    B(simcell).vep = vep_shift;
    B(simcell).verh = verh_shift;
    B(simcell).mshd = mshd_shift;
    B(simcell).msrh = msrh_shift;
    B(simcell).err = err_shift;
    B(simcell).g = g_shift;
    B(simcell).theta = theta_shift;
    B(simcell).x = x_shift;
    B(simcell).y = y_shift;
    
    % save every 10th iteration
    if mod(simcell,10) == 0
        save('D:\Data\External Data\Simulation-FP\ego100FP.mat', 'B', '-v7.3');
    end
end
warning('on','all')