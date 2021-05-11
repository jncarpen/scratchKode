data = dsetRH;

count = 1;
for sess = 1:length(data)
    P = data(sess).P;
    for u = 1:length(data(sess).unit)
        isSignificant = data(sess).unit(u).sigRH;
        isRH = data(sess).unit(u).out.measures.TS.RH;
        if isSignificant == 1 & isRH >.65
            ST = data(sess).unit(u).ST;
            pref_theta = mod(data(sess).unit(u).out.model.fitParams.thetaP,360);
            pref_rp = [data(sess).unit(u).out.model.fitParams.xref, ...
                data(sess).unit(u).out.model.fitParams.yref];
            [tc, binCtrs] = egoBearing(P, ST, pref_rp, pref_rp, 15, 'False', 'deg');
            whichTC{count} = tc;
            whichSess(count) = sess;
            whichUnit(count) = u;
            whichTheta(count) = pref_theta;
            whichRPx(count) = pref_rp(1);
            whichRPy(count) = pref_rp(2);
            count = count + 1;
        end
        
    end
end

% put the theta values in order
clear whichTC_sorted
[~, orderTheta_idx] = sort(whichTheta, 'ascend');
whichTC_sorted = whichTC(orderTheta_idx);

% put sorted tc values into a matrix
clear tcMatrix
for j = 1:length(whichTC_sorted)
    tcMatrix(j,:) = zscore(whichTC_sorted{j});
end