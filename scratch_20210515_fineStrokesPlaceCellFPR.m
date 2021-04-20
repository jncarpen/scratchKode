
sesslist = [1:length(OF), floor(linspace(4,84,14))];


%% MODEL (RH)
for i = 1:100
    sigx(i) = place100(i).param.sigma(1);
    sigy(i) = place100(i).param.sigma(2);
    ctrx(i) = place100(i).param.ctr(1);
    ctry(i) = place100(i).param.ctr(2);
    A(i) = place100(i).param.A;
    sessNum(i) = sesslist(i);
    rhFPHere(i) = rhFP(i);
    hdFPHere(i) = hdFP(i);
    numSpk(i) = length(place100(i).ST);
    
    g_fit(i) = place100(i).out.model.fitParams.g;
    theta_fit(i) = place100(i).out.model.fitParams.thetaP;
    x_fit(i) = place100(i).out.model.fitParams.xref;
    y_fit(i) = place100(i).out.model.fitParams.yref;
    d(i) = sqrt(((75-x_fit(i)).^2) + ((75-y_fit(i)).^2));
    dfield(i) = sqrt(((ctrx(i)-x_fit(i).*15).^2) + ((ctry(i)-y_fit(i).*15).^2));
    
    vep(i) = place100(i).out.measures.VE.place;
    verh(i) = place100(i).out.measures.VE.RH;
    ms_rh(i) = place100(i).out.measures.TS.RH;
    ms_hd(i) = place100(i).out.measures.TS.HD;
end


thing = dfield.*15;
thing1 = thing(rhsig_idx);
thing2 = thing(nsig_idx);
histogram(thing2, 10,'FaceColor', 'k'); hold on;
histogram(thing1, 10,'FaceColor', 'r');
l = legend({'not significant', 'significant'});
set(gca, 'FontSize', 15); box off;



% plot
colorStuff = zeros(100,1);
colorStuff(rhsig_idx) = 1;
s=scatter(x_fit.*15, ctrx, [40], colorStuff, 'filled');
s.MarkerFaceAlpha = .75;
colormap('copper')
set(gca, 'FontSize', 15); box off;
title('Place Cells (Simulated)')
xlabel('x-fit'); 
ylabel('y-fit');







