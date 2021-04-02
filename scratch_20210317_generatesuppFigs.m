    pop = ego100;
    % fp = ego100FP;

    % pull out model fit parameters
    for i = 1:100
%         g(i) = pop(i).out.model.fitParams.g;
%         xref(i) = pop(i).out.model.fitParams.xref*15;
%         yref(i) = pop(i).out.model.fitParams.yref*15;
        theta(i) = mod(pop(i).out.model.fitParams.thetaP,360);
        % ctrx(i) = pop(i).param.ctr(1);
        % ctry(i) = pop(i).param.ctr(2);

        param_theta(i) = mod(pop(i).param.theta,360);
    %     param_x(i) = pop(i).param.rp(1);
    %     param_y(i) = pop(i).param.rp(2);
    %     
%         tsrh(i) = pop(i).out.measures.TS.RH;
    %     tshd(i) = pop(i).out.measures.TS.HD;

%         verh(i) = pop(i).out.measures.VE.RH;
    %     vep(i) = pop(i).out.measures.VE.place;

    %     verh_shuff(i) = mean(fp(i).B.verh, 'omitnan');
    %     vep_shuff(i) = mean(fp(i).B.vep, 'omitnan');
    %     
    %     tsrh_shuff(i) = mean(fp(i).B.msrh, 'omitnan');
    %     tshd_shuff(i) = mean(fp(i).B.mshd, 'omitnan');
    end
    
    clear g i pop theta xref yref


% get distance between refpoint & predicted (in cm)
d = sqrt((param_x-xref).^2 + (param_y-yref).^2);
dt = abs(param_theta-theta); 

% g parameter (histogram)
histogram(g, 20, 'FaceColor', 'k');
xlabel('model-predicted g (amplitude)');
ylabel('units (%)');
set(gca,'FontSize', 15); box off;

% theta parameter (histogram)
histogram(theta, 20, 'FaceColor', 'k');
xlabel('model-predicted theta');
ylabel('units (%)');
set(gca,'FontSize', 15); box off;

% place field ctr v. predicted
close all;
s = scatter(ctrx, xref, [30], 'k', 'filled');
s.MarkerFaceAlpha = .4;
xlabel('place field center (x)');
ylabel('model-predicted refpoint (x)');
set(gca,'FontSize', 15); box off;


%% plot all units

for i = 1:100
    venow = pop(i).out.measures.VE.RH*100;
    if venow > 55
        subplot(1,2,1)
        pathPlot_hd(pop(i).param.P, pop(i).ST, pop(i).param.Z);
        subplot(1,2,2)
        plotMe(pop(i).out); hold on;
        scatter(pop(i).out.model.fitParams.xref, pop(i).out.model.fitParams.yref, ...
            [30], 'k', 'filled')
        title([num2str(i), ', VE: ',num2str(pop(i).out.measures.VE.RH*100)]);
        pause; clf;
    end
end

close all;
i = 98;
figure(1)
pathPlot_hd(pop(i).param.P, pop(i).ST, pop(i).param.Z);
scatter(pop(i).param.rp(1), pop(i).param.rp(2), [40], 'k', 'filled');

figure(2)
plotMe(pop(i).out); hold on;
scatter(pop(i).out.model.fitParams.xref, pop(i).out.model.fitParams.yref, ...
            [30], 'k', 'filled')

disp(['VE: ', num2str(pop(i).out.measures.VE.RH*100), '%'])
disp(['MS: ', num2str(pop(i).out.measures.TS.RH)])
disp(['theta: ', num2str(pop(i).param.theta)])

% plot the distributions
x = MS_ego';
vals = histcounts(x, 15);
pd = fitdist(vals','normal');
y = pdf(pd,vals);
plot(vals, y,'LineWidth',2)

figure;
hist{1,1} = histfit(theta_ego, 10, 'kernel');
figure;
hist{1,2} = histfit(theta_ring, 10, 'kernel');
figure;
hist{1,3} = histfit(theta_pe, 10, 'kernel');
figure;
hist{1,4} = histfit(theta_place, 10, 'kernel'); title('p')
figure;
hist{1,5} = histfit(theta_hd, 10, 'kernel');
figure;
hist{1,6} = histfit(theta_ph, 10, 'kernel');

for i = 1:6
    xd{1,i} = hist{1,i}(2).XData;
    yd{1,i} = hist{1,i}(2).YData;
end

figure; hold on;
for i = 1:6
    plot(xd{1,i}, yd{1,i}/100, 'LineWidth', 1);
end
set(gca,'FontSize', 15); box off;
xlabel('variance explained');
ylabel('pdx')
legend({'E', 'ED', 'PE', 'P', 'HD', 'PHD'})


savepath = 'D:\Data\Project Data\Simulation-POP2\g.mat';
save(savepath, 'g_ego', 'g_ring', 'g_pe', 'g_place', 'g_hd', 'g_ph', '-v7.3')



%%
x_ego(x_ego>450) = 500; x_ego(x_ego<-300) = -350;
y_ego(y_ego>450) = 500; y_ego(y_ego<-300) = -350;

x_ring(x_ring>450) = 500; x_ring(x_ring<-300) = -350;
y_ring(y_ring>450) = 500; y_ring(y_ring<-300) = -350;

x_pe(x_pe>450) = 500; x_pe(x_pe<-300) = -350;
y_pe(y_pe>450) = 500; y_pe(y_pe<-300) = -350;

x_place(x_place>450) = 500; x_place(x_place<-300) = -350;
y_place(y_place>450) = 500; y_place(y_place<-300) = -350;

x_hd(x_hd>450) = 500; x_hd(x_hd<-300) = -350;
y_hd(y_hd>450) = 500; y_hd(y_hd<-300) = -350;

x_ph(x_ph>450) = 500; x_ph(x_ph<-300) = -350;
y_ph(y_ph>450) = 500; y_ph(y_ph<-300) = -350;


figure; hold on;
s1 = scatter(x_ego, y_ego, [5], 'k', 'filled');
s2 = scatter(x_ring, y_ring, [5], [0 .8 .8], 'filled');
s3 = scatter(x_pe, y_pe, [5], 'o', 'filled');
s4 = scatter(x_place, y_place, [5], 'b', 'filled');
s5 = scatter(x_hd, y_hd, [5], [0 .8 0],'filled');
s6 = scatter(x_ph, y_ph, [5], 'r', 'filled');
xline(0); xline(150); yline(0); yline(150);
xlabel('model-predicted refpoint (x)');
ylabel('model-predicted refpoint (y)');
set(gca,'FontSize', 15); box off;
xlim([-300 450]); ylim([-300 450]);


