%% plot GLM results for example units
clear all;
loadpath = 'D:\Data\Project Data\GLM-exampleunits\glm-egodist270.mat';
load(loadpath);

%% plot real tuning curves
fig = figure('units','normalized','outerposition',[0 .5 1 .5]);
set(gcf,'color','w');
r = 1; c = 6; maxdist = 120;
lw = 1;

% plot the tuning curves (real data)
figure(1); set(gcf,'color','w');
subplot(r,c,1)
imagesc(glmout.curve.pos); colorbar
axis off; pbaspect([1 1 1]);
title('position (data)', 'FontWeight', 'normal')
set(gca,'FontSize',15);

subplot(r,c,2)
imagesc(reshape(glmout.response.pos,20,20)); colorbar
axis off; pbaspect([1 1 1]);
title('position (model)', 'FontWeight', 'normal')
set(gca,'FontSize',15);

subplot(r,c,3)
plot(rad2deg(glmout.xvec.hd),glmout.curve.hd,'k','linewidth',lw); hold on;
plot(rad2deg(glmout.xvec.hd),glmout.response.hd,'b','linewidth',lw);
box off; pbaspect([1 1 1]);
axis([0 360 -inf inf]); xticks([0 90 180 270]);
xlabel('angle (deg)')
title('head direction', 'FontWeight', 'normal')
set(gca,'FontSize',15);

subplot(r,c,4)
plot(glmout.xvec.speed,glmout.curve.speed,'k','linewidth',lw); hold on;
plot(glmout.xvec.speed,glmout.response.speed,'b','linewidth',lw); hold on;
box off; pbaspect([1 1 1]);
xlabel('running speed')
axis([0 50 -inf inf])
title('speed', 'FontWeight', 'normal')
set(gca,'FontSize',15);

subplot(r,c,5)
plot(rad2deg(glmout.xvec.hd),glmout.curve.ego,'k','linewidth',lw); hold on;
plot(rad2deg(glmout.xvec.hd),glmout.response.ego,'b','linewidth',lw);
xlabel('angle (deg)')
axis([0 360 -inf inf]); xticks([0 90 180 270]);
box off; pbaspect([1 1 1]);
title('egobearing', 'FontWeight', 'normal')
set(gca,'FontSize',15);

subplot(r,c,6)
plot(glmout.xvec.dist,glmout.curve.dist,'k','linewidth',lw); hold on;
plot(glmout.xvec.dist,glmout.response.dist,'b','linewidth',lw);
xlabel('dist (cm)')
axis([0 maxdist -inf inf])
box off; pbaspect([1 1 1]);
title('distance to refpoint', 'FontWeight', 'normal')
set(gca,'FontSize',15);
