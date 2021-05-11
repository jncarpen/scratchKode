%% analyze results from the RH data for Jan Sigurd's data

% pull out pre-calculated information*
clear vep sigRH verh mshd msrh g theta xref yref dis
count = 1;
 for sess = 1:length(dsetRH)
     for u = 1:length(dsetRH(sess).unit)
         vep(count) = dsetRH(sess).unit(u).out.measures.VE.place;
         sigRH(count) = dsetRH(sess).unit(u).sigRH;
         verh(count) = dsetRH(sess).unit(u).out.measures.VE.RH;
         mshd(count) = dsetRH(sess).unit(u).out.measures.TS.HD;
         msrh(count) = dsetRH(sess).unit(u).out.measures.TS.RH;
         g(count) = dsetRH(sess).unit(u).out.model.fitParams.g;
         theta(count) = dsetRH(sess).unit(u).out.model.fitParams.thetaP;
         xref(count) = dsetRH(sess).unit(u).out.model.fitParams.xref.*15;
         yref(count) = dsetRH(sess).unit(u).out.model.fitParams.yref.*15;
         dis(count) = sqrt(((xref(count))-75).^2+((yref(count))-75).^2);
         sessNow(count) = sess;
         unitNow(count) = u;
         siContentSig(count) = dsetSI(sess).unit(u).contSig;
         siRateSig(count) = dsetSI(sess).unit(u).rateSig;
         siRate(count) = dsetSI(sess).unit(u).irate;
         siContent(count) = dsetSI(sess).unit(u).icontent;
         count = count + 1;
         
     end
 end
 sig_idx = find(sigRH==1);
 nsig_idx = find(sigRH==0);
 
 
%  
%% PLOT 1D HISTOGRAMS
thing = rbar;
thing1 = thing(sig_in);
thing2 = thing(nsig_idx);
thing0 = thing(sig_out);
nBins = 30;
binrng = min(thing):(max(thing)-min(thing))/nBins:max(thing);  

counts_shuffle = histc(rbar_shuff, binrng); hold on;
b3= bar(binrng, counts_shuffle, 'FaceColor','r');
b3.FaceAlpha =.5;

counts0 = histc(thing0, binrng);
counts1 = histc(thing1, binrng);
counts2 = histc(thing2, binrng);
counts3 = counts0 + counts1 + counts2;
figure(1)
b1 = bar(binrng, counts3, 'k');
hold on
b3 = bar(binrng, counts0, 'FaceColor', [.6 .6 .6]);
b2 = bar(binrng, counts1, 'b');
legend('Data 1', 'Data 2');
set(gca, 'FontSize', 15); box off;
l = legend({'shuff', 'not tuned', 'tuned (allo)', 'tuned (ego)'});
% l.Location = 'northeastoutside';
xlabel('r (model)');
ylabel('# neurons');



% %% PLOT 2D HISTOGRAMS
customColor = flipud(hot);
customColor = customColor(10:length(customColor)-10,:);
nBins = 20;
[N, xEdges, yEdges, binX, binY] = histcounts2(ctrMassX,ctrMassY,nBins);
N(N==0)=nan;
imagescwithnan(N,customColor,[1 1 1]); colorbar
set(gca, 'FontSize', 15); box off;     
%  
%  
%  %% TOTAL # SIGNIFICANT UNITS
%  bigCount = 0;
%  totalCount = 1;
%  for sess = 1:length(dsetRH)
%      smallCount = 0;
%      for u = 1:length(dsetRH(sess).unit)
%          signow = dsetRH(sess).unit(u).sigRH;
%          totalCount = totalCount + 1;
%          if signow == 1
%              bigCount = bigCount + 1;
%              smallCount = smallCount + 1;
%          end
%      end
%      cntPerSess(sess) = smallCount;
%      totalPerSess(sess) = length(dsetRH(sess).unit);
%  end
%  cntPerSess = cntPerSess';
%  totalPerSess = totalPerSess';
%  percentSess = 100.*(cntPerSess./totalPerSess);