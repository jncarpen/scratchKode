%% Bimodality statistic
% load data
pathname = 'D:\Data\Project Data\Simulation-POP2\place100.mat';
load(pathname);

%% data distribution
% modulation strength distribution
for i = 1:100
    ms(i) = place100(i).out.measures.TS.RH;
end

% histogram (each bin should have >5 counts)
nbins = 10;
totalcnt = length(ms);
[bincnt,edges,binnum] = histcounts(ms, linspace(0,1,nbins+1));
% bin lower, center, upper limits
xl = edges(1:end-1); xu = edges(1:end-1);
xc = (diff(edges)/2) + edges(1:end-1);
% probability distribution
prop = (bincnt./totalcnt).*100;

% set mean, sd, c
mean = .5; sd = 1; c = 0;

%% predicted unimodal distribution
xa = abs(xl-xc);
xb = abs(xu-xc);
h1 = (.398942/sd)*exp(-(((xl-mean).^2)/(2*sd.^2)));
h2 = (.398942/sd)*exp(-(((xc-mean).^2)/(2*sd.^2)));
h3 = (.398942/sd)*exp(-(((xu-mean).^2)/(2*sd.^2)));
preduni = (.5.*(h1+h2).*xa+.5.*(h2+h3).*xb)+ c;


%% predicted bimodal distribution
mean1=.2; mean2=.80; sd1=.015; sd2=.09; ratio=.3; c=0;
h1 = (.398942/sd1)*exp(-(((xl-mean1).^2)./(2*sd1.^2)));
h2 = (.398942/sd1)*exp(-(((xc-mean1).^2)./(2*sd1.^2)));
h3 = (.398942/sd1)*exp(-(((xu-mean1).^2)./(2*sd1.^2)));
h4 = (.398942/sd2)*exp(-(((xl-mean2).^2)./(2*sd2.^2)));
h5 = (.398942/sd2)*exp(-(((xc-mean2).^2)./(2*sd2.^2)));
h6 = (.398942/sd2)*exp(-(((xu-mean2).^2)./(2*sd2.^2)));
predbi = ratio .*(.5*(h1+h2).*xa + .5*(h2+h3).*xb)...
    +(1-ratio).*(.5*(h4+h5).*xa + .5*(h5+h6).*xb) + c;
plot(xc,predbi); hold on; plot(xc, prop./100);








