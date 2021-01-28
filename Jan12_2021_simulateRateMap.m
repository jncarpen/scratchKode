%% Simulate ratemap, scratch code.
% J. Carpenter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a ratemap
i = 10;
xpos = optiData{1,i}.position.x;
ypos = optiData{1,i}.position.y;

% pull out mu (place field center of mass)
xmu = -0.5;
ymu = -0.33;

% generate 2d gaussian kernel
sigma = 0.15; nBins = 50;
x = linspace(min(xpos),max(xpos),nBins);
[Y,X] = meshgrid(x,x);
f = exp(-((X-xmu).^2 + (Y-ymu).^2)./(2*(sigma^2)));
f = f / sum(f(:));
g = (1/(2*pi*(sigma^2))).* f;

% plot
figure; set(gcf,'color','w');
subplot(1,2,1)
surf(f)
pbaspect([1 1 1])

subplot(1,2,2)
peakRate = nanmax(nanmax(f));
rate_map_title = strcat('max: ', sprintf('%.2f',peakRate));
plot.colorMap(f);
pbaspect([1 1 1])
c = colorbar; c.FontName = 'Helvetica'; c.FontSize = 15;
colormap(gca,'jet')
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Helvetica', 'FontSize', 15, 'FontWeight', 'normal');
box off

clear test
test = zeros(nBins,nBins).*NaN;
for XX = 1:nBins
    for YY = 1:nBins
        lambda = f(XX,YY)*10000;
        
        test(XX,YY) = poissrnd(lambda);
        
    end
end
imagesc(test)
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write a function to simulate a gaussian location preference. this is for 
% a sum of gaussians, so we can create bimodal distributions if we want.
def sum_of_gaussians(a, positions, c):
    f = lambda x: sum(a * np.exp( -(x-b)**2 / (2 * c**2)) for b in positions)
    return f
    
% sum of three gaussians
amp0 = 13721;
amp1=13721;
amp2=14753;




x_1 = 2
x_2 = 2
y_1 = -1
y_2 = -1

o_x_1 = 1
o_y_1 = 2
o_x_2 = 2
o_y_2 = 1

%// Define grid of points
[X,Y] = meshgrid(-10:0.01:10, -10:0.01:10);

%// Define parameters for each Gaussian
x_1 = 2;
y_1 = 2;
x_2 = -1;
y_2 = -1;

o_x_1 = 1;
o_y_1 = 2;
o_x_2 = 2;
o_y_2 = 1;

%// Define constant A... let's just assume 1 for both
A = 1;

%// Generate Gaussian values for both Gaussians
f1 = A*exp( -( ((X - x_1).^2 / (2*o_x_1^2)) + ((Y - y_1).^2 / (2*o_y_1^2)) ) );
f2 = A*exp( -( ((X - x_2).^2 / (2*o_x_2^2)) + ((Y - y_2).^2 / (2*o_y_2^2)) ) );

%// Add them up
f = f1 + f2;

%// Show the results
mesh(X, Y, f);
%// Label the axes
xlabel('x');
ylabel('y');
zlabel('z');
view(-130,50); %// For a better view
colorbar; %// Add a colour bar for good measure







    
    
    
    
    