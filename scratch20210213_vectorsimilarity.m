%% February 13, 2020
% Compute a similarity matrix for a vector field
% Similarity measure from publication:
% Similarity Measure for Vector Field Learning (Hongyu Li,I-Fan Shen 2006)

%% calculate vector field
dir = out.measures.mu.HD;
scale = out.measures.MVL.HD;
scale(isnan(scale))=0;
xpos = out.info.bin.X;
ypos = out.info.bin.Y;

% directional components
uraw = cos(dir * pi/180); 
vraw = sin(dir * pi/180);
[sz_i, sz_j] = size(uraw);

% scaled vectors
u = uraw.*sf; 
v = vraw.*sf;

%% 1. alpha: euclidean distance
% calculate [Euclidean] distance term
dist_term = exp(-sqrt(((xpos(5,5)-xpos(5,6)).^2)+((ypos(5,5)-ypos(5,6)).^2)));

%% 2. beta: angular difference
vi = [u(5,5) v(5,5)];
vj = [u(5,6) v(5,6)];

% take dot product of vi.vj
dot_vivj = dot(vi,vj);

% calculate the norm of vi and vj
norm_vi = norm(vi);
norm_vj = norm(vj);

% calcuate angular term
ang_term = exp(1-(dot_vivj./(norm_vi * norm_vj)));

%% 3. gamma: magnitude difference
% calculate magnitude term
mag_term = exp(-(norm(norm_vi-norm_vj)));

%% 4. vector pair similarity
% impose req. that all components contribute equally
alpha = 1/3;
beta = 1/3;
gamma = 1/3;

% calculate similarity value (range 0 to 1)
sim_pipj = alpha*dist_term + beta*ang_term + gamma*mag_term;

%% sliding window approach

% x contains values to be filtered
x = xpos; y = ypos;
% set horizontal and vertical limits
hl = 2; vl = 2;
% check size of input
[row col space]=size(x);
% make an output vector of same size and type as input
sim_mean=0*x; 
sim_std=0*x; 


for i=1:row % loop thru rows
    for j=1:col % loop thru values one by one
        if i==1 & j==1 % first row, first col
            sim = zeros(3,1);
            sim(1) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
            sim(2) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
            sim(3) = vsim(x(i,j), x(i+1,j+1), y(i,j), y(i+1,j+1), u(i,j), u(i+1,j+1), v(i,j), v(i+1,j+1));
        elseif i==1 & j==col % first row, last col
            sim = zeros(3,1);
            sim(1) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
            sim(2) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
            sim(3) = vsim(x(i,j), x(i+1,j-1), y(i,j), y(i+1,j-1), u(i,j), u(i+1,j-1), v(i,j), v(i+1,j-1));
        elseif i==row & j==1 % last row, first col
            sim = zeros(3,1);
            sim(1) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
            sim(2) = vsim(x(i,j), x(i-1,j+1), y(i,j), y(i-1,j+1), u(i,j), u(i-1,j+1), v(i,j), v(i-1,j+1));
            sim(3) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
        elseif i==row & j==col % last row, last col
            sim = zeros(3,1);
            sim(1) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
            sim(2) = vsim(x(i,j), x(i-1,j-1), y(i,j), y(i-1,j-1), u(i,j), u(i-1,j-1), v(i,j), v(i-1,j-1));
            sim(3) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
        elseif i==1 & j~=1 & j~=col % first row, not first or last col
            sim = zeros(5,1);
            sim(1) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
            sim(2) = vsim(x(i,j), x(i+1,j+1), y(i,j), y(i+1,j+1), u(i,j), u(i+1,j+1), v(i,j), v(i+1,j+1));
            sim(3) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
            sim(4) = vsim(x(i,j), x(i+1,j+1), y(i,j), y(i+1,j+1), u(i,j), u(i+1,j+1), v(i,j), v(i+1,j+1));
            sim(5) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
        elseif i==row & j~=1 & j~=col % last row, not first or last col
            sim = zeros(5,1);
            sim(1) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
            sim(2) = vsim(x(i,j), x(i-1,j-1), y(i,j), y(i-1,j-1), u(i,j), u(i-1,j-1), v(i,j), v(i-1,j-1));
            sim(3) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
            sim(4) = vsim(x(i,j), x(i-1,j+1), y(i,j), y(i-1,j+1), u(i,j), u(i-1,j+1), v(i,j), v(i-1,j+1));
            sim(5) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
        elseif j==1 & i~=1 & i~=row % not first or last row, first col
            sim = zeros(5,1);
            sim(1) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
            sim(2) = vsim(x(i,j), x(i-1,j+1), y(i,j), y(i-1,j+1), u(i,j), u(i-1,j+1), v(i,j), v(i-1,j+1));
            sim(3) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
            sim(4) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
            sim(5) = vsim(x(i,j), x(i+1,j+1), y(i,j), y(i+1,j+1), u(i,j), u(i+1,j+1), v(i,j), v(i+1,j+1));
        elseif j==col & i~=1 & i~=row % not first or last row, last col
            sim = zeros(5,1);
            sim(1) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
            sim(2) = vsim(x(i,j), x(i-1,j-1), y(i,j), y(i-1,j-1), u(i,j), u(i-1,j-1), v(i,j), v(i-1,j-1));
            sim(3) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
            sim(4) = vsim(x(i,j), x(i+1,j-1), y(i,j), y(i+1,j-1), u(i,j), u(i+1,j-1), v(i,j), v(i+1,j-1));
            sim(5) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
        else
            sim = zeros(8,1);
            sim(1) = vsim(x(i,j), x(i-1,j-1), y(i,j), y(i-1,j-1), u(i,j), u(i-1,j-1), v(i,j), v(i-1,j-1));
            sim(2) = vsim(x(i,j), x(i-1,j), y(i,j), y(i-1,j), u(i,j), u(i-1,j), v(i,j), v(i-1,j));
            sim(3) = vsim(x(i,j), x(i-1,j+1), y(i,j), y(i-1,j+1), u(i,j), u(i-1,j+1), v(i,j), v(i-1,j+1));
            sim(4) = vsim(x(i,j), x(i,j+1), y(i,j), y(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1));
            sim(5) = vsim(x(i,j), x(i+1,j+1), y(i,j), y(i+1,j+1), u(i,j), u(i+1,j+1), v(i,j), v(i+1,j+1));
            sim(6) = vsim(x(i,j), x(i+1,j), y(i,j), y(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j));
            sim(7) = vsim(x(i,j), x(i+1,j-1), y(i,j), y(i+1,j-1), u(i,j), u(i+1,j-1), v(i,j), v(i+1,j-1));
            sim(8) = vsim(x(i,j), x(i,j-1), y(i,j), y(i,j-1), u(i,j), u(i,j-1), v(i,j), v(i,j-1));
        end
        nanmax(sim)
        sim_std(i,j) = std(sim, 'omitnan');
        sim_mean(i,j) = mean(sim, 'omitnan');
    end
end









