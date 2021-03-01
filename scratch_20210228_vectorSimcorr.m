%% FEB 28
dir = out.measures.mu.HD;
sf = out.measures.MVL.HD;
sf(isnan(sf))=0;
xpos = out.info.bin.X;
ypos = out.info.bin.Y;

% directional components
uraw = cos(dir * pi/180); 
vraw = sin(dir * pi/180);
[sz_i, sz_j] = size(uraw);

% scaled vectors
u = uraw.*sf; 
v = vraw.*sf;

% linearize
ulin = u(:); vlin = v(:);
xlin = xpos(:); ylin = ypos(:);

clear corrmat
for a = 1:length(u(:))
    for b = 1:length(u(:))
        u1 = ulin(a); u2 = ulin(b);
        v1 = vlin(a); v2 = vlin(b);
        x1 = xlin(a); x2 = xlin(b);
        y1 = ylin(a); y2 = ylin(b);
        [sim_pipj] = vsim(x1, x2, y1, y2, u1, u2, v1, v2);
        corrmat(a,b) = sim_pipj;
        
    end
end

a = 5; b = 5;
for distin = 1:10
    for distout = 1:10
        binsInX = find(abs([1:10]-a)<distin);
        binsInY = find(abs([1:10]-b)<distin);
        for i = 1:length(binsInX)
            sim = vsim(u(binsInX(i)),u(binsInX(j),
        
        binsOutX = find(abs([1:10]-a)>=distin);
        binsOutY = find(abs([1:10]-b)>=distin);
        
    end
end




