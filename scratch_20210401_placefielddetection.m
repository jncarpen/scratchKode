% add BNT to the path
addpath(genpath('C:\Users\17145\OneDrive - NTNU\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001'))
% choose a session & cell

count = 1;
for sess = 1:86
    P = dset(sess).P; P(:,2:end) = P(:,2:end).*100;
    for u = 1:length(dset(sess).unit)
        % grab position and spike data
        ST = dset(sess).unit(u).ST;
        % ratemap
        map = analyses.map(P, ST);
        % place field detection
        [fieldsMap, fields] = analyses.placefield(map, 'pos', P, ...
            'minBins', 15,'minPeak', 2);
        subplot(1,2,1)
        pathPlot_hd(P, ST, get_hd(P));
        subplot(1,2,2)
        imagesc(flipud(fieldsMap)); pbaspect([1 1 1]); colorbar;
        pause; clf;
        
        fieldsMapArray{1,count} = fieldsMap;
        fieldsArray{1,count} = fields;
        
        % border score
        borderScore = analyses.borderScore(map.z, fieldsMap, fields);
        bs(count) = borderScore;
        sessnow(count) = sess;
        unitnow(count) = u;
        count = count + 1;
    end
end