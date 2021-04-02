pop = hd100;

% pull out predicted refpoints
for i = 1:100
    x(i) = pop(i).out.model.fitParams.xref*15;
    y(i) = pop(i).out.model.fitParams.yref*15;
end

% find local refpoints
for i = 1:100
    xx = x(i); yy = y(i);
    if xx<150 & xx>0 & yy<150 & yy>0
        local(i) = 1;
    else
        local(i) = 0;
    end
end