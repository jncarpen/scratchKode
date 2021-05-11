ctrX = 75; ctrY = 75;
for i = 1:length(xref)
    x = xref(i); y = yref(i);
    if x>300 | x<-300 | y>300 | y<-300
        if x<-300 & y>-300 & y<300 % x=-300 (vert)
            xLine = -300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            yLine = (m*xLine)+b;
        elseif x>-300 & x<300 & y<-300
            yLine = -300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            xLine = (yLine-b)/m;
        elseif x>300 & y>-300 & y<300
            xLine = 300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            yLine = (m*xLine)+b;
        elseif x>-300 & x<300 & y >300
            yLine = 300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            xLine = (yLine-b)/m;
        else 
            xLine = x; yLine = y;
        end
        xrefAdj(i) = xLine; yrefAdj(i) = yLine;
    end
end