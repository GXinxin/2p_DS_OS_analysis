function C = myDrawCircle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
C = plot(xunit, yunit, 'color', [129 129 129]/255, 'lineWidth', 2);
