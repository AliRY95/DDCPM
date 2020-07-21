function PLOT(Idx, XGrid, YGrid, xxx)
theta = 0:0.01:2*pi;
x = cos(theta); y = sin(theta);
plot(XGrid(Idx), YGrid(Idx), xxx); axis equal; hold on;
plot(x, y, 'b')
end