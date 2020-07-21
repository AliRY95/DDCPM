function [X, Y] = CLOSESTPOINT(XGrid, YGrid)
[theta, ~] = cart2pol(XGrid, YGrid);
[X, Y] = pol2cart(theta, 1);
end