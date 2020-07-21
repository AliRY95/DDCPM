clearvars; clc;
addpath(fullfile(pwd, 'ClassDef'));
addpath(fullfile(pwd, 'CP Functions'));
addpath(fullfile(pwd, 'Plot Functions'));
%% CURVE & PARAMETERS
Dimension = 2; InterpolationOrder = 4;
dx = 0.05; dy = 0.05;
c = 1;

NoOverlap = 5;

f = @(x) cos(x);
ExactSolution = @(x) cos(x)/2;
%% BAND & GRID
lambda = sqrt((Dimension-1)*(InterpolationOrder+1)^2/4 + (InterpolationOrder+3)^2/4)*dx;
r2 = (1 - lambda)^2;
R2 = (1 + lambda)^2;
temp_x = -1 - 5*dx:dx:1 + 5*dx; Nx = length(temp_x);
temp_y = 1 + 5*dy:-dy:-1 - 5*dy; Ny = length(temp_y);
[XGrid, YGrid] = meshgrid(temp_x, temp_y);

GlobalMask = (XGrid.^2 + YGrid.^2 >= r2).*(XGrid.^2 + YGrid.^2 <= R2);
% LocalMasks = {GlobalMask};
% LocalMasks = {GlobalMask.*(XGrid >= 0), GlobalMask.*(XGrid < 0)};
LocalMasks = {GlobalMask.*(YGrid >= 0), GlobalMask.*(YGrid <= 0)};
% LocalMasks = {GlobalMask.*(YGrid + XGrid >= 0), GlobalMask.*(YGrid + XGrid < 0)};

clear lambda r2 R2 temp_x temp_y;
%% SETTING-UP DOMAINS & SUBDOMAINS
Domain = cDomain;
Domain = Domain.SETUP(XGrid, YGrid, GlobalMask, c, f);
NoSubdomains = length(LocalMasks);
Subdomain(NoSubdomains) = cSubDomain;
for i = 1:NoSubdomains
    Subdomain(i) = Subdomain(i).SETUP(XGrid, YGrid, GlobalMask, LocalMasks{i}, NoOverlap, InterpolationOrder, c);
end
% CHECK OVERLAP
%% Initialization
u0 = zeros(length(Domain.ActiveNodesIdx), 1);
%% EXACT SOLUTION
xx = [Domain.Node.CPx];
yy = [Domain.Node.CPy];
[theta, ~] = cart2pol(xx', yy');
uExact = ExactSolution(theta);
%% RAS METHOD
u = u0; temp = 1; counter = 0;
while (abs(norm(temp, inf)) >= 1e-12 && counter <= 500)
    counter = counter + 1
    temp = 0;
    r = Domain.RHSMatrix - Domain.OperatorMatrix*u;
    for i = 1:NoSubdomains
        Local_r = Subdomain(i).OverlappingExtensionOperator * r;
        Local_r(Subdomain(i).LocalActiveBCNodesIdx) = 0;
        LocalCorrection = Subdomain(i).LocalOperatorMatrix\Local_r;
        temp = temp + Subdomain(i).DisjointExtensionOperator'*LocalCorrection;
    end
    disp(norm(temp, inf));
    u = u + temp;
    
    uu = Domain.InterpolationMatrix * u;
    err(counter) = norm(abs(uu - uExact), inf);
end
%% RESULTS
figure(1);
semilogy(1:counter, err, 'b');
xlabel Iteration; ylabel ||u-u_{exact}||_{\infty};
% title("Unit Circle(N_s=2, N_o=5)");
hold on;
    
figure(2);
uu = Domain.InterpolationMatrix * u;
plot3(xx, yy, uu, 'k.', 'MarkerSize', 16);
xlabel x; ylabel y; zlabel u;
title("DDCPM Solution");
axis equal; grid on

figure(3);
plot3(xx, yy, uExact, 'k.', 'MarkerSize', 16);
xlabel x; ylabel y; zlabel u;
title("Exact Solution");
axis equal; grid on

figure(4);
PLOT(Domain.NodesIdx, XGrid, YGrid, 'k.');
title("Computational Domain");