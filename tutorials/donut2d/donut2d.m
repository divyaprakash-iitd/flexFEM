clear; clc; close all;

%% Add path
libpath = '../../'; % Replace with the path where you have downloaded the source code
addpath(libpath);

%% Description: Solves a demo problem using FEM

%% Define the elastic constants
K   = 10000;
co  = 5000;
E   = co*(2*co+3*K)/(co+K);
nu  = K/(2*(co+K));

%% Geometry and FEM Mesh

% Import gmsh mesh
donut2d_mesh;

% Get the points
points = msh.POS;
[M,x,y] = deal(msh.TRIANGLES(:,1:3), points(:,1), points(:,2));

% From the gmsh script, annular_mesh.geo, we know that there are 
% 3 physical groups
ocTag = 1; % Tag-1: Outer circumference (oc)
icTag = 2; % Tag-2: Inner circumference (ic)
asTag = 3; % Tag-3: The annular surface (as)

ocLinesId = msh.LINES(:,end) == ocTag;
icLinesId = msh.LINES(:,end) == icTag;

ocPointsId = unique(msh.LINES(ocLinesId,1:2));
icPointsId = unique(msh.LINES(icLinesId,1:2));

nNodes = msh.nbNod;

% Create masks
ocMask = false(nNodes,1);
ocMask(ocPointsId) = true;

icMask = false(nNodes,1);
icMask(icPointsId) = true;

% Select nodes based on angle threshold
nangle = atand(y./x);
sangle = 45;

% Final mask to apply boundary forces
leftMask = x < 0 & abs(nangle) <= sangle & ocMask;
rightMask = x > 0 & abs(nangle) <= sangle & ocMask;

%% Forces at the nodes
FN      = zeros(nNodes,2);

%% Run time iterations
niter   = 50000;
dt      = 0.001;
mu      = 1000;
filenumber = 1;

%% Applied Force
Fboundary        = 500;          % Tip force
direction   = 1;            % 1: Lateral (x-direction); 2: Axial (y-direction)  

%% Plot controls
plotinterval = niter/100;

fig = figure();
fig.Position = [675, 549, 570, 413];
xlabel('x [m]')
ylabel('y [m]')
hold on

%Store simulation images in a directory
folder = 'images';

% Check if the folder exists, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end

%% Simulation loop
% Create festruct
festruct = create_festruct(M,x,y,FN,co,K);
iter    = 1;
while iter <= niter

    % Initialize forces at the nodes
    FApplied = zeros(nNodes,2);

    % Add boundary force
    if iter > niter/2
        Fboundary = 0;
    end
    FApplied(leftMask,direction) = -Fboundary;
    FApplied(rightMask,direction) = Fboundary;
    FN = FApplied;

    % Feed in force field to festruct
    festruct.FN = FN;

    % Calculate total forces
    festruct = calculate_forces(festruct);

    % Read updated values from festruct
    x = festruct.x(:,1);
    y = festruct.x(:,2);
    FN = festruct.FN;

    % Calculate the deformation
    for inode = 1:nNodes
        festruct.x(inode,1) = festruct.x(inode,1) + dt/mu * FN(inode,1);
        festruct.x(inode,2) = festruct.x(inode,2) + dt/mu * FN(inode,2);
    end

    % Plot data points
    if mod(iter,plotinterval) == 0 && iter ~= niter
        x = festruct.x(:,1);
        y = festruct.x(:,2);
        cla;
        triplot(M,x,y,'r','DisplayName','Deformed')
        h1 = quiver(x(leftMask),y(leftMask),FApplied(leftMask,1),FApplied(leftMask,2),'Color','k');
        quiver(x(rightMask),y(rightMask),FApplied(rightMask,1),FApplied(rightMask,2), 'Color','k')
        h2 = quiver(x(leftMask),y(leftMask),FN(leftMask,1),FN(leftMask,2), 'Color','b');
        quiver(x(rightMask),y(rightMask),FN(rightMask,1),FN(rightMask,2), 'Color','b')

        daspect([1,1,1])
        xlim([-1.7,1.7])
        ylim([-1,1])
        axis equal
        title(sprintf('Iteration = %d',iter))
        lgd = legend([h1, h2], 'F_{Applied}', 'F_{Total}');
        lgd.Location = 'southeast';

        filename = sprintf('%s/image_%02d.png', folder, filenumber);
        exportgraphics(gcf,filename,'Resolution',300);
        filenumber = filenumber + 1;
    end
    iter = iter + 1;
end

