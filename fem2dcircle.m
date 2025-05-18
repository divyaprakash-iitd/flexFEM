clear; clc; close all;

%% Description: Solves a beam problem using FEM

%% Define the elastic constants
K   = 10000;
co  = 5000;
E   = co*(2*co+3*K)/(co+K);
nu  = K/(2*(co+K));

%% Geometry and FEM Mesh
d   = 1;        % Unit depth for 2D cases
% Import gmsh mesh
annular_mesh;
% Get the points
points = msh.POS;
[M,x,y] = deal(msh.TRIANGLES(:,1:3), points(:,1), points(:,2));
% From the gmsh script we know that there are 3 physical groups
ocTag = 1; % Tag-1: Outer circumference (oc)
icTag = 2; % Tag-2: Inner circumference (ic)
asTag = 3; % Tag-3: The annular surface (as)

ocLinesId = msh.LINES(:,end) == ocTag;
icLinesId = msh.LINES(:,end) == icTag;

ocPointsId = unique(msh.LINES(ocLinesId,1:2));
icPointsId = unique(msh.LINES(icLinesId,1:2));

nNodes = msh.nbNod;

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

% Storing the original configuration
xorg = x;
yorg = y;

%% Calculate the shape function coefficients
nnode   = size(M,2);     % No. of nodes in the finite element (triangle)
ncor    = 2;            % No. of coordinates
nelem   = size(M,1);    % No. of elements in the FE mesh
[a,b,Aelem] = shapefunctioncoefficients(M,x,y);

%% Forces
nglobal = numel(x); % Total number of nodes
FN      = zeros(nglobal,2);

%% Run time iterations
niter   = 50000;
dt      = 0.001;
mu      = 1000;
filenumber = 1;

%% Applied Force
Ftip        = 500;          % Tip force
direction   = 1;            % 1: Lateral (x-direction); 2: Axial (y-direction)  

%% Plot controls
plotinterval = niter/100;

fig = figure();
fig.Position = [675, 549, 570, 413];
triplot(M,x,y)
xlabel('x [m]')
ylabel('y [m]')
hold on

%% Store simulation images in a directory
folder = 'images';

% Check if the folder exists, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end


%% Time loop
tol     = 1e-7;
iter    = 1;
while iter <= niter
    FApplied = zeros(nglobal,2);
    %% Apply tip force
    if iter > niter/2
        Ftip = 0;
    end

    FApplied(leftMask,direction)   = -Ftip;
    FApplied(rightMask,direction)   = Ftip;

    FN = FApplied;

    %% Calculate the deformation tensor
    % New locations of the nodes
    % Transpose is done below since the reshape function executes columnwise
    xt = reshape(x(M)',1,nnode,nelem);
    yt = reshape(y(M)',1,nnode,nelem);
    xvec = cat(1,xt,yt);
    
    % There will be one Deformation tensor per element
    D = pagemtimes(xvec,pagetranspose(b));
    
    %% Calculate the Invariant Derivative term
    dIdx = 2*pagemtimes(D,b);
        
    %% Calculate the Jacobian Derivative term
    mf = ones(2,2,nelem); % Multiplication factor
    mf(1,2,:) = -1; mf(2,1,:) = -1;
    dJdx = pagemtimes(mf.*fliplr(flip(D)), b);
    
    %% Calculate the Jacobian
    J = zeros(1,nelem);
    for ielem = 1:nelem
        J(ielem) = det(D(:,:,ielem));
    end
    
    %% Calculate the forces at each node
    F = zeros(size(b));
    for ielem = 1:nelem
        F(:,:,ielem) = -1/2 * co * Aelem(ielem) * d *... 
                (   dIdx(:,:,ielem) - 2/J(ielem) * ... 
                    (1 - K/co*log(J(ielem))) * dJdx(:,:,ielem)  );
    end
    
    %% Sum up the forces at each node   
    FF = pagetranspose(F);
    Fval=reshape(permute(FF,[2 1 3]),size(FF,2),[])';
    Fid = reshape(M',nelem*nnode,[]);
    
    for idindex = 1:numel(Fid)
        for icor = 1:ncor
            FN(Fid(idindex),icor) = FN(Fid(idindex),icor) + Fval(idindex,icor);
        end
    end

    %% Calculate the deformation
    for inode = 1:nglobal
        x(inode) = x(inode) + dt/mu * FN(inode,1);
        y(inode) = y(inode) + dt/mu * FN(inode,2);
    end

    %% Plot data points
    if mod(iter,plotinterval) == 0 && iter ~= niter
        cla;
        triplot(M,x,y,'r','DisplayName','Deformed')
        h1 = quiver(x(leftMask),y(leftMask),FApplied(leftMask,1),FApplied(leftMask,2),'Color','k');
        quiver(x(rightMask),y(rightMask),FApplied(rightMask,1),FApplied(rightMask,2), 'Color','k')
        h2 = quiver(x(leftMask),y(leftMask),FN(leftMask,1),FN(leftMask,2), 'Color','b');
        quiver(x(rightMask),y(rightMask),FN(rightMask,1),FN(rightMask,2), 'Color','b')

        title(sprintf('Iteration = %d',iter))
        lgd = legend([h1, h2], 'F_{Applied}', 'F_{Total}');
        lgd.Location = 'southeast';
        daspect([1,1,1])
        xlim([-1.7,1.7])
        ylim([-1,1])
        axis equal
        filename = sprintf('%s/image_%02d.png', folder, filenumber);
        exportgraphics(gcf,filename,'Resolution',300);
        filenumber = filenumber + 1;
    end
    iter = iter + 1;
end