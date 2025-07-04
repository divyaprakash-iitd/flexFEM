clear; clc; close all;


%% Add path
libpath = '../../'; % Replace with the path where you have downloaded the source code
addpath(libpath);

%% Description: Solves a demo problem using FEM (3D)

%% Define the elastic constants
K   = 100000;
co  = 50000;
E   = co*(2*co+3*K)/(co+K);
nu  = K/(2*(co+K));

% Import mesh generated in gmsh
ellipsoidal_shell_mesh;
nodes = msh.POS'; 
M = msh.TETS; 
% (assuming you haven't checked 'export all M', else it will be a column of zeros.
M = M(:,1:4);
nelem = size(M,1);
ncor = size(nodes,1);
nnode = size(nodes,2);

%% Generate Patches Masks
[x,y,z] = deal(nodes(1,:)',nodes(2,:)',nodes(3,:)'); % Note the transpose

% From the gmsh script
ocTag = 1; % Tag-1: Outer circumference (oc)
ocLinesId = msh.TRIANGLES(:,end) == ocTag;
ocPointsId = unique(msh.TRIANGLES(ocLinesId,1:3));
nNodes = msh.nbNod;

% Create masks
ocMask = false(nNodes,1);
ocMask(ocPointsId) = true;
surface_idx = find(ocMask);
theta = deg2rad(45);              % cone half-angle
% Direction vectors
top_dir = [0; 0; 1];
bottom_dir  = [0; 0; -1];

top_patch = get_patches(nodes,surface_idx,theta,top_dir);
bottom_patch = get_patches(nodes,surface_idx,theta,bottom_dir);

% Generate shape function coefficients
[a,b,Velem] = shapefunctioncoefficients(M,x,y,z);

% Run time iterations
niter   = 15000;
dt      = 0.001;
mu      = 1000;

% Store simulation images in a directory
folders = {'image', 'vtk'};  % Use a cell array of strings

% Check if each folder exists, if not, create it
for i = 1:length(folders)
    if ~exist(folders{i}, 'dir')
        mkdir(folders{i});
    end
end

%% Time Loop
iter    = 1;
fmag    = 500;
fig = figure;
fig.Position = [1 1 960 961];
fId = 1;

% Initialize forces on the nodes
FN = zeros(3,nnode);
festruct = create_festruct(M,x,y,z,FN,co,K);

while iter <= niter
    % Initialize forces on the nodes
    FN = zeros(3,nnode);
    if iter < niter/2

        % Apply stretching forces
        FN(3,top_patch) = fmag;
        FN(3,bottom_patch) = -fmag;

        % Apply twisting forces
        % Get forces 
        f_theta = azimuthal_force_direction(x(bottom_patch), y(bottom_patch), z(bottom_patch), fmag);
        FN(:,bottom_patch) = f_theta';
        
        f_theta = azimuthal_force_direction(x(top_patch), y(top_patch), z(top_patch), fmag);
        FN(:,top_patch) = -f_theta';

    end
    
    festruct.FN = FN;
    festruct = calculate_forces(festruct);
    

    % Read updated values from festruct
    x = festruct.x(:,1);
    y = festruct.x(:,2);
    z = festruct.x(:,3);
    FN = festruct.FN;


    % Calculate the deformation
    for inode = 1:nnode
        festruct.x(inode,1) = festruct.x(inode,1) + dt/mu * FN(1,inode);
        festruct.x(inode,2) = festruct.x(inode,2) + dt/mu * FN(2,inode);
        festruct.x(inode,3) = festruct.x(inode,3) + dt/mu * FN(3,inode);
    end
    
    if (mod(iter,50) == 0)
        fname = sprintf("image/img_%05d.png",fId);
        
        XX = [x(:) y(:) z(:)];
        tetramesh(M,XX,'FaceColor','b')
        camlight
        lim = 1.5;
        xlim([-lim,lim])
        ylim([-lim,lim])
        zlim([-lim,lim])
        daspect([1,1,1])
        title(sprintf("Iteration = %d",iter))
        %view(0, 0)
        pause(0.1)

        saveas(gcf,fname)
       
        % Save vtk files for visualization using paraview 
        filename = sprintf('vtk/img_%05d.vtk',fId);
        data_struct.type = 'vector';
        data_struct.name = 'Force';
        data_struct.data = FN';
        data_title = 'Force';
        vtx_coord = [x(:), y(:), z(:)];
        flipped = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following files needs to be downloaded and added to MATLAB's path or kept in the same directory as this code.
% If VTK files are not needed then comment out the following line.
        % Shawn Walker (2025). write BINARY VTK file for tetrahedral grid with scalar and vector data (https://www.mathworks.com/matlabcentral/fileexchange/58002-write-binary-vtk-file-for-tetrahedral-grid-with-scalar-and-vector-data), MATLAB Central File Exchange. Retrieved June 30, 2025. 
        stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,vtx_coord,M,data_struct,flipped);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fId = fId + 1;

    end

    iter = iter + 1;
end

function f_theta = azimuthal_force_direction(x,y,z,F0)
% azimuthal_force_direction_array - Computes azimuthal force vectors for an array of positions.
    x = x(:); y = y(:); z = z(:);
    % Compute radius in xy-plane
    r = sqrt(x.^2 + y.^2);

    % Handle scalar or vector F0
    if isscalar(F0)
        F0 = F0 * ones(size(r));
    end

    % Initialize force array
    f_theta = zeros(size(x,1),3);

    % Avoid division by zero: r == 0 (point on z-axis)
    nonzero = r > 1e-12;

    % Compute azimuthal direction only for valid points
    f_theta(nonzero, 1) = -F0(nonzero) .* y(nonzero) ./ r(nonzero);
    f_theta(nonzero, 2) =  F0(nonzero) .* x(nonzero) ./ r(nonzero);
    % f_theta(:,3) is already zero

%     % Optional: warn about singular cases
%     if any(~nonzero)
%         warning('Some positions lie on the z-axis; azimuthal direction undefined at those points.');
%     end
end

function patch = get_patches(nodes,surface_idx,theta,direction)
    cos_theta = cos(theta);
    
    % Get surface node vectors
    vecs = nodes(:, surface_idx);     % 3 × N
    
    % Normalize them
    vecs_norm = vecs ./ vecnorm(vecs);

    % Find those within the cone
    dot_product = vecs_norm' * direction;
    
    patch = surface_idx(dot_product > cos_theta);

end
