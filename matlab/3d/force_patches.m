clear; clc; close all;

% Import mesh generated in gmsh
surface_sphere
nodes = msh.POS'; 
elements = msh.TETS; 
% (assuming you haven't checked 'export all elements', else it will be a column of zeros.
elements = elements(:,1:4);
[x,y,z] = deal(nodes(1,:),nodes(2,:),nodes(3,:));

% From the gmsh script, annular_mesh.geo, we know that there are 
% 3 physical groups
ocTag = 1; % Tag-1: Outer circumference (oc)

ocLinesId = msh.TRIANGLES(:,end) == ocTag;

ocPointsId = unique(msh.TRIANGLES(ocLinesId,1:3));

nNodes = msh.nbNod;

% Create masks
ocMask = false(nNodes,1);
ocMask(ocPointsId) = true;

surface_idx = find(ocMask);
[xs,ys,zs] = deal(x(surface_idx),y(surface_idx),z(surface_idx));

theta = deg2rad(45);              % cone half-angle
   
% Direction vectors
right_dir = [1; 0; 0];
left_dir  = [-1; 0; 0];

right_patch = get_patches(nodes,surface_idx,theta,right_dir);
left_patch = get_patches(nodes,surface_idx,theta,left_dir);
    
figure;
% scatter3(x, y, z, 5, 'k'); 
hold on;
scatter3(x(right_patch), y(right_patch), z(right_patch), 30, 'r', 'filled');
scatter3(x(left_patch),  y(left_patch),  z(left_patch),  30, 'b', 'filled');
axis equal;
% legend('All nodes', 'Right patch', 'Left patch');
xlabel('x'); ylabel('y'); zlabel('z')

% Create patches for applying torque
% Get the top and bottom patches first
top_dir = [0; 0; 1];
bottom_dir  = [0; 0; -1];

top_patch = get_patches(nodes,surface_idx,theta,top_dir);
bottom_patch = get_patches(nodes,surface_idx,theta,bottom_dir);

scatter3(x(top_patch), y(top_patch), z(top_patch), 30, 'magenta', 'filled');
scatter3(x(bottom_patch),  y(bottom_patch),  z(bottom_patch),  30, 'black', 'filled');

% Get forces on the top
F0 = 1;
f_theta = azimuthal_force_direction(x(top_patch), y(top_patch), z(top_patch), F0);
quiver3(x(top_patch),y(top_patch),z(top_patch),f_theta(:,1)',f_theta(:,2)',f_theta(:,3)')

% Get forces on the bottom
F0 = -1;
f_theta = azimuthal_force_direction(x(bottom_patch), y(bottom_patch), z(bottom_patch), F0);
quiver3(x(bottom_patch),y(bottom_patch),z(bottom_patch),f_theta(:,1)',f_theta(:,2)',f_theta(:,3)')


% Forces for apply compression
fleft = ones(size(x,2),3);
quiver3(x(left_patch),y(left_patch),z(left_patch),fleft(left_patch,1)',0*fleft(left_patch,2)',0*fleft(left_patch,3)')

fright = ones(size(x,2),3);
quiver3(x(right_patch),y(right_patch),z(right_patch),-1*fright(right_patch,1)',0*fright(right_patch,2)',0*fright(right_patch,3)')

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


function f_theta = azimuthal_force_direction(x,y,z,F0)
% azimuthal_force_direction_array - Computes azimuthal force vectors for an array of positions.
    
    x = x(:); y = y(:); z = z(:);
    % Compute radius in xy-plane
    r = sqrt(x.^2 + y.^2)

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

    % Optional: warn about singular cases
    if any(~nonzero)
        warning('Some positions lie on the z-axis; azimuthal direction undefined at those points.');
    end
end
