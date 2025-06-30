function festruct = calculate_forces(festruct)
    
    M = festruct.M;
    FN = festruct.FN;
    co = festruct.co;
    K = festruct.K;
    b = festruct.b;
    [x,y,z] = deal(festruct.x(:,1), festruct.x(:,2), festruct.x(:,3));
    Velem = festruct.Velem;

    nelem = size(M,1);  

    % Calculate the deformation gradient tensor, D
    D = zeros(3,3,nelem);
    for ielem = 1:nelem
        xmat = [x(M(ielem,:))';y(M(ielem,:))';z(M(ielem,:))'];
        belem = b(:,:,ielem);
        D(:,:,ielem) = xmat * belem;
    end
    
    % Calculate the Jacobian of transformation
    J = zeros(1,nelem);
    for ielem = 1:nelem
        J(ielem) = det(D(:,:,ielem));
    end
    
    % Calculate the derivative of the first invariant
    dIdx = zeros(3,4,nelem);
    for ielem = 1:nelem
        dIdx(:,:,ielem) = 2 * D(:,:,ielem) * b(:,:,ielem)';
    end
    
    % Calculate the derivative of the Jacobian
    dJdx = zeros(3,4,nelem);
    for ielem = 1:nelem
        DC = cofactor(D(:,:,ielem)); % Cofactor matrix of the Deformation Gradient
        dJdx(:,:,ielem) = DC * b(:,:,ielem)';
    end
    
    % Calculate the forces at all the nodes
    F = zeros(3,4,nelem);
    for ielem = 1:nelem
        F(:,:,ielem) = (-1/2)*co*Velem(ielem) * ...
                            (dIdx(:,:,ielem) - (2/J(ielem))*(1-(K/co)*log(J(ielem)))*dJdx(:,:,ielem));
    end
    
    % Add up the forces in each node
    for ielem = 1:nelem
        for inode = 1:4
            nodeId = M(ielem,inode);
            FN(:,nodeId) = FN(:,nodeId) + F(:,inode,ielem);
        end
    end
    festruct.FN = FN;
end