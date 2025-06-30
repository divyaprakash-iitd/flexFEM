function [a,b,Aelem] = shapefunctioncoefficients(M,x,y)
    nnode   = size(M,2);    % No. of nodes in the finite element (triangle)
    ncor    = 2;            % No. of coordinates
    nelem   = size(M,1);    % No. of elements in the FE mesh
    
    % Shape function coefficients
    a = zeros(1,nnode,nelem);
    b = zeros(ncor,nnode,nelem);
    
    % Go over each of the elements
    % For each element a matrix is created, the cofactors of which give us 
    % the coefficients of the shape function
    A = ones(3,3);
    Aelem = zeros(1,nelem);
    for ielem = 1:nelem
        A(1,:) = x(M(ielem,:));
        A(2,:) = y(M(ielem,:));
        Aelem(ielem) = abs(det(A)); % Area of the element (Always positive)
        % Go over each node of an element
        % Calculate the coefficients for each node (a and b)
        C = cofactor(A);
        b(:,:,ielem) = C(1:end-1,:);
        a(:,:,ielem) = C(end,:);

        % Divide by area of the element
        b(:,:,ielem) = b(:,:,ielem)/Aelem(ielem);
        a(:,:,ielem) = a(:,:,ielem)/Aelem(ielem);
    end
end