function [a,b,Velem] = shapefunctioncoefficients(M,x,y,z)
    nnode   = size(M,2);    % No. of nodes in the finite element (triangle)
    ncor    = 3;            % No. of coordinates
    nelem   = size(M,1);    % No. of elements in the FE mesh
    
    % Shape function coefficients
    a = zeros(nnode,1,nelem);
    b = zeros(nnode,ncor,nelem);
    
    % Go over each of the elements
    % For each element a matrix is created, the cofactors of which give us 
    % the coefficients of the shape function
    A = ones(4,4);
    Velem = zeros(1,nelem);
    for ielem = 1:nelem
        A(:,1) = 1;
        A(:,2) = x(M(ielem,:));
        A(:,3) = y(M(ielem,:));
        A(:,4) = z(M(ielem,:));

        Velem(ielem) = (1/6)*abs(det(A)); % Volume of the element (Always positive)
        % Go over each node of an element
        % Calculate the coefficients for each node (a and b)
        C = cofactor(A);
        b(:,:,ielem) = (1/6)*C(:,2:end);
        a(:,:,ielem) = (1/6)*C(:,1);

        % Divide by the volume of the element
        b(:,:,ielem) = b(:,:,ielem)/Velem(ielem);
        a(:,:,ielem) = a(:,:,ielem)/Velem(ielem);
    end
end