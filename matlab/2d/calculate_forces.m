function festruct = calculate_forces(festruct)
    
    d = 1; % Since it's 2D, unit depth is unity aka 1
    M = festruct.M;
    x = festruct.x(:,1);
    y = festruct.x(:,2);
    xorg = festruct.xorg(:,1);
    yorg = festruct.xorg(:,2);
    FN = festruct.FN;
    co = festruct.co;
    K = festruct.K;
    b = festruct.b;
    Aelem = festruct.Aelem;

    nNodes   = festruct.nNodes;     % No. of nodes in the finite element (triangle)
    nElem   = festruct.nElem;    % No. of elements in the FE mesh
    ncor    = 2;            % No. of coordinates


    %% Calculate the deformation tensor
    % New locations of the nodes
    % Transpose is done below since the reshape function executes columnwise
    xt = reshape(x(M)',1,nNodes,nElem);
    yt = reshape(y(M)',1,nNodes,nElem);
    xvec = cat(1,xt,yt);
    
    % There will be one Deformation tensor per element
    D = pagemtimes(xvec,pagetranspose(b));
    
    %% Calculate the Invariant Derivative term
    dIdx = 2*pagemtimes(D,b);
        
    %% Calculate the Jacobian Derivative term
    mf = ones(2,2,nElem); % Multiplication factor
    mf(1,2,:) = -1; mf(2,1,:) = -1;
    dJdx = pagemtimes(mf.*fliplr(flip(D)), b);
    
    %% Calculate the Jacobian
    J = zeros(1,nElem);
    for ielem = 1:nElem
        J(ielem) = det(D(:,:,ielem));
    end
    
    %% Calculate the forces at each node
    F = zeros(size(b));
    for ielem = 1:nElem
        F(:,:,ielem) = -1/2 * co * Aelem(ielem) * d *... 
                (   dIdx(:,:,ielem) - 2/J(ielem) * ... 
                    (1 - K/co*log(J(ielem))) * dJdx(:,:,ielem)  );
    end
    
    %% Sum up the forces at each node   
    FF = pagetranspose(F);
    Fval=reshape(permute(FF,[2 1 3]),size(FF,2),[])';
    Fid = reshape(M',nElem*nNodes,[]);
    
    for idindex = 1:numel(Fid)
        for icor = 1:ncor
            FN(Fid(idindex),icor) = FN(Fid(idindex),icor) + Fval(idindex,icor);
        end
    end
    festruct.FN = FN;
end