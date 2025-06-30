function festruct = create_festruct(M,x,y,z,FN,co,K)
    %% Construct the festruct
    festruct.M = M;
    festruct.x = [x, y, z];
    festruct.xorg = [x, y, z];
    festruct.FN = FN;
    festruct.co = co;
    festruct.K = K;
    festruct.nNodes = size(M,2);
    festruct.nElem = size(M,1);
    
    % Calculate shape coefficients
    [~,b,Velem] = shapefunctioncoefficients(festruct.M,festruct.xorg(:,1),festruct.xorg(:,2),festruct.xorg(:,3));
    festruct.b = b;
    festruct.Velem = Velem;
end