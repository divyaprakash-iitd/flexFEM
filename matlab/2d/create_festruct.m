function festruct = create_festruct(M,x,y,FN,co,K)
    %% Construct the festruct
    festruct.M = M;
    festruct.x = [x, y];
    festruct.xorg = [x, y];
    festruct.FN = FN;
    festruct.co = co;
    festruct.K = K;
    festruct.nNodes = size(M,2);
    festruct.nElem = size(M,1);
    
    % Calculate shape coefficients
    [~,b,Aelem] = shapefunctioncoefficients(festruct.M,festruct.xorg(:,1),festruct.xorg(:,2));
    festruct.b = b;
    festruct.Aelem = Aelem;
end