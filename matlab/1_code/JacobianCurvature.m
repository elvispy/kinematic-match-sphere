function curvatureMatrix = JacobianCurvature(eta, zk, dr, newCPoints)
    % THis function calculates the jacobian of the non lineal system
    % (which happens to be associated only to curvature)
    if size(eta, 2) > 1
        eta = eta';
    end
    Ntot = length(eta);
    sizze = Ntot - newCPoints;
    
    new_etas = [eta((newCPoints + 1):end); 0];
    if newCPoints > 0
        new_etas = [zk - sqrt(1 - dr^2 * (newCPoints - 1)^2); new_etas];
    end
    A1 = auxiliaryfunction1(new_etas, dr);
    A2 = auxiliaryfunction2(new_etas, dr, newCPoints);
    A3 = auxiliaryfunction3(new_etas, dr);
    A4 = auxiliaryfunction4(new_etas, dr, newCPoints);
   
    mainDiagonal = -2 * A1;
    upperDiagonal = A1 + A2 - A3 - A4; 
    lowerDiagonal = A1 - A2 + A3 + A4;
    if newCPoints == 0
        lowerDiagonal = [0; lowerDiagonal];
        mainDiagonal = [-4/(dr^2); mainDiagonal];
        upperDiagonal = [4/(dr^2); upperDiagonal];
    end
    curvatureMatrix = sparse([(Ntot + 1):(Ntot + sizze), (Ntot + 1):(Ntot + sizze-1), (Ntot + 2):(Ntot + sizze), Ntot + 1], ...
        [1:sizze, 2:sizze, 1:(sizze-1), 2*Ntot - newCPoints + 1], ...
        [mainDiagonal; upperDiagonal(1:(end-1)); ...
        lowerDiagonal(2:end); lowerDiagonal(1)], 2 * Ntot - newCPoints + 2, 2 * Ntot - newCPoints + 2);
    %zkterm = lowerDiagonal(1);
end

function result = auxiliaryfunction1(new_etas, dr)
    % New etas has all the information about the membrane to calculate the
    % curvature. (the first elements is possibly a function of z_k and the
    % last element is zero.
    n = length(new_etas);
    diffs2 = zeros(n-2, 1);
    for ii = 2:(n - 1)
        diffs2(ii - 1) = new_etas(ii + 1) - new_etas(ii - 1);
    end
     AUX = (1+ (diffs2/(2*dr)).^2);
     result = 1./(dr^2 * AUX.^(3/2));
     assert(size(result, 2) == 1);
end


function result = auxiliaryfunction2(new_etas, dr, newCPoints)
    % This will calculate the first auxiliary function
    % (the first element is possibly a function of z_k and the
    % last element is zero.
    n = length(new_etas);
    diffs2 = zeros(n-2, 1);
    for ii = 2:(n - 1)
        diffs2(ii - 1) = new_etas(ii+1) - new_etas(ii - 1);
    end
    AUX = (1+ (diffs2/(2*dr)).^2);
    if newCPoints == 0
        newCPoints = 1;
    end
    result = 1./(2 * dr^2 * (newCPoints:(n - 3  + newCPoints))' .* sqrt(AUX));
    assert(size(result, 2) == 1);
end

function result = auxiliaryfunction3(new_etas, dr)
    n = length(new_etas);
    diffs = zeros(n - 2, 1);
    for ii = 2:(n - 1)
        diffs(ii - 1) = new_etas(ii + 1) - new_etas(ii - 1);
    end
    NUM = 6*diffs .* diff(new_etas, 2);
    DEN = (1+ (diffs/(2*dr)).^2);
    result = NUM./(8* dr^4 * DEN.^(5/2));
    assert(size(result, 2) == 1);
end


function result = auxiliaryfunction4(new_etas, dr, newCPoints)
    n = length(new_etas);
    diffs = zeros(n - 2, 1);
    for ii = 2:(n - 1)
        diffs(ii-1) = new_etas(ii+1) - new_etas(ii - 1);
    end
    NUM = diffs.^2;
    DEN = (1+ (diffs/(2*dr)).^2);
    if newCPoints == 0
        newCPoints = 1;
    end
    result = NUM./(8* dr^4 * (newCPoints:(n - 3 + newCPoints))' .* DEN.^(3/2));
    assert(size(result, 2) == 1);
end