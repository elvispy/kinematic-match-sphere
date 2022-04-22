function fc = curvature(eta, z_k, dr, newCPoints)
    % This function computes the contribution of the full curvature 
    % eta is a vector which has all the information about the position of the membrane
    % and dr is the spatial step
    if size(eta, 2) > 1
        eta = eta';
    end
    Ntot = length(eta); eta = eta((newCPoints+1):end); 
    if newCPoints > 0
        eta = [z_k - sqrt(1 - dr^2 * (newCPoints - 1)^2); eta];
    end
    fc = (2 * dr^2) * ones(newCPoints, 1); % dr^2 will dissappear  at line x+11
    eta = [eta; 0];
    n = length(eta);
    diffs2 = zeros(n-2, 1);
    for ii = 2:(n-1)
        diffs2(ii-1) = eta(ii+1) - eta(ii-1);
    end
    
    num1 = diff(eta, 2); 
    AUX = (1+ (diffs2/(2*dr)).^2);
    fc = [fc; num1./(AUX.^(3/2)) ...
        + diffs2 ./ (2* (max(newCPoints, 1):(Ntot-1))' .* sqrt(AUX))]/dr^2;
    if newCPoints == 0
        fc = [(4*eta(2) - 4 * eta(1))/dr^2; fc];
    end
    fc = [zeros(Ntot - newCPoints, 1); fc; zeros(2, 1)];
end