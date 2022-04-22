function jacobian_pieces = returnLinearizedMatrices(previous_jacobian_pieces, Ntot, dr, Mr)

% count how many entries are there
empty_indexes = sum(~cellfun(@isempty, previous_jacobian_pieces));
jacobian_pieces = previous_jacobian_pieces;

I_M = speye(Ntot);
Z_M = sparse(Ntot, Ntot);

A_prime = 2 * eye(Ntot);
A_prime(1, 1) = 4;
A_prime(1, 2) = -4;

for i = 2:Ntot
   A_prime(i, i-1) = -(2*i-3)/(2*i-2);
   if i < Ntot
       A_prime(i, i+1) = -(2*i-1)/(2*i-2);
   end
end
A_prime = (1 / dr^2) * A_prime;

A_2 = [-I_M; Z_M; sparse(2, Ntot)];
A_3 = [ Z_M; I_M; sparse(2, Ntot)];

for newCPoints = empty_indexes:(empty_indexes + 10)
    i = zeros(2*Ntot - newCPoints + 2, 1);
    j = zeros(2*Ntot - newCPoints + 2, 1);
    
    myA_prime = A_prime;
    myA_prime(1:newCPoints, :) = zeros(newCPoints, Ntot);
    myA_prime = sparse(myA_prime);
    % A_prime must be cutted, exact curvature. 
    myA_1 = [Z_M; myA_prime; sparse(2, Ntot)];
    
    residual_1 = myA_1(:, 1:newCPoints) + ...
        [eye(newCPoints); zeros(2*Ntot - newCPoints + 2, newCPoints)];
    % We add this term because in line 30
    % we didnt add a diagonal term
    myA_1 = myA_1(:, (newCPoints+1):end);
    
    i(1:(Ntot - newCPoints)) = 1:(Ntot - newCPoints);
    j(1:(Ntot - newCPoints)) = 1:(Ntot - newCPoints);
    
    myA_2 = A_2;
    residual_2 = myA_2(:, 1:newCPoints);
    myA_2 = myA_2(:, (newCPoints+1):end);
   
    i((Ntot - newCPoints + 1):(2*Ntot - 2*newCPoints)) = ...
        (Ntot + 1):(2*Ntot - newCPoints);
    j((Ntot - newCPoints + 1):(2*Ntot - 2*newCPoints)) = ...
        (Ntot -newCPoints + 1):(2*Ntot - 2*newCPoints);
    
    i((2*Ntot-2*newCPoints + 1):(2*Ntot - newCPoints)) = ...
        (Ntot - newCPoints + 1):Ntot;
    j((2*Ntot-2*newCPoints + 1):(2*Ntot - newCPoints)) = ...
        (2*Ntot - newCPoints + 2) * ones(newCPoints, 1);
    
    myA_3 = A_3(:, 1:newCPoints);
    myA_3(end, :) = -dr^2 * Mr * sparse(int_vector(newCPoints));
    
    myA_4 = [sum(residual_1, 2), sum(residual_2, 2)]; % Residuals dont have ones!
    myA_4(end-1, end) = -1;
    myA_4 = sparse(myA_4);
    myA = [myA_1, myA_2, myA_3, myA_4];
    myA = myA((newCPoints+1):end, :);
    
    i(end - 1) = 2 * Ntot - newCPoints + 1;
    j(end - 1) = 2 * Ntot - newCPoints + 1;
    i(end) = 2 * Ntot - newCPoints + 2;
    j(end) = 2 * Ntot - newCPoints + 2;
    
    %TODO: Update field names.
    jacobian_pieces{newCPoints + 1} = struct();
    jacobian_pieces{newCPoints+1}.A = myA;
    jacobian_pieces{newCPoints+1}.residual_1 = residual_1;
    jacobian_pieces{newCPoints+1}.ones = sparse(i, j, ones(length(j), 1));
    
end
end

function S = int_vector(n)
    % The integration vector without dhe dr^2 factor that goes into the
    % Pk block
    if n == 0
        S = zeros(1, 0);
    elseif n == 1
        S = [pi/12];
    else
        S = (pi) * ones(1, n);
        S(1) = 1/3 * S(1);
        for ii = 2:(n-1)
            S(ii) = 2*(ii-1) * S(ii);
        end
        S(n) = (3/2 * (n - 1) - 1/4) * S(n);
    end
end