function [probableNextConditions, errortan] = eulerLinearized(newCPoints, max_nb_cPoints, ...
                currentConditions, ~, dt, dr, Fr, Ntot, jacobian_pieces, ~)
    % Lets check the number of contact points
    if newCPoints < 0 || newCPoints > max_nb_cPoints
        errortan = Inf;
        probableNextConditions = struct('Eta_k', NaN, ...
            'u_k', NaN, 'z_k', NaN, 'v_k', NaN, 'P_k', NaN, ...
            'dt', NaN);
    else
        %% Building R and R_prime
        R = [zeros(Ntot - newCPoints, 1); (-Fr * dt) * ones(Ntot, 1); 0 ; -Fr * dt];
        for ii = (Ntot - newCPoints + 1):Ntot
            R(ii) = R(ii) + 2 * dt;
        end
        f = @(x) sqrt(1-dr^2 * x.^2);
        residual_1 = jacobian_pieces{newCPoints + 1}.residual_1;
        if newCPoints > 0
            residual_1(Ntot + newCPoints + 1, end) = ...
                dt * residual_1(Ntot + newCPoints + 1, end); 
        end
        R_prime = residual_1 * (f(0:(newCPoints-1))');
        R_prime = R_prime((newCPoints+1):end);
        
        %% Solving the system
        Eta_k_current = currentConditions.Eta_k; 
        u_k_current   = currentConditions.u_k;   
        z_k_current   = currentConditions.z_k; 
        v_k_current   = currentConditions.v_k;
        x = [Eta_k_current;u_k_current; z_k_current;v_k_current]; x = x((newCPoints+1):end);
        results = (dt * jacobian_pieces{newCPoints + 1}.A + jacobian_pieces{newCPoints + 1}.ones) ...
            \(x + R + R_prime);

        z_k_probable = results(end-1);
        v_k_probable = results(end);
        Eta_k_probable = [z_k_probable - (f(0:(newCPoints-1))'); results(1:(Ntot-newCPoints))];

        u_k_probable = v_k_probable * ones(newCPoints, 1);
        u_k_probable = [u_k_probable; results((Ntot-newCPoints + 1):(2*Ntot - 2 * newCPoints))];
        P_k_probable = results((2*Ntot - 2 * newCPoints + 1):(2*Ntot - newCPoints));
        
        %% Calculating errortan
        errortan = 0;
        
        if newCPoints ~= 0
            g = @(x) x/sqrt(1-x^2); ppt = (newCPoints-1)*dr + dr/2;
            approximateSlope = (Eta_k_probable(newCPoints+1)-Eta_k_probable(newCPoints))/dr;
            exactSlope = g(ppt);
            errortan = (exactSlope) - (approximateSlope);
        end
        
        % Check whether the sphere and the membrane intersect
        for ii = (newCPoints + 1):max_nb_cPoints
            if Eta_k_probable(ii) > z_k_probable - f(ii-1)
                errortan = Inf;
                break;
            end
        end % end for
        
        % Return structure
        probableNextConditions = struct('Eta_k', Eta_k_probable, ...
            'u_k', u_k_probable, 'z_k', z_k_probable, 'v_k', v_k_probable, ...
            'P_k', P_k_probable, 'dt', dt);
    end % end main if
end % end function