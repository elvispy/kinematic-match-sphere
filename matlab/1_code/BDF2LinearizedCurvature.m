function [probableNextConditions, errortan] = BDF2LinearizedCurvature(newCPoints, max_nb_cPoints, ...
                currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, ~)
    % BDF2LinearizedCurvature Advances the inear system one time step with
    % linearized curvature and BDF2 method.
    
    
    if newCPoints < 0 || newCPoints > max_nb_cPoints
        errortan = Inf;
        probableNextConditions = struct('Eta_k', NaN, ...
            'u_k', NaN, 'z_k', NaN, 'v_k', NaN, 'P_k', NaN, ...
            'dt', NaN);
    else        
        
        %% Building the left side of the system
        % Gravity Influence 
        gravityInfluence = [zeros(Ntot - newCPoints, 1); (-Fr * dt) * ones(Ntot, 1); 0 ; -Fr * dt];
        % Exact curvature Term
        for ii = (Ntot - newCPoints + 1):Ntot
            gravityInfluence(ii) = gravityInfluence(ii) + 2 * dt;
        end
        f = @(x) sqrt(1-dr^2 * x.^2);
        residual_1 = jacobian_pieces{newCPoints + 1}.residual_1;
        if newCPoints > 0
            residual_1(Ntot + newCPoints + 1, end) = ...
                dt * residual_1(Ntot + newCPoints + 1, end); 
        end
        sphereHeightTerm = residual_1 * (f(0:(newCPoints-1))');
        sphereHeightTerm = sphereHeightTerm((newCPoints+1):end);
        assert(sum(abs(sphereHeightTerm) > 1e-10) == (newCPoints > 0));
        
        %% Solving the system
        Eta_k_current = currentConditions.Eta_k; Eta_k_previous = previousConditions.Eta_k;
        u_k_current   = currentConditions.u_k;   u_k_previous   = previousConditions.u_k;
        z_k_current   = currentConditions.z_k;   z_k_previous   = previousConditions.z_k;
        v_k_current   = currentConditions.v_k;   v_k_previous   = previousConditions.v_k;
        
        x_current  = [Eta_k_current;  u_k_current;  z_k_current;  v_k_current];  x_current  =  x_current((newCPoints + 1):end);
        x_previous = [Eta_k_previous; u_k_previous; z_k_previous; v_k_previous]; x_previous = x_previous((newCPoints + 1):end);
        
        rk = dt/currentConditions.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        probable_conditions = (dt * jacobian_pieces{newCPoints + 1}.A + ak * jacobian_pieces{newCPoints + 1}.ones) \ ...
            ( - bk * x_current - ck * x_previous + gravityInfluence + sphereHeightTerm);
        z_k_probable = probable_conditions(end-1);
        v_k_probable = probable_conditions(end);
        Eta_k_probable = [z_k_probable - (f(0:(newCPoints-1))'); probable_conditions(1:(Ntot-newCPoints))];

        u_k_probable = v_k_probable * ones(newCPoints, 1);
        u_k_probable = [u_k_probable; probable_conditions((Ntot-newCPoints + 1):(2*Ntot - 2 * newCPoints))];
        P_k_probable = probable_conditions((2*Ntot - 2 * newCPoints + 1):(2*Ntot - newCPoints));
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
        
        probableNextConditions = struct('Eta_k', Eta_k_probable, ...
            'u_k', u_k_probable, 'z_k', z_k_probable, 'v_k', v_k_probable, ...
            'P_k', P_k_probable, 'dt', dt);
    end %end main if
      
end
