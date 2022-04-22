% BDF2FullCurvature Solve next step of motion with dt time step
%   This function will solve the motion of the sphere with dt time step and
%   BDF2 Method on the temporal variable. Since PDE is nonlinear, Newton
%   Method will be used to find the solution. This implementation uses full
%   curvature and simulates vacuum conditions.

function [probableNextConditions, errortan] = BDF2FullCurvature(newCPoints, max_nb_cPoints, ...
                currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, ~)
            
    NEWTON_METHOD_TRESHOLD = 1e-8; % Treshold to determine if solution converged
    
    % Lets check the number of contact points
    if newCPoints < 0 || newCPoints > max_nb_cPoints
        errortan = Inf;
        probableNextConditions = struct('Eta_k', NaN, ...
            'u_k', NaN, 'z_k', NaN, 'v_k', NaN, 'P_k', NaN, ...
            'dt', NaN);
    else 
        %% Building R and R_prime
        gravityInfluence = [zeros(Ntot - newCPoints, 1); (-Fr * dt) * ones(Ntot, 1); 0; -Fr * dt]; %as placed in RHS of system
        f = @(x) sqrt(1-dr^2 * x.^2);
        
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
        
        A = dt * jacobian_pieces{newCPoints + 1}.A + ak * jacobian_pieces{newCPoints + 1}.ones;
        independent_term =  - bk * x_current - ck * x_previous;
        W = [Eta_k_current((newCPoints+1):end); u_k_current((newCPoints+1):end); 2*ones(newCPoints, 1); z_k_current;v_k_current]; % we cut the first equations
        % To artificially enter to the next while loop
        probable_conditions = W + 1;
        %Newton's method while loop
        flag = 1;
        while norm(probable_conditions - W) >= NEWTON_METHOD_TRESHOLD && flag < 1000
            flag = flag + 1;
            W = probable_conditions;
            z_k_probable = W(end - 1);
            Eta_k_probable = [z_k_probable - (f(0:(newCPoints-1))'); W(1:(Ntot-newCPoints))];
            probable_conditions = W - (A - dt * JacobianCurvature(Eta_k_probable, z_k_probable, dr, newCPoints))\ ...
                (A*W - dt * curvature(Eta_k_probable, z_k_probable, dr, newCPoints) - gravityInfluence - independent_term);
        end
        
        %% Post processing
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
    end % end main if
            
end