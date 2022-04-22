using .problemStructsAndFunctions: JacobianPiece, ProblemConditions
using SparseArrays

function returnMethod(method::String)::Function
    if occursin(r"linear"i, method)
        if occursin(r"BDF"i, method)
            return BDF2Linearized
        else
            return eulerLinearized;
        end
    else
        return eulerNewton;
    end
end

function eulerLinearized(new_nb_cpoints::Int64, max_nb_cpoints::Int64, 
        current_conditions::ProblemConditions, previous_conditions::ProblemConditions, 
        maximum_time_step::Float64, dt::Float64, 
        dr::Float64, Fr::Float64, Ntot::Int64, jacobian_pieces::Array{JacobianPiece}
    )::Tuple{ProblemConditions, Float64}

    if new_nb_cpoints < 0 || new_nb_cpoints > max_nb_cpoints
        errortan = Inf;
        probable_next_conditions = ProblemConditions(
            [NaN], [NaN], NaN, NaN, [NaN], NaN
        )
    else
        # Building R and R_prime
        gravity_influence = [zeros((Ntot - new_nb_cpoints, )); (-Fr * dt) * ones((Ntot, )); 0; (-Fr * dt)];
        exact_curvature = [zeros((Ntot - new_nb_cpoints, )); 2*dt * ones((new_nb_cpoints, )); 
                zeros((Ntot - new_nb_cpoints + 2, ))];
        f(x) = @. sqrt(1-dr^2 * x^2);

        residual_1 = zeros((2 * Ntot + 2, new_nb_cpoints));
        if new_nb_cpoints > 0
            val = dt * jacobian_pieces[new_nb_cpoints+1].residual[Ntot + new_nb_cpoints + 1, end];
            residual_1[Ntot + new_nb_cpoints + 1, end] = val;                
        end
        
        residual = residual_1 * f(0:(new_nb_cpoints - 1));
        residual = residual[(new_nb_cpoints + 1):end];

        ## Solvint the system
        η_k_current = current_conditions.η_k; 
        u_k_current = current_conditions.u_k;   
        z_k_current = current_conditions.z_k; 
        v_k_current = current_conditions.v_k;
        x = [η_k_current;u_k_current; z_k_current;v_k_current]; x = x[(new_nb_cpoints+1):end];
        if maximum_time_step == dt
            matrix = jacobian_pieces[new_nb_cpoints + 1].full_matrix;
        else
            matrix = (dt * jacobian_pieces[new_nb_cpoints + 1].time_dependent + jacobian_pieces[new_nb_cpoints + 1].constant)
        end
        results =  matrix\(x + gravity_influence + exact_curvature + residual);

        z_k_probable = results[end-1];
        v_k_probable = results[end];
        η_k_probable = [z_k_probable .- (f(0:(new_nb_cpoints-1))); results[1:(Ntot-new_nb_cpoints)]];

        u_k_probable = v_k_probable * ones((new_nb_cpoints, ));
        u_k_probable = [u_k_probable; results[(Ntot-new_nb_cpoints + 1):(2*Ntot - 2 * new_nb_cpoints)]];
        P_k_probable = results[(2*Ntot - 2 * new_nb_cpoints + 1):(2*Ntot - new_nb_cpoints)];

        # Calculating errortan
        errortan = 0;

        if new_nb_cpoints != 0
            g(x) = x/sqrt(1-x^2); ppt = (new_nb_cpoints - 1) * dr + dr/2;
            approximated_slope = (η_k_probable[new_nb_cpoints + 1] - η_k_probable[new_nb_cpoints])/dr;
            exact_slope = g(ppt);
            errortan = exact_slope - approximated_slope;
        end

        # Check whether the sphere and the membrane intersect
        for ii = (new_nb_cpoints + 1):max_nb_cpoints
            if η_k_probable[ii] > z_k_probable - f(ii-1);
                errortan = Inf;
                break;
            end
        end

        # Return structure
        probable_next_conditions = ProblemConditions(
            η_k_probable, 
            u_k_probable, 
            z_k_probable, 
            v_k_probable, 
            P_k_probable, 
            dt
        );
    end
    return probable_next_conditions, errortan
end


function BDF2Linearized(new_nb_cpoints::Int64, max_nb_cpoints::Int64, 
    current_conditions::ProblemConditions, previous_conditions::ProblemConditions, 
    maximum_time_step::Float64, dt::Float64, 
    dr::Float64, Fr::Float64, Ntot::Int64, jacobian_pieces::Array{JacobianPiece}
)::Tuple{ProblemConditions, Float64}

if new_nb_cpoints < 0 || new_nb_cpoints > max_nb_cpoints
    errortan = Inf;
    probable_next_conditions = ProblemConditions(
        [NaN], [NaN], NaN, NaN, [NaN], NaN
    )
else
    # Building R and R_prime
    gravity_influence = [zeros((Ntot - new_nb_cpoints, )); (-Fr * dt) * ones((Ntot, )); 0; (-Fr * dt)];
    exact_curvature = [zeros((Ntot - new_nb_cpoints, )); 2*dt * ones((new_nb_cpoints, )); 
            zeros((Ntot - new_nb_cpoints + 2, ))];
    f(x) = @. sqrt(1-dr^2 * x^2);

    residual_1 = zeros((2 * Ntot + 2, new_nb_cpoints));
    if new_nb_cpoints > 0
        val = dt * jacobian_pieces[new_nb_cpoints+1].residual[Ntot + new_nb_cpoints + 1, end];
        residual_1[Ntot + new_nb_cpoints + 1, end] = val;                
    end
    
    residual = residual_1 * f(0:(new_nb_cpoints - 1));
    residual = residual[(new_nb_cpoints + 1):end];

    ## Solving the system
    η_k_current = current_conditions.η_k; η_k_previous = previous_conditions.η_k;
    u_k_current = current_conditions.u_k; u_k_previous = previous_conditions.u_k;
    z_k_current = current_conditions.z_k; z_k_previous = previous_conditions.z_k;
    v_k_current = current_conditions.v_k; v_k_previous = previous_conditions.v_k;

    x_current  = [η_k_current;  u_k_current;  z_k_current;  v_k_current];  x_current  =  x_current[(new_nb_cpoints + 1):end];
    x_previous = [η_k_previous; u_k_previous; z_k_previous; v_k_previous]; x_previous = x_previous[(new_nb_cpoints + 1):end];

    rk = dt/current_conditions.dt;
    ak = (1+2*rk)/(1+rk);
    bk = -(1+rk);
    ck = rk^2/(1+rk);

    matrix = (dt * jacobian_pieces[new_nb_cpoints + 1].time_dependent + ak * jacobian_pieces[new_nb_cpoints + 1].constant)
    
    results =  matrix\(-bk * x_current + -ck * x_previous +  gravity_influence + exact_curvature + residual);

    z_k_probable = results[end-1];
    v_k_probable = results[end];
    η_k_probable = [z_k_probable .- (f(0:(new_nb_cpoints-1))); results[1:(Ntot-new_nb_cpoints)]];

    u_k_probable = v_k_probable * ones((new_nb_cpoints, ));
    u_k_probable = [u_k_probable; results[(Ntot-new_nb_cpoints + 1):(2*Ntot - 2 * new_nb_cpoints)]];
    P_k_probable = results[(2*Ntot - 2 * new_nb_cpoints + 1):(2*Ntot - new_nb_cpoints)];

    # Calculating errortan
    errortan = 0;

    if new_nb_cpoints != 0
        g(x) = x/sqrt(1-x^2); ppt = (new_nb_cpoints - 1) * dr + dr/2;
        approximated_slope = (η_k_probable[new_nb_cpoints + 1] - η_k_probable[new_nb_cpoints])/dr;
        exact_slope = g(ppt);
        errortan = exact_slope - approximated_slope;
    end

    # Check whether the sphere and the membrane intersect
    for ii = (new_nb_cpoints + 1):max_nb_cpoints
        if η_k_probable[ii] > z_k_probable - f(ii-1);
            errortan = Inf;
            break;
        end
    end

    # Return structure
    probable_next_conditions = ProblemConditions(
        η_k_probable, 
        u_k_probable, 
        z_k_probable, 
        v_k_probable, 
        P_k_probable, 
        dt
    );
end
return probable_next_conditions, errortan
end


function eulerNewton(new_nb_cpoints::Int64, max_nb_cpoints::Int64, 
    current_conditions::ProblemConditions, maximum_time_step::Float64, dt::Float64, 
    dr::Float64, Fr::Float64, Ntot::Int64, jacobian_pieces::Array{JacobianPiece}
    )::Tuple{ProblemConditions, Float64}
    return Array{ProblemConditions}(undef, 1);
end
