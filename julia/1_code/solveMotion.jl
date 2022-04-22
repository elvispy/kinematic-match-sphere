"""
    Author: Elvis Aguero.
    Date: 13th April, 2022.
    solveMotion
"""
#Includes
using Printf;
using Plots;
default(legend = false);
using JLD2, FileIO;
using Dates;
using CSV;
using DataFrames;
include("./problemStructsAndFunctions.jl");
using .problemStructsAndFunctions # Specialized structs and functions to be used in this code

# Main function
function solve_motion(; # <== Keyword arguments!   
     rS::Float64 = 2.38,       Tm::Float64 = 107.0,            R_f::Float64 = NaN,         η_k::Vector{Float64} = Float64[],
     mS::Float64 = NaN,         g::Float64 = 9.80665e-3,       z_k::Float64 = Inf,          u_k::Vector{Float64} = Float64[],    
      μ::Float64 = 0.3,         N::Int64 = 100, save_after_contact::Bool = false,           P_k::Vector{Float64} = Float64[],     
plotter::Bool = false, recorded_time::Float64 = 0.0,     file_name::String="historial.csv", v_k::Float64 = -0.6312, 
 method::String = "EulerLinearized",                   export_data::Bool = true,            ρ_a::Float64 = 1.225e-3,
 simul_time::Float64 = NaN
 )

    if isnan(simul_time) == true
        simul_time = 20; 
    end
    if isnan(R_f); R_f = 52.5/rS; end
    if isnan(mS);  mS = 3.25 * 4 * pi * (rS.^3) / 3; end
    # Constants definitions
    Fr = (μ * rS * g)/Tm;
    Mr = μ * (rS^2) / mS;
    D = ρ_a * rS / μ;

    # Units
    length_unit   = rS; # Spatial unit of measurement (in mm)
    velocity_unit = sqrt(Tm/μ); # Velocity in mm/mS
    time_unit     = length_unit/velocity_unit; # Temporal unit of measurement (in ms)
    pressure_unit = μ * length_unit / time_unit^2; # Unit of pressure (in mg/(ms^2 * mm))

    # Numerical constants
    Ntot = ceil(Int64, R_f * N); # Total number of non-trivial points in the spatial coordinate
    dr = 1/N; # Spatial step

    ## Initial conditions
    P_k = P_k/pressure_unit; # Current vector of pressures. Only non trivial points saved (Empty vector == membrane and ball should not be touching)
    contact_points = length(P_k); # Initial number of contact points
    max_nb_contact_points = N + 1; # Maximum allowed number of contact points

    # Set η_k initial conditions
    if length(η_k) != Ntot
      η_k = initialConditionsLinearized(Fr, dr, Ntot); # TODO: Define this method and treat other cases (BDF2, Full curvature, etc)
    else
      if size(η_k, 2) > 1; η_k = vec(transpose(η_k)); end 
      η_k = η_k/length_unit;
    end

    # Set u_k initial conditions
    if Ntot == length(u_k)
        if size(u_k, 2) > 1; u_k = transpose(u_k); end
        u_k = u_k/velocity_unit;
    else
        u_k = zeros((Ntot, )); # Initial vertical velocity of the membrane (a.k.a. d/dt Eta_k = u_k)    
    end

    # Set sphere's height initial conditions
    if z_k == Inf
        z_k = η_k[1] + 1; # Curr1ent position of the center of the ball (Dimensionless)
    else
        @assert z_k >= η_k[1] + 1 "Sphere is not in a valid position!"
        z_k = z_k/length_unit;
    end


    # Zero height to measure contact time experimentally
    ZERO_HEIGHT = η_k[1];

    # To store our zero potential: place of the center of membrane at impact
    zero_potential = NaN;

    probable_next_conditions = Array{ProblemConditions}(undef, 5);
    if recorded_time > 0; dt = recorded_time/time_unit; end

    # Set sphere's velocity initial condition
    v_k = v_k/velocity_unit;

    # TODO: Implement other cases (i.e. BDF2, Newton, etc)    
    if recorded_time == 0; dt = 1/( N*time_unit); end
    current_conditions = ProblemConditions(η_k, u_k, z_k, v_k, P_k, dt);
    previous_conditions = retrievePreviousConditions(current_conditions, method); # TODO: implement one step method for BDF2.

    final_time = simul_time/time_unit; # Maximum time (dimensionless) to be simulated.
    initial_time = 0/time_unit; # Initial time (dimensionless)
    maximum_time_step = dt; # Interval of time to be recorded in output matrix.
    t = initial_time; # t is the variable that will keep track of the time
    current_index = 2; # Index of current matrix to save 

    initial_velocity = v_k * velocity_unit;
    impact_velocity = NaN;

    # Choose appropriate solver (Euler, BDF2, Newton, etc)
    getNextStep = returnMethod(method);

    returnMatrices = matrixGenerationMethod(method)
    jacobian_pieces = Array{JacobianPiece}(undef, max_nb_contact_points + 1);
    jacobian_pieces = returnMatrices(jacobian_pieces, Ntot, dr, Mr, maximum_time_step);
    
    # For preallocation efficiency
    maximum_index = ceil(Int64, (final_time - initial_time)/maximum_time_step) + 4;
    number_of_extra_indexes = 0;

    ## Flow control variables
    contact_time = 0; lab_contact_time = 0;
    maximum_deflection = 0;
    first_fall_flag = true; # Boolean to record maximum deflection time
    contact_flag = false; # To record first contact time
    velocity_out_flag = false; lab_velocity_out_flag = false; # To check if mechanical energy out was recorded
    iii = 0; jjj = 0; # Indexes to keep track  of subdivisions of maximum_time_step
    resetter = 0; # This variable controls how fast dt can grow
    
    # Check file existence, create .csv
    A = match(r"[\\\/]([^\\\/]+)\.jl$", PROGRAM_FILE);
    if A !== nothing
        A = A[1];
        SAVE_PATH = "./2_pipeline/$A/out";
    else
        SAVE_PATH = "./2_pipeline/default/out";
    end
    if isdir(SAVE_PATH) == false; mkpath(SAVE_PATH); end
    if isfile("$SAVE_PATH/$file_name") == false
        headers_data_frame = DataFrame(
            ID  = [],
            vi = [],
            surfacetension = [],
            radius = [],
            cTime = [],
            maxDeflection = [],
            coefOfRestitution = [],
            Density = [],
            labcTime = [],
            labcoefOfRestitution = []
        )
        CSV.write("$SAVE_PATH/$file_name", headers_data_frame)
    end
   
    if plotter == true
        plot_width = ceil(Int64, 3 * N);
        xplot = LinRange(0, plot_width/N, plot_width);
        ηX = [-xplot[end:-1:2]; xplot] * length_unit;
        # ηU  = zeros(1, 2 * plot_width - 1);
        #step = ceil(Int64, N/15);

        θ = LinRange(0, 2*π, 500)
        circleX = length_unit * sin.(θ);
        circleY = length_unit * cos.(θ);
        # For animation:
        # https://stackoverflow.com/questions/58103056/erasing-previous-data-plots-in-julia-plots-jl-gr-backend

    end

    ## Preparing post-processing

    # Preallocate variables that will be exported
    recorded_z = zeros((maximum_index, )); recorded_z[1] = z_k * length_unit;
    recorded_η = zeros(maximum_index, Ntot); recorded_η[1, :] = η_k * length_unit;
    recorded_u = zeros(maximum_index, Ntot); recorded_u[1, :] = u_k * velocity_unit;
    recorded_P = zeros(maximum_index, max_nb_contact_points);
    recorded_v = zeros((maximum_index, )); recorded_v[1] = v_k * velocity_unit;
    recorded_t = zeros((maximum_index, )); recorded_t[1] = initial_time * time_unit;

    # Coefficient of restitution
    mechanical_energy_in = NaN;
    mechanical_energy_out = NaN; lab_mechanical_energy_out = NaN;

    indexes_to_save = zeros(Int64, (maximum_index, )); indexes_to_save[1] = 1;
    current_to_save = 2;

    ## Main loop
    while (t <= final_time)
        errortan = Inf * ones((5, ));
        recalculate = false;

        # If empty, allocate more space for matrices
        if isassigned(jacobian_pieces, contact_points + 3) == false
            jacobian_pieces = returnMatrices(jacobian_pieces, Ntot, dr, Mr, maximum_time_step);
        end

        # First, we try to solve with the same number of contact points
        probable_next_conditions[3], errortan[3] = getNextStep(contact_points, max_nb_contact_points,
            current_conditions, previous_conditions, maximum_time_step, dt, dr, Fr, Ntot, jacobian_pieces);

        if abs(errortan[3]) < 1e-8 # If almost no error, we accept the solution
            previous_conditions = current_conditions;
            current_conditions  = probable_next_conditions[3];

        else # If there is some error, we try with diffferent contact points
            # Lets try with one more point
            probable_next_conditions[4], errortan[4] = getNextStep(contact_points + 1, max_nb_contact_points,
                current_conditions, previous_conditions, maximum_time_step, dt, dr, Fr, Ntot, jacobian_pieces);
            # Lets try with one point less
            probable_next_conditions[2], errortan[2] = getNextStep(contact_points - 1, max_nb_contact_points,
                current_conditions, previous_conditions, maximum_time_step, dt, dr, Fr, Ntot, jacobian_pieces);
            
            if (abs(errortan[3]) > abs(errortan[4]) ||  (abs(errortan[3]) > abs(errortan[2])))
                if abs(errortan[4]) <= abs(errortan[2])
                    # Now lets check with one more point to be sure
                    _, errortan[5] = getNextStep(contact_points + 2, max_nb_contact_points, 
                        current_conditions, previous_conditions, maximum_time_step, dt, dr, Fr, Ntot, jacobian_pieces);

                    if abs(errortan[4]) < abs(errortan[5])
                        # Accept new data 
                        previous_conditions = current_conditions;
                        current_conditions  = probable_next_conditions[4];
                        contact_points      = contact_points + 1;
                    else
                        recalculate = true;
                    end
                else
                    # now lets check if errortan is good enough with one point less
                    _, errortan[1] = getNextStep(contact_points - 2, max_nb_contact_points,
                    current_conditions, previous_conditions, maximum_time_step, dt, dr, Fr, Ntot, jacobian_pieces);

                    if abs(errortan[2]) < abs(errortan[1])
                        # Accept new data
                        previous_conditions = current_conditions;
                        current_conditions  = probable_next_conditions[2];
                        contact_points      = contact_points - 1;
                    else
                        recalculate = true;
                    end 

                end # End of errortan[4] <= errortan[2]
            else # The same number of contact points may be the best
                if errortan[3] == Inf # ==> All errors are infinity
                    recalculate = true;
                else
                    # Accept new data
                    previous_conditions = current_conditions;
                    current_conditions  = probable_next_conditions[3];
                end
            end
        end

        if recalculate == true
            dt = dt/2;
            # Refine time step in index notation 
            iii = iii + 1; jjj = 2 * jjj;
        else
            t = t + dt; jjj = jjj + 1;
            if mod(jjj, 2) == 0 && resetter > 0
                jjj = div(jjj, 2); 
                iii = iii - 1;
                # Increase time step
                dt = 2 * dt;
                # Decrease the number of time you can make dt bigger
                resetter = resetter - 1;
            end

            # Update Indexes if necessary and store results into a matrix
            if current_index > maximum_index
                recorded_z = [recorded_z; zeros(number_of_extra_indexes, 1)];
                recorded_η = [recorded_η; zeros(number_of_extra_indexes, Ntot)];
                recorded_u = [recorded_u; zeros(number_of_extra_indexes, Ntot)];
                recorded_P = [recorded_P; zeros(number_of_extra_indexes, max_nb_contact_points)];
                recorded_v = [recorded_v; zeros(number_of_extra_indexes, 1)];
                recorded_t = [recorded_t; zeros(number_of_extra_indexes, 1)];
                maximum_index = maximum_index + number_of_extra_indexes;
                number_of_extra_indexes = 0;
            end

            # Stored data
            recorded_z[current_index] = current_conditions.z_k * length_unit;
            recorded_η[current_index, :] = current_conditions.η_k * length_unit;
            recorded_u[current_index, :] = current_conditions.u_k * velocity_unit;
            recorded_P[current_index, 1:length(current_conditions.P_k)] = current_conditions.P_k * pressure_unit;
            recorded_v[current_index] = current_conditions.v_k * velocity_unit;
            recorded_t[current_index] = t * time_unit;
            current_index = current_index + 1; # Point to the next space in memory 

            # If we are in a multiple of rt, reset indexes
            if jjj == 2^iii
                jjj = 0;
                resetter = 1;
                indexes_to_save[current_to_save] = current_index - 1;
                current_to_save = current_to_save + 1;
            else
                number_of_extra_indexes = number_of_extra_indexes + 1;
            end

            if plotter == true
                # Plot current conditions
                tt = @sprintf(" t = %.5f (ms), CP = %g, vz = %.5f (mm/ms), centre = %.5f (mm)", 
                    (t*time_unit), contact_points, current_conditions.v_k * velocity_unit,
                    current_conditions.η_k[1] * length_unit
                );
                plot(circleX, circleY .+ recorded_z[current_index - 1],
                    xlabel = "r/R",
                    ylabel = "z/R",
                    xlims = [-plot_width*length_unit/N, plot_width*length_unit/N],
                    ylims = [-plot_width*length_unit/N + length_unit/2, plot_width*length_unit/N + length_unit/2],
                    linewidth = 2.5,
                    size = (860, 720),
                    location = (5, 5),
                    title = tt
                );
                η_k_plot = current_conditions.η_k;
                ηY = [η_k_plot[plot_width:-1:2]; η_k_plot[1:plot_width]] * length_unit;
                plot!(ηX, ηY, linedwidth = 3.5, color = :gray);

                gui();
            else
                @printf("t = %.10f (ms), CP = %g, vz = %.10f (mm/ms), centre = %.10f (mm) \n", 
                    t * time_unit, contact_points, current_conditions.v_k * velocity_unit, current_conditions.z_k * length_unit);
            end

            ####
            ####
            # ANALISE CONTACT TIME, MAXIMUM DEFLECTION, COEF OF RESTITUTION
            if (contact_flag == false && contact_points > 0) # if contact began, start counting
                contact_flag = true;
                maximum_deflection = recorded_z[current_index - 2, 1]; # Store position for later use, remember, it has dimensions!
                impact_velocity = round(recorded_v[current_index - 2], digits = 10); # Store initial impact velocity.
                zero_potential = recorded_η[current_index - 2, 1] + rS; # Store our zero-Potential 
                mechanical_energy_in = 1/2 * mS * ((recorded_v[current_index - 2])^2); # Mechanical energy at contact time

            elseif (contact_flag == true && (velocity_out_flag == false || lab_velocity_out_flag == false)) # Record last contact time
                if recorded_z[current_index - 1] > length_unit * (ZERO_HEIGHT + 1)
                    if lab_velocity_out_flag == false
                        lab_mechanical_energy_out = 1/2 * mS * ((recorded_v[current_index - 2])^2) +
                           mS * g * (recorded_z[current_index - 2] - zero_potential); # Lab Mechanical Energy Out;
                        lab_velocity_out_flag = true;
                    end
                elseif lab_velocity_out_flag == false
                    lab_contact_time = lab_contact_time + dt;
                end
                
                if contact_points == 0
                    if velocity_out_flag == false
                        mechanical_energy_out =   1/2 * mS * ((recorded_v[current_index - 2])^2) +
                        mS * g * (recorded_z[current_index - 2] - zero_potential); # Mechanical Energy Out;
                        velocity_out_flag = true;
                    end
                elseif velocity_out_flag == false
                    contact_time = contact_time + dt;
                end
            end

            if save_after_contact == false && recorded_u[current_index - 1] < 0 &&
                velocity_out_flag == true && lab_velocity_out_flag == true
                break; # end simulation if contact ended.
            end

            if (contact_flag == true && recorded_v[current_index - 1] >= 0) # Record maximum deflection
                if first_fall_flag == true # If velocity has changed sign, record maximum
                    maximum_deflection = maximum_deflection - recorded_z[current_index - 2];
                    first_fall_flag = false;
                    println("t = $(t*time_unit)");
                    println("vk = $(current_conditions.v_k * velocity_unit)");
                    println("contact_points = $contact_points");
                end
            end

            if recorded_v[current_index - 1] < 0 && first_fall_flag == false
                if velocity_out_flag     == false; contact_time     = Inf; end
                if lab_velocity_out_flag == false; lab_contact_time = Inf; end

                if save_after_contact == false; break; end
            end
        end
    end # end main while loop

    ###
    # Post processing
    ###

    # Round and dimentionalize contact time
    contact_time            = round(contact_time * time_unit, digits = 10);
    lab_contact_time        = round(lab_contact_time * time_unit, digits = 10);
    maximum_deflection      = round(maximum_deflection, digits = 10); # Round maximum deflection 
    maximum_time_step       = maximum_time_step * time_unit;
    coef_of_restitution     = mechanical_energy_out/mechanical_energy_in;
    lab_coef_of_restitution = lab_mechanical_energy_out/mechanical_energy_in;

    if export_data == true
        # Cut data
        indexes_to_save = indexes_to_save[1:(current_to_save - 1)];
        recorded_z = recorded_z[indexes_to_save];
        recorded_η = recorded_η[indexes_to_save, :];
        recorded_u = recorded_u[indexes_to_save, :];
        recorded_P = recorded_P[indexes_to_save];
        recorded_v = recorded_v[indexes_to_save];
        recorded_t = recorded_t[indexes_to_save, :];

        # TODO: Convert memory consuming data into single precision

        # Save data into a .jdl2 file
        exported_data_name = Dates.format(now(), "yyyymmddHHMMss");
        
        save("$SAVE_PATH/$exported_data_name.jld2", 
            "impact_velocity"       , impact_velocity,
            "contact_time"          , contact_time,
            "lab_contact_time"      , lab_contact_time, 
            "maximum_deflection"    , maximum_deflection,
            "recorded_η"            , recorded_η,
            "recorded_u"            , recorded_u,
            "recorded_z"            , recorded_z,
            "recorded_P"            , recorded_P,
            "recorded_t"            , recorded_t,
            "recorded_v"            , recorded_v,
            "Tm"                    , Tm,
            "mu"                    , μ,
            "rS"                    , rS,
            "recorded_time"         , recorded_time,
            "initial_velocity"      , initial_velocity, 
            "N"                     , N,
            "mS"                    , mS, 
            "coef_of_restitution"   , coef_of_restitution,
           "lab_coef_of_restitution", lab_coef_of_restitution
        );
        
        data_to_write = DataFrame(
            ID                   = exported_data_name,
            vi                   = impact_velocity,
            surfaceTension       = Tm,
            redius               = rS,
            cTime                = contact_time,
            maxDeflection        = maximum_deflection,
            coefOfRestitution    = coef_of_restitution,
            Density              = 3*mS/(4*pi*rS^3),
            labcTime             = lab_contact_time, 
            labcoefOfRestitution = lab_coef_of_restitution
        );
        CSV.write("$SAVE_PATH/$file_name", data_to_write, append = true);
      
    end
    @printf("Radius: %6.6g (mm)\n", rS);
    @printf("Tm: %0.5g \n", Tm);
    @printf("Contact time: %6.6g (ms)\n", contact_time);
    @printf("Contact time (according to lab procedures): %6.6g (ms)\n", lab_contact_time);
    @printf("Maximum deflection: %6.6g (mm)\n", maximum_deflection);
    @printf("Coefficient of Restitution squared: %6.6g \n", coef_of_restitution);
    @printf("Coefficient of Restitution squared (lab): %6.6g \n", lab_coef_of_restitution);
    @printf("dt: %6.6g (ms)\n", dt * time_unit);
    @printf("dr: %6.6g (mm)\n", dr * length_unit);
    println("-----------------------------");
    nothing
end # End main function declaration

function get_current_wd()::String


end

if (abspath(PROGRAM_FILE) == @__FILE__) || occursin(r"debugger"i, PROGRAM_FILE)
    solve_motion(plotter = true);
end