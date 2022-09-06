% Author: Elvis Agüero
% email: elvisavfc65@gmail.com
% Date: March 31th, 2022.

function solveMotion2_1(options)
    %solveMotion_2.1
    %Tries to solve the kinematic match between a rigid spherical object
    %and an axissymmetric membrane in vacuum conditions, Newton's Method and 
    % BDF2 finite differences.
    
    %% Handling default values
    arguments
        % All values are assumed to be given with dimensions.
       
        %Radius of the sphere (in mm)
        options.rS (1, 1) {mustBePositive} = 2.38
        %Tension of the material (in mg / ms^2)
        options.Tm (1, 1) {mustBePositive} = 107
        %Number of raddi in half a length of the membrane (dimensionless)
        options.R_f (1, 1) {mustBePositive} = inf
        %Density of membrane per unit of area (mg/mm^2) (Sara Wrap membrane)
        options.mu (1, 1) {mustBePositive} = 0.3 
        %Mass of the ball (in mg) (3.25 is the ball's density in mg/mm3)
        options.mS (1, 1) {mustBePositive} = inf
        %Gravity of earth (in mm/ms^2)
        options.g (1, 1) {mustBePositive} = 9.80665e-3 %(mm/ms2);
        
        options.z_k   (1, 1) {mustBeReal}    = inf % in mm
        options.v_k   (1, 1) {mustBeReal}    = -0.6312 % in mm/ms, as in the paper
        options.Eta_k (1, :) {mustBeNumeric} = []
        options.u_k   (1, :) {mustBeNumeric} = []
        options.rhoa  (1, 1) {mustBeReal}    = 1.225e-3 % FOR NOW ADIOMENSIONAL
        options.P_k   (1, :) {mustBeNumeric} = []
        
        options.SimulTime             (1, 1) {mustBePositive} = 15 % in ms
        options.RecordedTime          (1, 1)                  = 0 % in ms
        options.N                     (1, 1) {mustBeInteger}  = 100 % Adimensional
        options.FileName              (1, 1) string           = "historial.csv"
        options.plotter               (1, 1) logical          = false
        options.method                (1, 1) string {mustBeMember(options.method, ...
        ["EulerLinearized", "EulerFullCurvature", ...
        "BDF2Linearized", "BDF2FullCurvature"])}              = "EulerLinearized"
        options.save_after_contact_ended (1, 1) logical          = false
        options.exportData            (1, 1) logical          = true
    end
    
    if options.save_after_contact_ended == false
        options.SimulTime = 20; 
    end
    
     mu = options.mu; rS = options.rS; Tm = options.Tm; g = options.g;
    if isinf(options.R_f); R_f = 52.5/rS; else; R_f = options.R_f; end
    if isinf(options.mS); mS = 3.25 * 4 * pi * (rS.^3) / 3; else; mS = options.mS; end
    %% CONSTANTS' DEFINITIONS
    
    Fr = (mu*rS*g)/Tm; %Froude number
    Mr = mu*(rS^2)/mS; %Ratio of masses1
    D = options.rhoa * rS / mu;

    %Units (Just for the record)
    Lunit = rS; %Spatial unit of mesurement (in mm)
    Vunit = sqrt(Tm/mu); %Velocity in mm/ms
    Tunit = Lunit/Vunit; %TempSoral unit of measurement (in ms)
    Punit = mu * Lunit / Tunit^2; %Unit of pressure (in mg/(ms^2 * mm))

    %Numerical constants
    
    N = options.N; %Number of dx intervals in the spatial (radial) coordinate per unit length
    Ntot = ceil(R_f * N); %Total number of non-trivial points in the spatial coordinate
    dr = 1/N; % %Spatial step
    
    %% Initial Conditions

    P_k = options.P_k/Punit; %Current vector of pressures. Only non trivial points saved (Empty vector == membrane and ball should not be touching)
    cPoints = length(P_k); %Initial number of contact points
    
    % Set Eta_k initial contidion
    if length(options.Eta_k) ~= Ntot
        if contains(options.method, "curvature", 'IgnoreCase', true)
            Eta_k = initial_condition_full_curvature(Fr, dr, Ntot); 
        else
            Eta_k = initial_condition_linearized(Fr, dr, Ntot); 
        end
    else
        if size(options.Eta_k, 2) > 1; Eta_k = options.Eta_k'; else; Eta_k = options.Eta_k; end
        Eta_k = Eta_k/Lunit;
    end
    % Set u_k initial condition
    if Ntot == length(options.u_k)
        if size(options.u_k, 2) > 1; u_k = options.u_k'; else; u_k= options.u_k; end
        u_k = u_k/Vunit;
    else
        u_k = zeros(Ntot, 1); %Initial vertical velocity of the membrane (a.k.a. d/dt Eta_k = u_k)    
    end

    % set sphere's height initial condition
    if options.z_k == inf
        z_k = Eta_k(1) + 1; %Curr1ent position of the center of the ball (Dimensionless)
    else
        assert(options.z_k >= Eta_k(1) + 1, ...
            "Sphere is not in a valid position!");
        z_k = options.z_k/Lunit;
    end

    % set sphere's velocity initial condition
    v_k = options.v_k/Vunit; 
    
    % Zero height to measure contact time experimentally
    ZERO_HEIGHT = Eta_k(1);    
    
    probableNextConditions = cell(5, 1);
    if options.RecordedTime > 0; dt = options.RecordedTime/Tunit; end
    
    if isscalar(regexp(options.method, "bdf.*full", 'ignorecase'))
        if options.RecordedTime == 0; dt = 10/(N^2 * Tunit); end
        returnMatrices2_1 = @returnFullMatrices;
        getNextStep = @BDF2FullCurvature;
        previousConditions = struct('Eta_k', Eta_k, ...
            'u_k', u_k, 'z_k', z_k, 'v_k', v_k, 'P_k', P_k, ...
            'dt', dt);
        currentConditions = previousConditions; % TODO: Uncomment next line
        %[currentConditions, previousConditions.dt] = advanceOneStep(previousConditions, options.method);
    elseif isscalar(regexp(options.method, "bdf.*linear", 'ignorecase'))
        if options.RecordedTime == 0; dt = 1/(N * Tunit); end
        returnMatrices2_1 = @returnLinearizedMatrices;
        getNextStep = @BDF2LinearizedCurvature;
        previousConditions = struct('Eta_k', Eta_k, ...
            'u_k', u_k, 'z_k', z_k, 'v_k', v_k, 'P_k', P_k, ...
            'dt', dt);
        currentConditions = previousConditions; % TODO: Uncomment next line
        %[currentConditions, previousConditions.dt] = advanceOneStep(previousConditions, options.method);
    elseif isscalar(regexp(options.method, "euler.*linear", 'ignorecase'))
        if options.RecordedTime == 0; dt = 1/(N * Tunit); end
        
        returnMatrices2_1 = @returnLinearizedMatrices;
        getNextStep = @eulerLinearized;
        currentConditions = struct('Eta_k', Eta_k, ...
            'u_k', u_k, 'z_k', z_k, 'v_k', v_k, 'P_k', P_k, ...
            'dt', dt);
        previousConditions = struct('Eta_k', NaN, ...
            'u_k', NaN, 'z_k', NaN, 'v_k', NaN, 'P_k', NaN, ...
            'dt', NaN);
    elseif isscalar(regexp(options.method, "euler.*full", 'ignorecase'))
        if options.RecordedTime == 0; dt = 10/(N^2 * Tunit); end
        returnMatrices2_1 = @returnFullMatrices;
        getNextStep = @eulerNewton;
        currentConditions = struct('Eta_k', Eta_k, ...
            'u_k', u_k, 'z_k', z_k, 'v_k', v_k, 'P_k', P_k, ...
            'dt', dt);
        previousConditions = struct('Eta_k', NaN, ...
            'u_k', NaN, 'z_k', NaN, 'v_k', NaN, 'P_k', NaN, ...
            'dt', NaN);
    else
        error("Not a valid or implemented method: '%s'", options.method);
    end


    final_time = options.SimulTime/Tunit; %Maximum time (dimensionless) to be simulated.
    t_init = 0/Tunit; %Initial time (dimensionless)
    rt = dt; %Interval of time to be recorded in the matrix. (1/0.25 = 4 frames per ms)
    t = t_init; % t is the variable that will keep track of the time
    current_index = 2; %for saving matrix
    
    %For preallocation efficiency
    maximum_index = ceil((final_time-t_init)/rt) + 4;
    number_of_extra_indexes = 0;
    
    %% Flow control Variables
    
    cTime = 0; labcTime = 0;%To record contact time
    maxDef = 0; %to record max deflection;
    firstFallFlag = true; %Boolean to record maximum deflection time
    contactFlag = false; %To record first contact time
    velocityOutRecorded = false; labvelocityOutRecorded = false;  % To check if Em_out was recorded
    iii = 0; jjj = 0; % Indexes to keep track  of subdivisions of rt
    resetter = 0;
    max_nb_cPoints = N + 1; %Maximum number of allowed contact points
    plotter = options.plotter;
    %getNextStep = options.method;
    file_path = fullfile("../2_pipeline/", mfilename, "out");
    if exist(file_path, 'dir') ~= 7 % CHeck if folder simulations exists
        mkdir(file_path); % If not, create it
    end
    file_name = fullfile(file_path, options.FileName);
    if exist(file_name, 'file') == 0 %Check if .csv file exists
        writematrix(["ID", "vi", "surfaceTension", "radius", ...
            "cTime", "maxDeflection", "coefOfRestitution", "Density", ...
            "labcTime", "labcoefOfRestitution"], file_name, ...
            'WriteMode', 'append'); % If not, create it
    end
    
    %Now we save the sphere in an array (For plotting)
    width = min(Ntot, ceil(3 * N));
    xplot = linspace(0, width/N, width);
    EtaX = [-fliplr(xplot(2:end)), xplot] * Lunit;
    %EtaU = zeros(1, 2*width - 1);
    %step = ceil(N/15);
    if options.plotter == true
        close all;
        myFigure = figure('Position', [10 10 1080 960]);
        axis equal;
        set(gca, ...
            'Box'         , 'on'     , ...
            'FontSize'    , 16        , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.02 .02] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'YTick'       , -ceil(2*width*Lunit/N):1:ceil(2*width*Lunit/N), ...
            'XTick'       , -ceil(2*width*Lunit/N):1:ceil(2*width*Lunit/N));
        xlabel('$r/R$', 'Interpreter', 'latex');
        ylabel('$\frac{z}{R}\ \ $', 'Interpreter', 'latex', 'Rotation', 0, ...
            'FontSize', 24);
        
        xlim([-width*Lunit/N, width*Lunit/N]);
        ylim([-width*Lunit/N + Lunit/2, width*Lunit/N + Lunit/2]);
        %ylim([(z_k-2)*Lunit, (z_k+2)*Lunit]);
    end
    
    %% For post processing

    % Store recorded values with their dimensions.
    recorded_z    = zeros(maximum_index, 1);    recorded_z(1) = z_k * Lunit;
    recorded_Eta  = zeros(maximum_index, Ntot); recorded_Eta(1, :) = Eta_k * Lunit;
    recorded_u    = zeros(maximum_index, Ntot); recorded_u(1, :) = u_k * Vunit;
    recorded_P    = zeros(maximum_index, max_nb_cPoints); %recordedPk(1, :) = Pk;
    recorded_v    = zeros(maximum_index, 1);    recorded_v(1) = v_k * Vunit;
    recordedTimes = zeros(maximum_index, 1);    recordedTimes(1) =  t_init * Tunit;
    vars = {datestr(now, 'dddHHMMss'), 0, Tm, rS};
    
    % For coefficient of restitution
    Em_in = nan;
    Em_out = nan; labEm_out = nan;

    jacobian_pieces = cell(max_nb_cPoints + 1, 1);
    jacobian_pieces = returnMatrices2_1(jacobian_pieces, Ntot, dr, Mr);
    
    indexes_to_save = zeros(maximum_index, 1); indexes_to_save(1) = 1;
    current_to_save = 2;
    
     %% Main Loop
    while (t <= final_time)
        errortan = Inf * ones(1, 5);
        recalculate = false;

        % If empty, allocate more space for matrices.
        if isempty(jacobian_pieces{cPoints + 3}) == true
            jacobian_pieces = returnMatrices2_1(jacobian_pieces, Ntot, dr, Mr);
        end
        %First, we try to solve with the same number of contact points
        [probableNextConditions{3}, errortan(3)] = getNextStep(cPoints, max_nb_cPoints, ...
            currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, D);

        if abs(errortan(3)) < 1e-8 %If almost no error, we accept the solution
            previousConditions = currentConditions;
            currentConditions  = probableNextConditions{3};

        else % If there is some error, we try with diffferent contact points
            %Lets try with one more point
            [probableNextConditions{4}, errortan(4)] = getNextStep(cPoints + 1, max_nb_cPoints, ...
                currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, D);
            %Lets try with one point less
            [probableNextConditions{2}, errortan(2)] = getNextStep(cPoints - 1, max_nb_cPoints, ...
                currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, D);

            if (abs(errortan(3)) > abs(errortan(4)) || abs(errortan(3)) > abs(errortan(2)))
                if abs(errortan(4)) <= abs(errortan(2))
                    %Now lets check with one more point to be sure
                    [~, errortan(5)] = getNextStep(cPoints + 2, max_nb_cPoints, ...
               currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, D);

                    if abs(errortan(4)) < abs(errortan(5))
                        %Accept new data
                        previousConditions = currentConditions;
                        currentConditions = probableNextConditions{4};
                        cPoints = cPoints + 1;
                    else
                        %time step too big
                        recalculate = true;
                    end
                else
                    %now lets check if errortan is good enough with one point
                    %less
                    [~, errortan(1)] = getNextStep(cPoints-2, max_nb_cPoints, ...
                currentConditions, previousConditions, dt, dr, Fr, Ntot, jacobian_pieces, D);

                    if abs(errortan(2)) < abs(errortan(1))
                        %Accept new data
                        previousConditions = currentConditions;
                        currentConditions = probableNextConditions{2};
                        cPoints = cPoints - 1;
                    else
                        recalculate = true;
                    end
                end %End of (errortan(4) < errortan(2))

            else %the same number of contact points is best    
                if errortan(3) == Inf % ALl errors are infinity
                    recalculate = true;
                else
                    %Accept new data
                    previousConditions = currentConditions;
                    currentConditions = probableNextConditions{3};
                end
            end %
        end
        
        if recalculate == true
            dt = dt/2;
            % Refine time step in index notation
            iii = iii + 1; jjj = 2 * jjj;
        else
            t = t + dt; jjj = jjj + 1;
            if mod(jjj, 2) == 0 && resetter > 0
                jjj = jjj/2;
                iii = iii - 1;
                % Increase time step
                dt = 2 * dt;
                %Decrease the number of times you can reset the time
                resetter = resetter - 1;
            end
            
            %Update Indexes if necessary and store results into a matrix
            if current_index > maximum_index
                recorded_z   = [recorded_z;   zeros(number_of_extra_indexes, 1)];
                recorded_Eta   = [recorded_Eta;   zeros(number_of_extra_indexes, Ntot)];
                recorded_u   = [recorded_u;   zeros(number_of_extra_indexes, Ntot)];
                recorded_P    = [recorded_P;    zeros(number_of_extra_indexes, max_nb_cPoints)];
                recorded_v   = [recorded_v;   zeros(number_of_extra_indexes, 1)];
                recordedTimes = [recordedTimes; zeros(number_of_extra_indexes, 1)];
                maximum_index = maximum_index + number_of_extra_indexes;
                number_of_extra_indexes = 0;
            end
            
            % Store Data
            recorded_z(current_index) = currentConditions.z_k * Lunit;
            recorded_Eta(current_index, :) = currentConditions.Eta_k * Lunit;
            recorded_u(current_index, :) = currentConditions.u_k * Vunit;
            recorded_P(current_index, 1:length(currentConditions.P_k)) = currentConditions.P_k * Punit;
            recorded_v(current_index) = currentConditions.v_k * Vunit;
            recordedTimes(current_index) = t * Tunit;
            current_index = current_index + 1; % Set the current index to +1
            
            % If we are in a multiple of rt, reset indexes
            if jjj == 2^iii %Alternatively, 2*mod(t+dt/16, rt) < dt
                %disp(errortan);
                jjj = 0; resetter = 1;
                indexes_to_save(current_to_save) = current_index - 1;
                current_to_save = current_to_save + 1;
            else
                number_of_extra_indexes = number_of_extra_indexes + 1;
            end
            
            %%%%%%%%%%%% Plotting current conditions
            if plotter == true
                %PLOT 
                %figure(myFigure);
                cla; hold on; grid on;
                %plot(circleX*Lunit,(z_k+circleY)*Lunit,'k','Linewidth',1.5);
                rectangle('Position', [-Lunit, (currentConditions.z_k-1)*Lunit, 2*Lunit, 2*Lunit], ...
                    'Curvature', 1, 'Linewidth', 2);
                Eta_k = currentConditions.Eta_k;
                plot([-width*Lunit/N, width*Lunit/N], [0, 0], '--', 'Color', [.8 .8 .8], 'LineWidth', 2.5);
                EtaY = [flipud(Eta_k(2:width));Eta_k(1:width)]' * Lunit;
                %EtaV = [flipud(u_k(2:width));u_k(1:width)]' * Vunit;
                plot(EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
                
                plot(0, (currentConditions.z_k-1)*Lunit, 'Marker', 'o', 'MarkerSize', 12, ...
                    'LineWidth', 1, 'Color', [0 0 0]);
                plot(0, Eta_k(1)*Lunit, 'Marker', 'o', 'MarkerSize', 4, ...
                    'MarkerFaceColor', [.3 .3 .3], 'LineWidth', 1.5, 'Color', [.3 .3 .3]);
                if contactFlag == false
                    %quiver(EtaX(1:step:end), EtaY(1:step:end), EtaU(1:step:end), EtaV(1:step:end), 0);
                else
                    %quiver(EtaX(1:step:end), EtaY(1:step:end), EtaU(1:step:end), EtaV(1:step:end));
                end
               title(sprintf(" t = %-8.5f, CP = %-8g, \n v_k = %-8.5f, z_k = %-8.5f \n", t*Tunit, cPoints, ...
                        currentConditions.v_k*Vunit, currentConditions.z_k * Lunit),'FontSize',16);
                drawnow;
            else
                %clc;
                fprintf("t = %.10f, CP = %g, v_k = %.10f, z_k = %.10f \n", t*Tunit, cPoints, ...
                        currentConditions.v_k*Vunit, currentConditions.z_k * Lunit);
            end
            
            %%%%%%%%%%%%
            %%%%%%%%%%%%
            % ANALISE CONTACT TIME, MAXIMUM DEFLECTION, COEF OF RESTITUTION
            if (contactFlag == false && cPoints > 0) % if contact began, start counting
                contactFlag = true;
                maxDef = recorded_z(current_index - 2, 1); %Store position for later use, remember, it has dimensions!
                vars{2} = round(recorded_v(current_index - 2, 1), 10); % Store initial impact velocity
                zeroPotential = recorded_Eta(current_index - 2, 1) + rS; % Store our zero-Potential 
                Em_in = 1/2 * mS * ((recorded_v(current_index - 2))^2); % Mechanical Energy In;
                if options.plotter == true
                    %tt = text(-2.4, 2.7, '$a)$', 'Interpreter', 'latex', 'FontSize', 16);
                    %print(myFigure, '-depsc', '-r300', 'Graficos/ImminentContact.eps');
                    %delete(tt);
                end
            elseif (contactFlag == true && ...
                    (velocityOutRecorded == false || labvelocityOutRecorded == false)) % record last contact time
                if recorded_z(current_index - 1) > Lunit * (ZERO_HEIGHT + 1)
                    if labvelocityOutRecorded == false
                        labEm_out = 1/2 * mS * ((recorded_v(current_index - 2))^2) ...
                            + mS*g*(recorded_z(current_index - 2) - zeroPotential); % Mechanical Energy Out;
                        labvelocityOutRecorded = true;
                    end
                elseif labvelocityOutRecorded == false
                    labcTime = labcTime + dt;
                end
                if cPoints == 0 
                    if velocityOutRecorded == false
                        Em_out = 1/2 * mS * ((recorded_v(current_index - 2))^2) ...
                            + mS*g*(recorded_z(current_index - 2) - zeroPotential); % Mechanical Energy Out;
                        velocityOutRecorded = true;
                    end
                elseif velocityOutRecorded == false
                    cTime = cTime + dt;
                end
            end
            if options.save_after_contact_ended == false && recorded_u(current_index - 1, 1) < 0 &&...
                    velocityOutRecorded == true && labvelocityOutRecorded == true
                break; % end simulation if contact ended.
            end
            if (contactFlag == true && recorded_v(current_index - 1) >= 0) % Record maximum deflection
                if (firstFallFlag == true) %If velocity has changed sign, record maximum 
                    %deflection only once
                    maxDef = maxDef - recorded_z(current_index - 2);
                    firstFallFlag = false;
                    if options.plotter == true
%                         tt = text(-2.4, 2.7, '$b)$', 'Interpreter', 'latex', 'FontSize', 16);
%                         print(myFigure, '-depsc', '-r300', 'Graficos/MaximumDeflection.eps');
%                         delete(tt);
                    end
                end
            end
            if  recorded_v(current_index - 1) < 0 && firstFallFlag == false 
                if velocityOutRecorded    == false; cTime    = inf; end
                if labvelocityOutRecorded == false; labcTime = inf; end
                
                if options.save_after_contact_ended == false
                    break;
                end
            end
        end
    end % END WHILE
    
    
    %% POST PROCESSING
    %%%%%%%%%%%%
    if options.plotter == true
        delete(myFigure);
    end
    
    % Round and dimentionalize contact time
    cTime = round(cTime * Tunit, 10);
    labcTime = round(labcTime * Tunit, 10);
    maxDef = round(maxDef, 10); % Round maximum deflection
    v_k = vars{2}; % Initial velocity is already rounded
    rt = rt * Tunit;
    coefOfRestitution = Em_out/Em_in;
    labcoefOfRestitution = labEm_out/Em_in;
    
    if options.exportData == true 
        %%%
        %%% TODO LO QUE SE EXPORTA ES DIMENSIONAL
        %%%
        
        %Cut data
        indexes_to_save = indexes_to_save(1:(current_to_save-1));
        indexes_to_save = indexes_to_save(1:3:end); %% BEWARE
        recorded_z = recorded_z(indexes_to_save);
        recorded_Eta = recorded_Eta(indexes_to_save, :);
        recorded_u = recorded_u(indexes_to_save, :);
        recorded_P = recorded_P(indexes_to_save, :);
        recorded_v = recorded_v(indexes_to_save, :);
        recordedTimes = recordedTimes(indexes_to_save, :);
        
        % Convert memory-consuming data into single precision
        recorded_Eta = single(recorded_Eta);
        recorded_u = single(recorded_u);
        
        % Save into a .mat file
        save(fullfile(file_path, sprintf("%s.mat",vars{1})), ...
            'recorded_Eta', 'recorded_u', 'recorded_z', 'recorded_P', 'recordedTimes', ...
            'recorded_v', 'cTime', 'labcTime', 'maxDef', 'Tm', 'mu', 'rS', 'rt', 'v_k', 'N', 'mS', ...
            'coefOfRestitution', 'labcoefOfRestitution');
        writecell([vars, cTime, maxDef, coefOfRestitution, 3*mS/(4*pi*rS^3), ...
            labcTime, labcoefOfRestitution], file_name, 'WriteMode', 'append');
    end
    
    % Plot summary into matlab console.
    fprintf('Radius: %6.6g (mm)\n', rS);
    fprintf('Tm: %0.5g \n', Tm);
    fprintf('Contact time: %6.6g (ms)\n', cTime);
    fprintf('Contact time (according to lab procedures): %6.6g (ms)\n', labcTime);
    fprintf('Maximum deflection: %6.6g (mm)\n', maxDef);
    fprintf('Coefficient of Restitution squared: %6.6g \n', coefOfRestitution);
    fprintf('Coefficient of Restitution squared (lab): %6.6g \n', labcoefOfRestitution);
    fprintf('dt: %6.6g (ms)\n', dt*Tunit);
    fprintf('dr: %6.6g (mm)\n', dr * Lunit);
    disp('-----------------------------');
    
end % end function
    