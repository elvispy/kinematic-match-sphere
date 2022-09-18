% This script will sweep on viscocity to try to predict the correct value 
% of the viscocity of the membrane
export_file_path = "../2_pipeline/solveMotion2_1/out/viscocity_sweep.csv";
%n = 20;
%velo_samples = zeros(n, 1);
%coef_samples = zeros(n, 1);
exp_data = readtable("experimental_data.csv", 'PreserveVariableNames', true);
exp_data = exp_data(exp_data.Density == 3.25, :);
exp_data.labcoefOfRestitution = exp_data.labcoefOfRestitution.^2; % In simulations the coefficient of restitution is squared.
predictor = exp_data(end, :);
%exp_data = exp_data(1:(end-1),:); 
exp_data = exp_data(1:(2),:); 

file_path = fullfile("../2_pipeline/", mfilename, "out");
    if exist(file_path, 'dir') ~= 7 % CHeck if folder simulations exists
        mkdir(file_path); % If not, create it
    end
[fid, msg] = fopen(fullfile(file_path, "debug.txt"), "a+");

 % Intervalo de chute de viscocidade

baixa_viscocidade = 0;
alta_viscocidade  = 10;

for ii = 1:height(exp_data)
   solveMotion(exp_data.vi(ii), 3.25, baixa_viscocidade)  
end
erro_abaixo = get_simulation_error(height(exp_data), export_file_path, exp_data);

fprintf(fid, "---------------\n");
fprintf(fid, "Resumo\n");
fprintf(fid, "Viscocidade: %.5g\n", baixa_viscocidade);
fprintf(fid, "Erro Total: %.5f\n", erro_abaixo);
%fprintf(fid, "Alpha, TC, Delta");


for ii = 1:height(exp_data)
   solveMotion(exp_data.vi(ii), 3.25, alta_viscocidade) 
end
erro_acima  = get_simulation_error(height(exp_data), export_file_path, exp_data);

fprintf(fid, "---------------\n");
fprintf(fid, "Resumo\n");
fprintf(fid, "Viscocidade: %.5g\n", alta_viscocidade);
fprintf(fid, "Erro Total: %.5f\n", erro_acima);
%fprintf(fid, "Alpha, TC, Delta");


    
errmin = Inf;
viscmin = nan;
while alta_viscocidade - baixa_viscocidade > 0.05
    viscocidade_media = (baixa_viscocidade + alta_viscocidade)/2;
    for ii = 1:height(exp_data)
       solveMotion(exp_data.vi(ii), 3.25, viscocidade_media); 
    end
    [erro_medio, t]  = get_simulation_error(height(exp_data), export_file_path, exp_data);
    err_vec = [erro_medio, erro_abaixo, erro_acima];
    if erro_medio < errmin
        errmin = erro_medio;
        viscmin = viscocidade_media;
    end
    
    if erro_abaixo < erro_acima
        alta_viscocidade = viscocidade_media;
    else
        baixa_viscocidade = viscocidade_media;
    end
    
    
    fprintf(fid, "---------------\n");
    fprintf(fid, "Resumo\n");
    fprintf(fid, "Viscocidade: %.5g\n", viscocidade_media);
    fprintf(fid, "Erro Total: %.5f\n", erro_medio);
    fprintf(fid, "Alpha: %.5f, TC: %.5f, Delta: %.5f\n", t(1), t(2), t(3));

end

disp("RESULTADO FINAL");
fprintf("Menor erro: %.5f\n", errmin);
fprintf("Atingido com vu=%.5f\n", viscmin);



fprintf(fid, "Resultado final para rho = %.5g\n", 3.25);
fprintf(fid, "Menor erro: %.5f\n", errmin);
fprintf(fid, "Atingido com vu=%.5f\n", viscmin);
fprintf(fid, "-------------------------\n-------------------------\n-------------------------\n-------------------------\n");


function solveMotion(v_k, rho, vu)
    arguments
       v_k (1, 1) {mustBeReal}     = -0.1;  % Initial velocity, in mm/ms
       rho (1, 1) {mustBePositive} = 3.25; % Sphere's density, in g/cm3 = mg/mm^3
       vu  (1, 1) {mustBeNonnegative} = 0; % Membrane viscocity, in mg/ms
    end
    R_f = 105/2; % Membrane RADIUS (in mm)
    Tm = 107;  %Linear tension of the material (in N / m = kg*m/(s2*m) = mg / ms^2)
    rS = 4.76/2; %(both in mm)
    mass = @(rS, rho) rho * (4*pi/3) * rS^3; % in mg
    %rho = 3.25; 
    mu = 0.3; %Density of membrane per unit of area (mg/mm^2) (Sara Wrap membrane)
    
    solveMotion2_1(...
        'rS', rS, ...
        'Tm', Tm, ...
        'R_f', R_f/rS, ...
        'vu', vu, ...
        'mu', mu, ...
        'mS', mass(rS, rho), ...
        'v_k'     , -abs(v_k), ...
        'N'       , 25, ...
        'plotter' , false, ...
        'FileName', 'viscocity_sweep.csv', ...
        'export_data', true, ...
        'method', 'EulerLinearized', ...
        'save_after_contact_ended', false ...
        ); 
end

function [final_error, errors] = get_simulation_error(nb_of_simulations, export_file_path, experimental_data)
    % This function will return the error that is related to the last
    % simulations 
    valuesNum = readtable(export_file_path, 'PreserveVariableNames', true);
    %valuesExp = valuesExp(valuesExp.Viscocity == vu , :);
    %valuesExp = valuesExp(valuesExp.Density   == rho, :);
    valuesNum = valuesNum((end-nb_of_simulations+1):end, :);
    %valuesExp = sortrows(valuesExp, 'vi');
    %experimental_data = sortrows(experimental_data, 'vi');
    
    [loca, locb] = ismember(valuesNum{:, 'vi'}, experimental_data{:, 'vi'});
    if ~all(loca); error("Something not found!"); end
    total_error_tc    = 0; % Contact time
    total_error_delta = 0; % Maximum deflection
    total_error_alpha = 0; %Coefficient of restitution
    for ii = 1:height(valuesNum)
       if ~(experimental_data{locb(ii), 'labcoefOfRestitution'} - experimental_data{locb(ii), 'sigma_alpha'} ...
               <= valuesNum{ii, 'labcoefOfRestitution'} && ...
                  valuesNum{ii, 'labcoefOfRestitution'} <= ...
            experimental_data{locb(ii), 'labcoefOfRestitution'} + experimental_data{locb(ii), 'sigma_alpha'})
            total_error_alpha = total_error_alpha + (experimental_data{locb(ii), 'labcoefOfRestitution'} -valuesNum{ii, 'labcoefOfRestitution'})^2;
       end
       
       if ~(experimental_data{locb(ii), 'maxDeflection'} - experimental_data{locb(ii), 'sigma_delta'} ...
               <= valuesNum{ii, 'labcoefOfRestitution'} && ...
                  valuesNum{ii, 'labcoefOfRestitution'} <= ...
            experimental_data{locb(ii), 'maxDeflection'} + experimental_data{locb(ii), 'sigma_delta'})
            total_error_delta = total_error_delta + (experimental_data{locb(ii), 'maxDeflection'} -valuesNum{ii, 'maxDeflection'})^2;
       end
       
       if ~(experimental_data{locb(ii), 'labcTime'} - experimental_data{locb(ii), 'sigma_tc'} ...
               <= valuesNum{ii, 'labcoefOfRestitution'} && ...
                  valuesNum{ii, 'labcoefOfRestitution'} <= ...
            experimental_data{locb(ii), 'labcTime'} + experimental_data{locb(ii), 'sigma_tc'})
            total_error_tc = total_error_tc + (experimental_data{locb(ii), 'labcTime'} -valuesNum{locb(ii), 'labcTime'})^2;
       end
       
    end
    
    final_error = total_error_alpha + total_error_tc + total_error_delta;
    errors = [total_error_alpha, total_error_tc, total_error_delta];
    
    
end
