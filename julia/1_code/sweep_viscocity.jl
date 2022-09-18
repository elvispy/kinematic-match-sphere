include("solveMotion.jl");
using .kinematic_match
using CSV, DataFrames, Printf;

# Defining functions
function solveMotion_special(v_k::Float64, rho::Float64, vu::Float64, export_path::String)

    R_f = 105/2; # Membrane RADIUS (in mm)
    Tm = 107.0; # Linear tension of the material in N / m)
    rS = 4.76/2; # in mm
    mass(rS, rho) = rho * 4 * pi/3 * rS^3;
    μ = 0.3;

    solveMotion(v_k = v_k, rS = rS, Tm = Tm, R_f = R_f/rS, vu = vu, μ = μ, mS = mass(rS, rho), N = 100, plotter = false, 
                file_name = export_path, export_data = true, method = "EulerLinearized", save_after_contact = false);

end

# defining functions
function get_simulation_error(nb_of_simulations, export_file_path, experimental_data)
    values_num = DataFrame(CSV.File(export_file_path));; # Read CSV
    values_num = values_num[(end - nb_of_simulations + 1):end, :];

    vec_exp = map(vi -> findfirst(x -> x == vi, experimental_data.vi), values_num.vi);
    if nothing in vec_exp; throw("Something not found in experimental data!"); end

    total_error_tc    = 0; # Contact time
    total_error_delta = 0; # Maximum deflection
    total_error_alpha = 0; # Coefficient of restitution

    for ii = 1:size(values_num, 1)
        if ~(experimental_data.labcoefOfRestitution[vec_exp[ii]] - experimental_data.sigma_alpha[vec_exp[ii]] 
                <= values_num.labcoefOfRestitution[ii] && 
                    values_num.labcoefOfRestitution[ii] <= 
                experimental_data.labcoefOfRestitution[vec_exp[ii]] + experimental_data.sigma_alpha[vec_exp[ii]])
                total_error_alpha = total_error_alpha + (experimental_data.labcoefOfRestitution[vec_exp[ii]] -values_num.labcoefOfRestitution[ii])^2;
        end
        
        if ~(experimental_data.maxDeflection[vec_exp[ii]] - experimental_data.sigma_delta[vec_exp[ii]] 
                <= values_num.labcoefOfRestitution[ii] && 
                    values_num.labcoefOfRestitution[ii] <= 
                experimental_data.maxDeflection[vec_exp[ii]] + experimental_data.sigma_delta[vec_exp[ii]])
                total_error_delta = total_error_delta + (experimental_data.maxDeflection[vec_exp[ii]] -values_num.maxDeflection[ii])^2;
        end
        
        if ~(experimental_data.labcTime[vec_exp[ii]] - experimental_data.sigma_tc[vec_exp[ii]] ...
                <= values_num.labcoefOfRestitution[ii] && 
                    values_num.labcoefOfRestitution[ii] <= 
                experimental_data.labcTime[vec_exp[ii]] + experimental_data.sigma_tc[vec_exp[ii]])
                total_error_tc = total_error_tc + (experimental_data.labcTime[vec_exp[ii]] -values_num.labcTime[vec_exp[ii]])^2;
        end
        
        end
        
        final_error = total_error_alpha + total_error_tc + total_error_delta;
        errors = [total_error_alpha, total_error_tc, total_error_delta];
        
    return final_error, errors
end


# This script will try to predict the correct value for the elastic modulus of a given 
# Dataset

export_file_path = "../2_pipeline/sweep_viscocity/out/viscocity_sweep.csv";
cd(dirname(@__FILE__)); 
exp_data = DataFrame(CSV.File("../0_data/manual/experimental_data.csv")); # We need to read experimental data and filter 
rename!(exp_data, Symbol.(replace.(string.(names(exp_data)), Ref(r"\s"=>""))))
#print(names(exp_data));
subset!(exp_data, :Density => d -> d.== 3.25); # Filter densities
exp_data.labcoefOfRestitution = exp_data.labcoefOfRestitution.^2;
#predictor = last(exp_data, 1);

exp_data = exp_data[5:5, :];

#= Trying to debug
file_path = @__FILE__;
cd("../2_pipeline/$");
=#

baixa_viscocidade =  0.0;
alta_viscocidade  = 10.0;

for ii = 1:size(exp_data, 1)
    solveMotion_special(exp_data.vi[ii], 3.25, baixa_viscocidade, export_file_path);
end
erro_abaixo, _ = get_simulation_error(size(exp_data, 1), export_file_path, exp_data);

@printf("------------\n");
@printf("Resumo\n");
@printf("Viscocidade \n")
@printf("Viscocidade: %.5g\n", baixa_viscocidade);
@printf("Erro Total: %.5f\n", erro_abaixo);


for ii = 1:size(exp_data, 1)
    solveMotion_special(exp_data.vi[ii], 3.25, alta_viscocidade, export_file_path);
end
erro_acima, _ = get_simulation_error(size(exp_data, 1), export_file_path, exp_data);

@printf("------------\n");
@printf("Resumo\n");
@printf("Viscocidade \n")
@printf("Viscocidade: %.5g\n", alta_viscocidade);
@printf("Erro Total: %.5f\n", erro_acima);

errmin = Inf;
viscmin = NaN;

while alta_viscocidade - baixa_viscocidade > 0.05
    viscocidade_media = (baixa_viscocidade + alta_viscocidade)/2;
    for ii = 1:size(exp_data, 1)
        solveMotion_special(exp_data.vi[ii], 3.25, viscocidade_media, export_file_path);
    end
    erro_medio, t = get_simulation_error(size(exp_data, 1), export_file_path, exp_data);
    err_vec = [erro_medio, erro_abaixo, erro_acima];
    if erro_medio < errmin
        global errmin = erro_medio;
        global viscmin = viscocidade_media;
    end

    if erro_abaixo < erro_acima
        global alta_viscocidade = viscocidade_media;
    else
        global baixa_viscocidade = viscocidade_media;
    end

    @printf("---------------\n");
    @printf("Resumo\n");
    @printf("Viscocidade: %.5g\n", viscocidade_media);
    @printf("Erro Total: %.5f\n", erro_medio);
    @printf("Alpha: %.5f, TC: %.5f, Delta: %.5f\n", t[1], t[2], t[3]);
end

println("RESULTADO FINAL");
@printf("Menor erro: %.5f\n", errmin);
@printf("Atingido com vu=%.5f\n", viscmin);

@printf("Resultado final para rho = %.5g\n", 3.25);
@printf("Menor erro: %.5f\n", errmin);
@printf("Atingido com vu=%.5f\n", viscmin);
@printf("-------------------------\n-------------------------\n-------------------------\n-------------------------\n");


